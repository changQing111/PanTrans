#!/usr/bin/env python3
"""Prepare, submit and account for immutable PanTrans experiments (Linux/Slurm)."""
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import uuid

BASE = Path('/data/changq/Project/pantrans_test')
PYTHON = '/data/changq/miniforge3/bin/python'
TOOLS = '/home/changq/my_data/bio_soft/cactus-bin-v2.9.8/bin'
TD = '/home/changq/my_data/bio_soft/TransDecoder-TransDecoder-v5.7.1'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save(path, data):
    path = Path(path)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    os.replace(f.name, path)


def capture(argv, cwd=None):
    p = subprocess.run(argv, cwd=cwd, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    return {'exit_code': p.returncode, 'output': p.stdout}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def identity(path):
    p = Path(path).resolve()
    s = p.stat()
    if not p.is_file() or not s.st_size:
        raise ValueError(f'Missing/empty input: {p}')
    return {'path': str(p), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}


def update_index(root):
    root = Path(root)
    columns = ['run_id', 'preset', 'status', 'pantrans_version', 'git_commit',
               'history_run_id', 'job_id', 'elapsed_seconds', 'gnu_time_max_rss_kib',
               'slurm_elapsed_seconds', 'slurm_step_max_rss_kib',
               'output_dir', 'purpose']
    with (root / '.index.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = []
        for p in sorted(root.glob('*/run.json')):
            d = json.loads(p.read_text())
            rows.append({k: d.get(k, '') for k in columns})
        with tempfile.NamedTemporaryFile(mode='w', dir=root, delete=False) as f:
            w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        os.replace(f.name, root / 'runs.tsv')


def persist(run, d):
    save(run / 'run.json', d)
    update_index(run.parent)


def prepare(a):
    root = a.runs.resolve(); root.mkdir(parents=True, exist_ok=True)
    run = root / (dt.datetime.now().strftime('%Y%m%dT%H%M%S') + '_' + a.preset + '_' + uuid.uuid4().hex[:8])
    base = a.base.resolve(); repo = a.repo.resolve()
    prefix = 'five_variety' if a.preset == 'construct5' else 'ten_variety'
    inputs = {}; args = ['append' if a.preset == 'append5plus5' else 'construct']
    def add(flag, path):
        inputs[flag] = identity(path); args.extend([flag, inputs[flag]['path']])
    history_id = None
    if a.preset == 'append5plus5':
        history = (a.history_dir or base / 'five_variety_construct_out8').resolve()
        if a.history_graph_cdna:
            add('--history-graph-cdna', a.history_graph_cdna)
        if (history.parent / 'run.json').exists():
            history_id = json.loads((history.parent / 'run.json').read_text())['run_id']
        for flag, name in [('--cdna', 'five_variety_unrenamed_cdna.refer.fasta'),
                           ('--history-gff3', 'five_variety_unrenamed.gff3'),
                           ('--history-graph', 'five_variety.graph.json')]:
            add(flag, history / name)
        for flag, name in [('--query-cdna', 'five_variety_append_mrna.fa'),
                           ('--gdna', 'five_variety_append_gene3kb_gdna.fa'),
                           ('--bed', 'ten_variety_append_genes.bed')]:
            add(flag, base / 'five_append' / name)
        # Track all files the graph loader implicitly opens, not just CLI inputs.
        g = json.loads((history / 'five_variety.graph.json').read_text())
        for key in ['bed', 'cdna', 'gdna', 'nodes', 'edges', 'filtered_bam']:
            entries = g[key] if isinstance(g[key], list) else [g[key]]
            for i, entry in enumerate(entries):
                p = Path(entry['path'])
                inputs[f'history:{key}:{i}'] = identity(p if p.is_absolute() else history / p)
    else:
        folder = base if a.preset == 'construct5' else base / 'ten_variety'
        for flag, suffix in [('--cdna', '_mrna.fa'), ('--gdna', '_gene3kb_gdna.fa'), ('--bed', '_genes.bed')]:
            add(flag, folder / (prefix + suffix))
        add('--main-chroms', folder / 'main_chroms')
        args += ['--reference', 'CS']
        if a.preset == 'construct5' and not a.full:
            add('--bam', base / 'five_variety_construct_out2/five_variety_cdna_align_gdna.bam')
    rescue = a.rescue if a.rescue is not None else a.preset != 'construct10'
    if rescue:
        args.append('--rescue')
    args += ['--threads', str(a.cpus), '--prefix', prefix, '--annotate-cds',
             '--transdecoder-parallel', '--transdecoder-dir', TD, '--output', str(run / 'results')]
    run.mkdir(); (run / 'logs').mkdir(); (run / 'snapshot').mkdir()
    shutil.copytree(repo / 'src', run / 'snapshot/src', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo / 'pyproject.toml', run / 'snapshot/pyproject.toml')
    shutil.copy2(__file__, run / 'runner.py')
    snapshot_files = {str(p.relative_to(run / 'snapshot')): digest(p)
                      for p in sorted((run / 'snapshot').rglob('*')) if p.is_file()}
    save(run / 'snapshot.sha256.json', snapshot_files)
    status = capture(['git', 'status', '--porcelain'], repo)
    (run / 'source.diff').write_text(capture(['git', 'diff', 'HEAD', '--binary'], repo)['output'])
    (run / 'source.status.txt').write_text(status['output'])
    env = dict(os.environ, PYTHONPATH=str(run / 'snapshot/src'), PYTHONDONTWRITEBYTECODE='1')
    code = 'import json,sys; from pantrans.main import read_parameters; print(json.dumps(vars(read_parameters(sys.argv[1:]))))'
    parsed = subprocess.run([PYTHON, '-c', code] + args, env=env, cwd=run,
                            capture_output=True, text=True, check=True)
    version = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)',
                        (run / 'snapshot/src/pantrans/version.py').read_text()).group(1)
    d = dict(schema_version=1, run_id=run.name, preset=a.preset, purpose=a.purpose,
             status='prepared', created_at=now(), pantrans_version=version,
             git_commit=capture(['git', 'rev-parse', 'HEAD'], repo)['output'].strip(),
             dirty=bool(status['output'].strip()), source_repo=str(repo),
             snapshot_files=snapshot_files, inputs=inputs, argv=args,
             resolved_parameters=json.loads(parsed.stdout), history_run_id=history_id,
             history_origin='tracked_run' if history_id else ('legacy_external' if a.preset == 'append5plus5' else None),
             reuse_mode='raw_bam' if '--bam' in inputs else ('historical_graph_and_final_cdna' if a.preset == 'append5plus5' else 'full'),
             output_dir=str(run / 'results'), resources_requested=dict(cpus=a.cpus, memory=a.memory, partition=a.partition, time_limit=a.time_limit),
             memory_note='GNU time MaxRSS and Slurm MaxRSS are not a measured simultaneous process-tree total; no cgroup peak collector configured.')
    for value in [a.partition, a.memory, a.time_limit]:
        if not re.fullmatch(r'[A-Za-z0-9_:.-]+', value):
            raise ValueError('Invalid Slurm option')
    header = f'''#!/bin/bash
#SBATCH --job-name=pt-{a.preset}
#SBATCH --partition={a.partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={a.cpus}
#SBATCH --mem={a.memory}
#SBATCH --time={a.time_limit}
#SBATCH --output={run}/logs/slurm-%j.out
#SBATCH --error={run}/logs/slurm-%j.err
set -euo pipefail
export PATH={TOOLS}:/data/changq/miniforge3/bin:$PATH
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
'''
    (run / 'job.slurm').write_text(header + shlex.join([PYTHON, str(run / 'runner.py'), 'execute', str(run)]) + '\n')
    collector = header.replace(f'--cpus-per-task={a.cpus}', '--cpus-per-task=1').replace(f'--mem={a.memory}', '--mem=1G').replace(f'--time={a.time_limit}', '--time=00:10:00')
    (run / 'collect.slurm').write_text(collector + shlex.join([PYTHON, str(run / 'runner.py'), 'collect', str(run), '--retry']) + '\n')
    (run / 'command.sh').write_text('#!/bin/bash\n' + shlex.join([PYTHON, '-m', 'pantrans'] + args) + '\n')
    persist(run, d)
    print(run)


def execute(run):
    d = json.loads((run / 'run.json').read_text())
    if d['status'] not in ['prepared', 'submitted']:
        raise ValueError('Run already started; prepare a fresh run directory')
    d.update(status='running', started_at=now(), job_id=os.getenv('SLURM_JOB_ID', d.get('job_id')),
             node=os.uname().nodename)
    persist(run, d)
    rc = 1
    try:
        for name, sha in d['snapshot_files'].items():
            if digest(run / 'snapshot' / name) != sha:
                raise ValueError('Code snapshot changed: ' + name)
        # Hash inputs at runtime; metadata was validated during preparation.
        hashes = {}
        for entry in d['inputs'].values():
            current = identity(entry['path'])
            if any(entry[k] != current[k] for k in ['size', 'mtime_ns']):
                raise ValueError('Input changed after preparation: ' + entry['path'])
            if entry['path'] not in hashes:
                hashes[entry['path']] = digest(entry['path'])
            entry['sha256'] = hashes[entry['path']]
        env = dict(os.environ, PYTHONPATH=str(run / 'snapshot/src'), PYTHONDONTWRITEBYTECODE='1',
                   PATH=TOOLS + ':/data/changq/miniforge3/bin:' + os.environ['PATH'])
        environment = {name: {'path': shutil.which(name, path=env['PATH'])} for name in ['minimap2', 'samtools']}
        for name, item in environment.items():
            item['version'] = capture([item['path'], '--version']) if item['path'] else None
        environment['python'] = capture([PYTHON, '--version'])
        environment['packages'] = capture([PYTHON, '-m', 'pip', 'freeze'])
        environment['transdecoder_files'] = {str(p): digest(p) for p in Path(TD).rglob('*') if p.is_file()}
        save(run / 'environment.json', environment)
        d['analysis_started_at'] = now(); persist(run, d)
        start = time.monotonic()
        with (run / 'logs/stdout.log').open('w') as stdout, (run / 'logs/stderr.log').open('w') as stderr:
            rc = subprocess.call(['/usr/bin/time', '-v', '-o', str(run / 'logs/time.txt'),
                                  PYTHON, '-m', 'pantrans'] + d['argv'], cwd=run, env=env, stdout=stdout, stderr=stderr)
        d['elapsed_seconds'] = time.monotonic() - start
        text = (run / 'logs/time.txt').read_text()
        match = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
        d['gnu_time_max_rss_kib'] = int(match.group(1)) if match else None
    except Exception as exc:
        d['error'] = repr(exc)
    finally:
        d.update(status='succeeded' if rc == 0 else 'failed', exit_code=rc, finished_at=now())
        d['output_files'] = [identity(p) for p in (run / 'results').rglob('*') if p.is_file() and p.stat().st_size]
        persist(run, d)
    return rc


def collect(run, retry=False):
    d = json.loads((run / 'run.json').read_text())
    job = d.get('job_id')
    if not job:
        raise ValueError('No job ID recorded')
    fields = 'JobID,State,ExitCode,ElapsedRaw,TotalCPU,AllocCPUS,ReqMem,MaxRSS,Submit,Start,End,NodeList'
    rows = []
    for attempt in range(12 if retry else 1):
        result = capture(['sacct', '-j', str(job), '--units=K', '-P', '--format=' + fields])
        (run / 'logs/sacct.txt').write_text(result['output'])
        if result['exit_code'] == 0:
            rows = list(csv.DictReader(result['output'].splitlines(), delimiter='|'))
            primary = next((r for r in rows if r['JobID'] == str(job)), {})
            if primary.get('State') and primary['State'] not in ['RUNNING', 'PENDING', 'COMPLETING']:
                break
        if retry and attempt < 11:
            time.sleep(10)
    d = json.loads((run / 'run.json').read_text())
    d['slurm_accounting'] = rows; d['accounting_collected_at'] = now()
    primary = next((r for r in rows if r['JobID'] == str(job)), {})
    state = primary.get('State', '')
    d['accounting_complete'] = bool(state and state not in ['RUNNING', 'PENDING', 'COMPLETING'])
    d['slurm_elapsed_seconds'] = primary.get('ElapsedRaw') or None
    peaks = []
    for row in rows:
        match = re.fullmatch(r'(\d+(?:\.\d+)?)K?', row.get('MaxRSS', ''))
        if match:
            peaks.append(float(match.group(1)))
    d['slurm_step_max_rss_kib'] = max(peaks) if peaks else None
    if d['accounting_complete'] and state != 'COMPLETED':
        d['status'] = state.lower(); d['slurm_exit_code'] = primary.get('ExitCode')
    elif state == 'COMPLETED' and d['status'] not in ['succeeded', 'failed']:
        d['status'] = 'completed_missing_runtime_record'
    persist(run, d)
    print(json.dumps({'status': d['status'], 'accounting_complete': d['accounting_complete']}))


def submit(run):
    d = json.loads((run / 'run.json').read_text())
    if d['status'] != 'prepared':
        raise ValueError('Only prepared runs can be submitted')
    # Hold main job until both job IDs are durably recorded, avoiding startup races.
    job = subprocess.check_output(['sbatch', '--hold', '--parsable', str(run / 'job.slurm')], text=True).strip().split(';')[0]
    d.update(job_id=job, status='submitted', submitted_at=now()); persist(run, d)
    try:
        collector = subprocess.check_output(['sbatch', '--parsable', '--dependency=afterany:' + job,
                                             str(run / 'collect.slurm')], text=True).strip().split(';')[0]
        d['collector_job_id'] = collector; persist(run, d)
        subprocess.run(['scontrol', 'release', job], check=True)
    except Exception:
        print(f'Job {job} remains held if release failed. Repair collector/release manually.', file=sys.stderr)
        raise
    print(json.dumps({'job_id': job, 'collector_job_id': collector}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('preset', choices=['construct5', 'append5plus5', 'construct10'])
    p.add_argument('--base', type=Path, default=BASE)
    p.add_argument('--repo', type=Path, default=Path('/data/changq/Project/PanTrans'))
    p.add_argument('--runs', type=Path, default=BASE / 'experiment_runs')
    p.add_argument('--purpose', required=True)
    p.add_argument('--history-dir', type=Path)
    p.add_argument('--history-graph-cdna', type=Path)
    p.add_argument('--full', action='store_true', help='construct5: do not reuse legacy raw BAM')
    p.add_argument('--rescue', action=argparse.BooleanOptionalAction, default=None)
    p.add_argument('--cpus', type=int, default=32)
    p.add_argument('--partition', default='fatPartition')
    p.add_argument('--memory', default='240G')
    p.add_argument('--time-limit', default='2-00:00:00')
    for action in ['execute', 'submit', 'collect']:
        p = sub.add_parser(action); p.add_argument('run', type=Path)
        if action == 'collect': p.add_argument('--retry', action='store_true')
    a = parser.parse_args()
    if a.action == 'prepare':
        if a.cpus < 2: parser.error('--cpus must be >= 2')
        prepare(a)
    elif a.action == 'execute': return execute(a.run.resolve())
    elif a.action == 'submit': submit(a.run.resolve())
    else: collect(a.run.resolve(), a.retry)
    return 0


if __name__ == '__main__':
    sys.exit(main())
