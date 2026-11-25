"""Content-addressed pre-dedup transfer evidence, with read-only parent caches."""
import collections
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import tempfile
import threading
import zlib
import Bio
import pysam

from .prepared_bundle import prepared_file


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def evidence_policy(functions):
    executable = shutil.which('minimap2')
    if executable is None:
        raise ValueError('minimap2 is required for annotation')
    return fingerprint(dict(schema=1, biopython=Bio.__version__, pysam=pysam.__version__, code=[inspect.getsource(fn) for fn in functions],
                            minimap2_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
                            minimap2_version=subprocess.check_output([executable, '--version'], text=True).strip()))


class EvidenceCache:
    """Stage writes on node-local storage; publish only after successful annotation."""
    def __init__(self, destination, sources=()):
        self.destination = Path(destination)
        self.sources = [Path(p) for p in sources]
        self.lock = threading.Lock()
        self.counts = collections.Counter()

    def __enter__(self):
        if self.destination.exists():
            raise FileExistsError(self.destination)
        self.temp = tempfile.TemporaryDirectory(prefix='pantrans-evidence-',
            dir=os.environ.get('PANTRANS_LOCAL_SCRATCH'))
        self.path = Path(self.temp.name)/'cache.sqlite'
        self.db = sqlite3.connect(self.path, check_same_thread=False)
        self.parents = []
        try:
            self.db.execute('CREATE TABLE evidence (key TEXT PRIMARY KEY, payload BLOB NOT NULL)')
            self.db.execute('PRAGMA user_version=1')
            for index,path in enumerate(self.sources):
                if not path.is_file():raise FileNotFoundError(path)
                local=Path(self.temp.name)/('parent-'+str(index)+'.sqlite')
                shutil.copyfile(path,local)
                # Never create missing source caches or write to historical results.
                db = sqlite3.connect(local.resolve().as_uri()+'?mode=ro', uri=True, check_same_thread=False)
                self.parents.append(db)
                if db.execute('PRAGMA user_version').fetchone()[0] != 1:
                    raise ValueError('Unsupported evidence cache: '+str(path))
                db.execute('SELECT key,payload FROM evidence LIMIT 0')
            return self
        except BaseException:
            self.__exit__(*__import__('sys').exc_info())
            raise

    def get(self, keys):
        found = {}
        with self.lock:
            for tid, key in keys.items():
                row = self.db.execute('SELECT payload FROM evidence WHERE key=?', (key,)).fetchone()
                if row is None:
                    for db in self.parents:
                        row = db.execute('SELECT payload FROM evidence WHERE key=?', (key,)).fetchone()
                        if row is not None:
                            self.db.execute('INSERT OR IGNORE INTO evidence VALUES (?,?)', (key,row[0]))
                            break
                if row is not None:
                    value=json.loads(zlib.decompress(row[0]))
                    if value['candidate'] is not None:
                        for name in ('exons','cds'):
                            value['candidate'][name]=[tuple(x) for x in value['candidate'][name]]
                    found[tid]=value
            self.counts['hits'] += len(found)
            self.counts['misses'] += len(keys)-len(found)
        return found

    def put(self, keys, values):
        with self.lock:
            self.db.executemany('INSERT OR REPLACE INTO evidence VALUES (?,?)',
                [(keys[tid],zlib.compress(json.dumps(value,sort_keys=True).encode(),1))
                 for tid,value in values.items()])
            self.counts['computed_pairs'] += len(values)

    def __exit__(self, exc_type, exc, tb):
        try:
            if exc_type is None:self.db.commit()
            self.db.close()
            for db in self.parents:db.close()
            if exc_type is None:
                temporary=self.destination.with_suffix('.tmp.sqlite')
                shutil.copyfile(self.path,temporary)
                os.replace(temporary,self.destination)
        finally:
            self.temp.cleanup()


def cached_cluster(cluster, prepared, models, gene_tids, work, handles, cache, policy, compute):
    hub=cluster[0]
    queries=[t for g in cluster[1:] for t in gene_tids[g] if models[t]['valid']]
    if not queries:
        return compute(cluster,prepared,models,gene_tids,work,fasta_handles=handles)
    gdna=handles.get(prepared_file(prepared,'gdna'))
    cdna=handles.get(prepared_file(prepared,'cdna'))
    target=fingerprint([policy,hub,gdna.fetch(hub).upper(),
                        [(t,models[t]) for t in gene_tids[hub]]])
    keys={t:fingerprint([target,t,models[t],cdna.fetch(t)]) for t in queries}
    evidence=cache.get(keys)
    missing=[t for t in queries if t not in evidence]
    _,candidates,checks,sams=compute(cluster,prepared,models,gene_tids,work,
                                    fasta_handles=handles,query_ids=missing)
    seeds=[m for m in candidates if m['source'] in gene_tids[hub]]
    transferred={m['source']:m for m in candidates if m['source'] not in gene_tids[hub]}
    sam_by_query={sam.split('\t',1)[0]:sam for sam in sams}
    values={t:dict(candidate=transferred.get(t),status=status,sam=sam_by_query.get(t))
            for _,t,status in checks}
    if set(values)!=set(missing):raise ValueError('Incomplete query evidence returned by annotation')
    cache.put(keys,values)
    evidence.update(values)
    # Rebuild in original query order, preserving duplicate biological evidence.
    return (hub, seeds+[evidence[t]['candidate'] for t in queries if evidence[t]['candidate'] is not None],
            [[hub,t,evidence[t]['status']] for t in queries],
            [evidence[t]['sam'] for t in queries if evidence[t]['sam'] is not None])
