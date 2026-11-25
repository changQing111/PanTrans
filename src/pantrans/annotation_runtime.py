"""Bounded annotation scheduling and thread-private indexed FASTA readers."""
import collections
import json
import os
import threading
import time
import shutil
import tempfile
from pathlib import Path
from contextlib import contextmanager

import pysam
from .prepared_bundle import prepared_file


@contextmanager
def local_annotation_inputs(prepared, output, stage):
    """Optional node-local copies for random cDNA/gDNA fetches during annotation."""
    started=time.perf_counter()
    sources={key:prepared_file(prepared,key) for key in ('cdna','gdna')}
    scratch=os.environ.get('PANTRANS_LOCAL_SCRATCH',tempfile.gettempdir())
    required=sum(path.stat().st_size for path in sources.values())
    if shutil.disk_usage(scratch).free < required+1024**3:
        raise OSError('Insufficient local scratch for annotation FASTAs; need '+str(required+1024**3)+' bytes')
    with tempfile.TemporaryDirectory(prefix='pantrans-fasta-',dir=scratch) as directory:
        local=Path(directory)
        files={}
        for key,path in sources.items():
            files[key]=key+'.fa'
            target=local/files[key]
            shutil.copyfile(path,target)
            index=Path(str(path)+'.fai')
            if index.is_file():shutil.copyfile(index,str(target)+'.fai')
            else:pysam.faidx(str(target))
        (local/'local_preprocess.json').write_text(json.dumps(dict(files=files)))
        (output/(stage+'.local_fasta.json')).write_text(json.dumps(dict(
            copied_bytes=required,copy_seconds=time.perf_counter()-started,
            sources={key:str(p.resolve()) for key,p in sources.items()}),indent=2))
        yield local


class FastaHandles:
    def __init__(self):
        self.local = threading.local()
        self.handles = []
        self.lock = threading.Lock()

    def get(self, path):
        if not hasattr(self.local, 'handles'):
            self.local.handles = {}
        key = str(path)
        if key not in self.local.handles:
            handle = pysam.FastaFile(key)
            self.local.handles[key] = handle
            with self.lock:
                self.handles.append(handle)
        return self.local.handles[key]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # The enclosing executor must finish before readers are closed.
        for handle in self.handles:
            handle.close()


def ordered_results(pool, fn, items, limit):
    """Refill a bounded queue while retaining deterministic output order."""
    if limit < 1:
        raise ValueError('Queue limit must be positive')
    items = iter(items)
    pending = collections.deque()
    for _ in range(limit):
        try:
            item = next(items)
        except StopIteration:
            break
        pending.append(pool.submit(fn, item))
    while pending:
        result = pending.popleft().result()
        try:
            item = next(items)
        except StopIteration:
            pass
        else:
            pending.append(pool.submit(fn, item))
        yield result


def write_progress(output, stage, completed, total, started):
    elapsed = time.monotonic() - started
    progress = dict(stage=stage, completed_clusters=completed, total_clusters=total,
                    elapsed_seconds=elapsed, time=time.time())
    temporary = output / 'progress.tmp.json'
    temporary.write_text(json.dumps(progress, indent=2))
    os.replace(temporary, output / 'progress.json')
    print(json.dumps(progress), flush=True)
