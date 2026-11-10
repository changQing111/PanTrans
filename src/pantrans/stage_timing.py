"""Persist stage timings as stages finish, including failed stages."""
from contextlib import contextmanager
import json
import os
import time


class StageTimer:
    def __init__(self, output):
        self.output=output
        self.stages={}

    @contextmanager
    def stage(self, name):
        started=time.perf_counter()
        print(json.dumps(dict(stage=name,event='start')),flush=True)
        status='completed'
        try:
            yield
        except BaseException:
            status='failed'
            raise
        finally:
            self.stages[name]=dict(elapsed_seconds=time.perf_counter()-started,status=status)
            if self.output.is_dir():
                path=self.output/'timings.tmp.json'
                path.write_text(json.dumps(self.stages,indent=2))
                os.replace(path,self.output/'timings.json')
            print(json.dumps(dict(stage=name,**self.stages[name])),flush=True)
