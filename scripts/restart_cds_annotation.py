"""Annotation-only restart; immutable graph inputs, persistent worker FASTA handles."""
import argparse,collections,concurrent.futures,gzip,hashlib,importlib.util,inspect,json,os,sqlite3,sys,threading,time,types
from pathlib import Path
import pysam
import tempfile, shutil, warnings


def load_pipeline(source):
 sys.path.insert(0,str(source/'src'))
 entry=source/'cds_pipeline.py'
 if not entry.is_file():entry=source/'src/pantrans/cds_construct.py'
 if not entry.is_file():entry=source/'scripts/experimental_cds_construct.py'
 spec=importlib.util.spec_from_file_location('original',entry);p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p);return p

class LocalTemporaryDirectory:
 def __init__(self,*args,**kwargs):
  # Each task has a private node-local directory; never use shared work for scratch.
  kwargs['dir']=os.environ.get('PANTRANS_LOCAL_SCRATCH','/tmp')
  self.path=tempfile.mkdtemp(*args,**kwargs)
 def __enter__(self):return self.path
 def __exit__(self,exc_type,exc,tb):
  try:shutil.rmtree(self.path)
  except OSError as error:
   warnings.warn('Temporary cleanup deferred: '+self.path+': '+str(error),RuntimeWarning)
  return False

class Handles:
 def __init__(self):self.local=threading.local();self.lock=threading.Lock();self.all=[];self.count=0
 def open(self,path,*args,**kwargs):
  if not hasattr(self.local,'cache'):self.local.cache={}
  key=str(path)
  if key not in self.local.cache:
   h=pysam.FastaFile(path,*args,**kwargs);self.local.cache[key]=h
   with self.lock:self.all.append(h);self.count+=1
  return Proxy(self.local.cache[key])
 def close(self):
  for h in self.all:
   if not h.closed:h.close()
class Proxy:
 def __init__(self,h):self.h=h
 def __enter__(self):return self
 def __exit__(self,*args):pass
 def __getattr__(self,k):return getattr(self.h,k)

class Pool:
 """Keep original writer order but refill bounded queue without batch barriers."""
 def __init__(self,workers):self.pool=concurrent.futures.ThreadPoolExecutor(max_workers=workers)
 def __enter__(self):return self
 def __exit__(self,*args):self.pool.shutdown(wait=True)
 def ordered(self,fn,items,limit):
  it=iter(items);pending=collections.deque()
  for _ in range(limit):
   try:pending.append(self.pool.submit(fn,next(it)))
   except StopIteration:break
  while pending:
   result=pending.popleft().result()
   try:pending.append(self.pool.submit(fn,next(it)))
   except StopIteration:pass
   yield result


def install(p,clusters,models,gt,output):
 p.tempfile=types.SimpleNamespace(TemporaryDirectory=LocalTemporaryDirectory,gettempdir=tempfile.gettempdir)
 handles=Handles();p.pysam=types.SimpleNamespace(**{k:getattr(pysam,k) for k in dir(pysam) if not k.startswith('__')});p.pysam.FastaFile=handles.open
 original=p.annotation_cluster;unions={}
 for stage,cs in clusters.items():
  for c in cs:unions.setdefault(c[0],set()).update(c)
 db=sqlite3.connect(output/'annotation_cache.sqlite',check_same_thread=False);db.execute('PRAGMA journal_mode=WAL');db.execute('CREATE TABLE cache(hub TEXT PRIMARY KEY,payload BLOB)');lock=threading.Lock();counts=collections.Counter()
 def cached(cluster,prepared,models,gene_tids,work,**kwargs):
  hub=cluster[0]
  with lock:row=db.execute('SELECT payload FROM cache WHERE hub=?',(hub,)).fetchone()
  if row:result=json.loads(gzip.decompress(row[0]));counts['cache_hits']+=1
  else:
   result=original([hub]+sorted(unions[hub]-{hub}),prepared,models,gene_tids,work)
   payload=gzip.compress(json.dumps(result).encode(),compresslevel=1)
   with lock:db.execute('INSERT INTO cache VALUES (?,?)',(hub,payload));db.commit();counts['computed_hubs']+=1
  h,candidates,audit,sams=result
  for m in candidates:
   m['exons']=[tuple(x) for x in m['exons']];m['cds']=[tuple(x) for x in m['cds']]
  qids={t for g in cluster[1:] for t in gene_tids[g] if models[t]['valid']};ids=qids|set(gene_tids[hub])
  return h,[m for m in candidates if m['source'] in ids],[a for a in audit if a[1] in qids],[s for s in sams if s.split('\t')[0] in qids]
 p.annotation_cluster=cached
 if hasattr(p,'ordered_results'):return handles,db,counts
 # Preserve the validated writer verbatim except scheduling and progress flushing.
 code=inspect.getsource(p.annotation_stage)
 old='''            for start in range(0,len(clusters),workers*4):
                futures=[pool.submit(annotation_cluster,c,prepared,models,gene_tids,work) for c in clusters[start:start+workers*4]]
                for future in futures:
                    hub,candidates,checks,sam_records=future.result();unique=deduplicate(candidates);stats['genes']+=bool(unique);stats['transcripts']+=len(unique)'''
 new='''            fn=lambda c: annotation_cluster(c,prepared,models,gene_tids,work)
            for index,result in enumerate(pool.ordered(fn,clusters,workers*4),1):
                if True:
                    hub,candidates,checks,sam_records=result;unique=deduplicate(candidates);stats['genes']+=bool(unique);stats['transcripts']+=len(unique)
                    if index % 1000 == 0:
                        for handle in (gff,cdna,sources,audit):handle.flush()
                        progress={'stage':stage,'completed_clusters':index,'total_clusters':len(clusters),'time':time.time()}
                        temp=output/'progress.tmp.json';temp.write_text(json.dumps(progress));os.replace(temp,output/'progress.json')
                        print(json.dumps(progress),flush=True)'''
 assert old in code;code=code.replace(old,new).replace('ThreadPoolExecutor(max_workers=workers)','RestartPool(workers)');p.RestartPool=Pool;exec(code,p.__dict__)
 return handles,db,counts


def main():
 a=argparse.ArgumentParser();a.add_argument('--source',type=Path,required=True);a.add_argument('--old-results',type=Path,required=True);a.add_argument('--prepared',type=Path,required=True);a.add_argument('--annotations',type=Path,help='Only required for schema-1 inputs');a.add_argument('--output',type=Path,required=True);a.add_argument('--workers',type=int,default=8);a.add_argument('--selection',type=Path);a.add_argument('--expected',type=Path);a.add_argument('--annotate-pre',action='store_true');args=a.parse_args();
 if args.workers<1:a.error('workers must be positive')
 stages=['pre','last'] if args.annotate_pre else ['last']
 args.output.mkdir(exist_ok=False);p=load_pipeline(args.source);start=time.monotonic()
 if args.selection:
  tasks=json.loads(args.selection.read_text())['tasks'];clusters={stage:[c for s,c in tasks if s==stage] for stage in stages};selected={g for cs in clusters.values() for c in cs for g in c}
 else:clusters={stage:[l.split() for l in (args.old_results/(stage+'.cluster')).open()] for stage in stages};selected=None
 checks={}
 for name in ['cds.raw.bam','cds.filtered.bam','graph.edges.tsv','pre.cluster','last.cluster','graph.summary.json']:
  path=args.old_results/name;assert path.is_file();st=path.stat();checks[name]={'path':str(path),'size':st.st_size,'mtime_ns':st.st_mtime_ns}
 (args.output/'reused_graph.json').write_text(json.dumps(checks,indent=2))
 # Rebuild only in-memory original-annotation metadata; graph_stage is never called.
 bundle=args.output/'metadata';bundle.mkdir();regions,models,gt=p.load_bundle(args.prepared,args.annotations.read_text().splitlines() if args.annotations else [],bundle,selected)
 handles,db,counts=install(p,clusters,models,gt,args.output)
 summaries={}
 try:
  if args.expected:
   expected=json.loads(args.expected.read_text());observed={};work=args.output/'checkwork';work.mkdir()
   def canon(result):
    h,c,au,s=result
    return {'hub':h,'models':sorted(p.deduplicate(c),key=lambda x:json.dumps(x,sort_keys=True)),'audit':sorted(au),'accepted_sam':sorted(s)}
   with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
    for stage,cs in clusters.items():
     futures=[(c[0],pool.submit(p.annotation_cluster,c,args.prepared,models,gt,work)) for c in cs]
     for h,f in futures:observed[stage+':'+h]=canon(f.result())
   assert json.loads(json.dumps(observed))==expected,'Biological result/SAM mismatch'
   (args.output/'equivalence.json').write_text(json.dumps({'passed':True,'tasks':len(observed),'reference':str(args.expected)}))
  for stage,cs in clusters.items():summaries[stage]=p.annotation_stage(stage,cs,args.prepared,models,gt,args.output,args.workers)
 finally:handles.close();db.close()
 for name,info in checks.items():
  st=Path(info['path']).stat();assert (st.st_size,st.st_mtime_ns)==(info['size'],info['mtime_ns']),name+' changed'
 (args.output/'summary.json').write_text(json.dumps({'annotation_only':True,'elapsed_seconds':time.monotonic()-start,'stages':summaries,'fasta_opens':handles.count,'cache':dict(counts),'original_representative_models_seeded':True},indent=2))
 print('ANNOTATION RESTART COMPLETE',flush=True)
if __name__=='__main__':main()
