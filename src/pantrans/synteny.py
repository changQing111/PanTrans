#!/usr/bin/env python3
"""Build all unordered variety pairs using indexed pre clusters and gene order.

Stages: index, shard, pair, collect. Reciprocal mode reads every directed edge
file once; pair tasks intersect sorted uint64 keys instead of rescanning TSVs.
The chaining and positional-support rules match the CS pilot.
"""
from __future__ import annotations

import argparse
from array import array
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import os
import re
from pathlib import Path
try:
    import resource
except ImportError:  # Windows can run the small correctness fixtures.
    resource = None
import sys
import time

import numpy as np
compiled_components = compiled_chains = None

VERSION = '1.0'
MAIN = {f'Chr{i}{s}' for i in range(1, 8) for s in 'ABD'}
META = np.dtype([('var','u1'),('chrom','u4'),('start','u8'),('end','u8'),('rank','u4'),('cluster','u4')])


def save(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def save_array(path, value):
    tmp = path.with_name(path.name + '.partial')
    with tmp.open('wb') as f:
        np.save(f, value, allow_pickle=False)
    tmp.replace(path)


def table(path, fields, rows):
    tmp = path.with_name(path.name + '.partial')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(rows)
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def measured(start):
    return dict(wall_seconds=time.monotonic()-start,
                max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else None)


def read_manifest(root):
    return json.loads((root/'inputs.json').read_text())


def load_index(root):
    m = read_manifest(root)
    names = (root/'index/gene_ids.txt').read_bytes().splitlines()
    meta = np.load(root/'index/genes.npy', mmap_mode='r', allow_pickle=False)
    copies = np.load(root/'index/copies.npy', mmap_mode='r', allow_pickle=False)
    return m, names, meta, copies


def index(args):
    start = time.monotonic()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root/'inputs.json').exists():
        raise ValueError('Index already completed; reuse it or choose a new output directory')
    prepared = args.prepared.resolve()
    from .prepared_bundle import manifest, prepared_file
    pm = manifest(str(prepared))
    if pm.get('schema_version') != 2:
        raise ValueError('synteny requires schema-2 prepared inputs')
    varieties = [v['variety'] for v in pm['varieties']]
    if len(varieties) < 2 or len(varieties) != len(set(varieties)) or len(varieties)>255:
        raise ValueError('synteny requires 2 to 255 distinct varieties')
    varieties = ['CS']+[v for v in varieties if v!='CS'] if 'CS' in varieties else varieties
    vmap = {v:i for i,v in enumerate(varieties)}
    bed = prepared_file(prepared, 'bed')
    pre = args.pre_cluster.resolve()
    edges = args.edges.resolve()
    if not edges.is_file():raise FileNotFoundError(edges)
    main = list(MAIN)
    if args.main_chroms:
        main = args.main_chroms.read_text().split()
        if not main or len(main)!=len(set(main)):raise ValueError('Main chromosome list must be nonempty and unique')
    ids, records, cmap, chroms = [], [], {}, []
    groups = defaultdict(list)
    with bed.open('rb') as f:
        for number,line in enumerate(f,1):
            parts=line.split()
            if len(parts)!=6:
                raise ValueError(f'BED line {number}: expected six columns')
            chrom,lo,hi,gene,var,strand=parts
            chrom,var=chrom.decode(),var.decode()
            a,b=int(lo),int(hi)
            if not 0<=a<b or strand not in (b'+',b'-'):
                raise ValueError(f'BED line {number}: invalid interval or strand')
            if chrom not in cmap:
                cmap[chrom]=len(chroms);chroms.append(chrom)
            if var not in vmap:raise ValueError('Unknown BED variety: '+var)
            i=len(ids);v=vmap[var];c=cmap[chrom]
            ids.append(gene);records.append((v,c,a,b,0,0));groups[v,c].append(i)
    if not ids:raise ValueError('Empty gene BED')
    if {int(r[0]) for r in records} != set(range(len(varieties))):
        raise ValueError('Every manifest variety must have BED genes')
    if len(ids)>=2**32:
        raise ValueError('Packed pair keys support fewer than 2**32 genes')
    genes = dict(zip(ids,range(len(ids))))
    if len(genes)!=len(ids):
        raise ValueError('Duplicate BED gene ID')
    meta=np.array(records,dtype=META);del records
    for members in groups.values():
        members.sort(key=lambda i:(int(meta[i]['start']),int(meta[i]['end']),ids[i]))
        meta['rank'][members]=np.arange(len(members),dtype=np.uint32)
    del groups
    counts=[np.zeros(len(varieties),dtype=np.uint32)]
    seen=np.zeros(len(ids),dtype=bool)
    with pre.open('rb') as f:
        for cid,line in enumerate(f,1):
            if not line.split() or any(g not in genes for g in line.split()):
                raise ValueError('Empty pre cluster or unknown gene at line '+str(cid))
            members=np.array([genes[g] for g in line.split()],dtype=np.uint32)
            if not len(members) or len(np.unique(members))!=len(members) or np.any(seen[members]):
                raise ValueError(f'Empty cluster or duplicate membership at pre line {cid}')
            seen[members]=True;meta['cluster'][members]=cid
            counts.append(np.bincount(meta['var'][members],minlength=len(varieties)).astype(np.uint32))
    if not np.all(seen):
        raise ValueError('BED and pre membership do not cover the same gene set')
    copies=np.array(counts,dtype=np.uint32);del counts,genes,seen
    idx=root/'index';idx.mkdir(exist_ok=True)
    (idx/'gene_ids.txt').write_bytes(b'\n'.join(ids)+b'\n')
    save_array(idx/'genes.npy',meta);save_array(idx/'copies.npy',copies)
    # O(number of clusters x varieties) sizing, without expanding high-copy pairs.
    potential=copies.astype(np.uint64).T @ copies.astype(np.uint64)
    seedcounts=copies.copy();seedcounts[seedcounts>args.max_copy]=0
    seedpotential=seedcounts.astype(np.uint64).T @ seedcounts.astype(np.uint64)
    pairs=[]
    for i,j in itertools.combinations(range(len(varieties)),2):
        pairs.append(dict(index=len(pairs),left=varieties[i],right=varieties[j],left_index=i,right_index=j,
            name=f'pair_{len(pairs)+1:04d}',co_membership_pairs=int(potential[i,j]),
            low_copy_co_membership_pairs=int(seedpotential[i,j])))
    m=dict(status='completed',version=VERSION,prepared=str(prepared),edges=str(edges),main_chroms=sorted(main),
           bed=str(bed),pre_cluster=str(pre),varieties=varieties,chromosomes=chroms,nodes=len(ids),
           clusters=len(copies)-1,pairs=pairs,candidate_mode=args.mode,
           parameters=dict(max_copy=args.max_copy,rank_distance=args.distance,min_anchors=args.min_anchors,
                           lift_distance=args.lift_distance),
           source_sha256=dict(bed=sha256(bed),pre_cluster=sha256(pre)),
           code_sha256=sha256(Path(__file__)),index_resources=measured(start),
           candidate_semantics='Same pre cluster AND direct edges in both directions' if args.mode=='reciprocal' else
                               'All same-pre-cluster low-copy pairs; high-copy candidates queried locally only')
    save(root/'inputs.json',m)
    table(root/'pairs.tsv',['index','left','right','co_membership_pairs','low_copy_co_membership_pairs'],
          ((p['index'],p['left'],p['right'],p['co_membership_pairs'],p['low_copy_co_membership_pairs']) for p in pairs))
    print(json.dumps(dict(stage='index',nodes=len(ids),clusters=m['clusters'],pairs=len(pairs),**measured(start))),flush=True)


def binary_lines(path, digest):
    """Hash large blocks while parsing lines, avoiding a second full file read."""
    with path.open('rb') as f:
        tail=b''
        for chunk in iter(lambda:f.read(8*1024*1024),b''):
            digest.update(chunk)
            lines=(tail+chunk).split(b'\n');tail=lines.pop()
            yield from lines
        if tail:
            yield tail


def shard(args):
    """Read the complete construct/append graph once, partition by target variety."""
    start=time.monotonic();root=args.output.resolve()
    m,names,meta,_=load_index(root)
    for label in ('bed','pre_cluster'):
        if sha256(Path(m[label])) != m['source_sha256'][label]:
            raise ValueError('Synteny input changed after indexing: '+label)
    shard_root=root/'shards';shard_root.mkdir(exist_ok=False)
    ids=dict(zip(names,range(len(names))))
    buffers={};counts=Counter();kept=Counter();h=hashlib.sha256()
    paths={}
    for pair_info in m['pairs']:
        for v in (pair_info['left_index'],pair_info['right_index']):
            directory=shard_root/str(v);directory.mkdir(exist_ok=True)
            key=(v,pair_info['name']);path=directory/(pair_info['name']+'.raw')
            buffers[key]=array('Q');paths[key]=path;path.touch(exist_ok=False)
    def flush(key):
        with paths[key].open('ab') as f:np.asarray(buffers[key],dtype='<u8').tofile(f)
        buffers[key]=array('Q')
    lookup={frozenset((p['left_index'],p['right_index'])):p['name'] for p in m['pairs']}
    for number,line in enumerate(binary_lines(Path(m['edges']),h),1):
        parts=line.split()
        if len(parts)!=2 or any(g not in ids for g in parts):
            raise ValueError('Unknown gene or malformed graph edge at line '+str(number))
        q,t=ids[parts[0]],ids[parts[1]]
        a,b=int(meta[q]['var']),int(meta[t]['var']);counts[b]+=1
        if a==b or meta[q]['cluster']!=meta[t]['cluster']:continue
        left,right=(q,t) if a<b else (t,q)
        key=(b,lookup[frozenset((a,b))])
        buffers[key].append((left<<32)|right);kept[b]+=1
        if len(buffers[key])>=65536:flush(key)
    for key in buffers:flush(key)
    for key in paths:
        v,pair_name=key;directory=shard_root/str(v)
        raw=directory/(pair_name+'.raw')
        keys=np.unique(np.fromfile(raw,dtype='<u8'))
        save_array(directory/(pair_name+'.npy'),keys)
        raw.unlink()
    for v,variety in enumerate(m['varieties']):
        save(shard_root/str(v)/'summary.json',dict(status='completed',variety=variety,
            input_edges=counts[v],kept_cross_variety_same_cluster_edges=kept[v]))
    save(root/'shard.summary.json',dict(status='completed',edge_sha256=h.hexdigest(),
        input_edges=sum(counts.values()),kept_edges=sum(kept.values()),**measured(start)))


def components(points,distance,minimum):
    """Bounded rank-neighbor connected components, no compiler/JCVI dependency."""
    points=sorted(points);parent=list(range(len(points)))
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    for i,(x,y,_) in enumerate(points):
        for j in range(i-1,-1,-1):
            if x-points[j][0]>distance:break
            if abs(y-points[j][1])<=distance:
                a,b=find(i),find(j)
                if a!=b:parent[b]=a
    groups=defaultdict(list)
    for i,p in enumerate(points):groups[find(i)].append(p)
    return [p for p in groups.values() if len(p)>1 and
            min(len({x for x,y,w in p}),len({y for x,y,w in p}))>=minimum]


def monotonic_chains(points,distance,minimum):
    """Exact pilot tie rules; the bounded x window avoids quadratic full scans."""
    remaining=set(points);result=[]
    while len(remaining)>=minimum:
        best=[];best_direction=1;p=sorted(remaining)
        for direction in (1,-1):
            scores=[1]*len(p);previous=[-1]*len(p)
            for i,(x,y,w) in enumerate(p):
                for j in range(i-1,-1,-1):
                    dx=x-p[j][0]
                    if dx>distance:break
                    dy=direction*(y-p[j][1])
                    if dx>0 and 0<dy<=distance and scores[j]+1>scores[i]:
                        scores[i]=scores[j]+1;previous[i]=j
            end=max(range(len(p)),key=lambda i:scores[i]);chain=[]
            while end>=0:chain.append(p[end]);end=previous[end]
            chain.reverse()
            if len(chain)>len(best):best=chain;best_direction=direction
        if len(best)<minimum:break
        result.append((best_direction,best));remaining.difference_update(best)
    return result


def find_position_support(x,y,key,block_points,spatial,lift):
    """Exact pilot rule after a seed-rank index narrows the blocks to inspect."""
    possible=set()
    for bucket in range((x-lift)//(lift+1),(x+lift)//(lift+1)+1):
        possible.update(spatial.get((key,bucket),()))
    choices=[]
    for bi in possible:
        bid,direction,points=block_points[key][bi]
        if x<points[0][0]-lift or x>points[-1][0]+lift:continue
        pos=bisect_left(points,(x,-1));neighbors=points[max(0,pos-1):pos+1]
        distances=[abs(x-u)+abs(y-w) for u,w in neighbors if abs(x-u)<=lift and abs(y-w)<=lift]
        if not distances:continue
        if 0<pos<len(points):
            if not min(points[pos-1][1],points[pos][1])<=y<=max(points[pos-1][1],points[pos][1]):continue
        elif pos==0:
            if direction*(y-points[0][1])>0:continue
        elif direction*(y-points[-1][1])<0:continue
        choices.append((min(distances),bid))
    choices.sort()
    if not choices:return None,False
    if len(choices)>1 and choices[0][0]==choices[1][0]:return None,True
    return choices[0],False


def pair(args):
    start=time.monotonic();root=args.output.resolve();m,names,meta,copies=load_index(root)
    pp=m['pairs'][args.index];lv,rv=pp['left_index'],pp['right_index'];left,right=pp['left'],pp['right']
    dest=root/'pairs'/pp['name'];dest.mkdir(parents=True,exist_ok=True)
    if (dest/'summary.json').exists():
        raise ValueError('Pair already completed')
    if m['candidate_mode']!='reciprocal':
        raise ValueError('This implementation currently supports reciprocal candidates only')
    for var in (lv,rv):
        if not (root/'shards'/str(var)/'summary.json').exists():raise ValueError('Required shard incomplete')
    read_start=time.monotonic()
    aa=np.load(root/'shards'/str(lv)/(pp['name']+'.npy'),mmap_mode='r')
    bb=np.load(root/'shards'/str(rv)/(pp['name']+'.npy'),mmap_mode='r')
    keys=np.intersect1d(aa,bb,assume_unique=True)
    del aa,bb
    # Keep the complete candidate universe compact. High-copy contig relations
    # can exceed 20 million per pair and must not become Python tuple objects.
    save_array(dest/'candidate_pairs.npy',keys)
    n_candidates=len(keys)
    left_ids=(keys>>np.uint64(32)).astype(np.uint32)
    right_ids=(keys&np.uint64(0xffffffff)).astype(np.uint32)
    all_clusters=meta['cluster'][left_ids]
    if not np.array_equal(all_clusters,meta['cluster'][right_ids]):raise ValueError('Cross-cluster candidate')
    placed_mask=np.isin(meta['chrom'][left_ids],[i for i,c in enumerate(m['chromosomes']) if c in m['main_chroms']])
    placed_mask &= np.isin(meta['chrom'][right_ids],[i for i,c in enumerate(m['chromosomes']) if c in m['main_chroms']])
    ca,cb=copies[all_clusters,lv],copies[all_clusters,rv]
    low_mask=(ca<=m['parameters']['max_copy'])&(cb<=m['parameters']['max_copy'])
    contig_count=int(n_candidates-np.count_nonzero(placed_mask))
    all_counts=dict(one_to_one=int(np.count_nonzero((ca==1)&(cb==1))),
                    low_copy=int(np.count_nonzero(low_mask&~((ca==1)&(cb==1)))),
                    high_copy=int(n_candidates-np.count_nonzero(low_mask)))
    del ca,cb,all_clusters,low_mask
    export_all=getattr(args,'export_candidates',False) or left=='CS'
    # Only main-chromosome candidates participate in block construction or
    # positional support. All other candidates remain in candidate_pairs.npy.
    select=np.arange(n_candidates) if export_all else np.flatnonzero(placed_mask)
    decoded=[(int(left_ids[i]),int(right_ids[i])) for i in select]
    del select,left_ids,right_ids,placed_mask,keys
    decoded.sort(key=lambda p:(names[p[0]],names[p[1]]))
    read_seconds=time.monotonic()-read_start
    def gid(i):return names[i].decode()
    p=m['parameters'];dist,minimum,lift=p['rank_distance'],p['min_anchors'],p['lift_distance']
    chroms=m['chromosomes'];rank=meta['rank'];chrom=meta['chrom'];cluster=meta['cluster']
    main=np.array([c in m['main_chroms'] for c in chroms],dtype=bool)
    bychrom=defaultdict(list);lookup={};candidates=[];tiers=Counter();seed_count=0
    for a,b in decoded:
        cid=int(cluster[a]);assert cid==cluster[b] and meta[a]['var']==lv and meta[b]['var']==rv
        na,nb=int(copies[cid,lv]),int(copies[cid,rv]);placed=bool(main[chrom[a]] and main[chrom[b]])
        low=max(na,nb)<=p['max_copy'];tier='one_to_one' if na==nb==1 else 'low_copy' if low else 'high_copy'
        row=(a,b,cid,na,nb,tier,placed,placed and low);candidates.append(row);tiers[tier]+=1
        key=chroms[int(chrom[a])],chroms[int(chrom[b])];x,y=int(rank[a]),int(rank[b])
        lookup[key,x,y]=row
        if placed and low:bychrom[key].append((x,y,1));seed_count+=1
    del decoded
    candidate_seconds=time.monotonic()-start-read_seconds
    chain_start=time.monotonic();blocks=[];anchors=[];used=set();bp=defaultdict(list);spatial=defaultdict(set);components_count=0
    for key,points in sorted(bychrom.items()):
        comps=(compiled_components(points,dist,minimum) if compiled_components else
                    components(points,dist,minimum));components_count+=len(comps)
        for comp in sorted(comps):
            for direction,chain in (compiled_chains(comp,dist,minimum) if compiled_chains else monotonic_chains(comp,dist,minimum)):
                bid=f'{left}_{right}_B{len(blocks)+1:06d}';ids=[]
                for x,y,_ in chain:
                    a,b,cid,na,nb,tier,placed,seed=lookup[key,x,y]
                    anchors.append((bid,a,b,key[0],key[1],x,y,cid,tier));used.add((a,b));ids.append((a,b,tier))
                wheat=[re.fullmatch(r'Chr([1-7])[ABD]',c) for c in key]
                kind=('corresponding_chromosome' if key[0]==key[1] else
                      'same_homoeologous_group' if all(wheat) and wheat[0][1]==wheat[1][1] else 'other_chromosome')
                blocks.append((bid,key[0],key[1],'+' if direction==1 else '-',len(chain),sum(t=='one_to_one' for a,b,t in ids),
                    min(int(meta[a]['start']) for a,b,t in ids),max(int(meta[a]['end']) for a,b,t in ids),
                    min(int(meta[b]['start']) for a,b,t in ids),max(int(meta[b]['end']) for a,b,t in ids),kind))
                bi=len(bp[key]);xy=[(x,y) for x,y,_ in chain];bp[key].append((bid,direction,xy))
                for x,y in xy:spatial[key,x//(lift+1)].add(bi)
    chain_seconds=time.monotonic()-chain_start
    lift_start=time.monotonic();supported=[];ambiguous=0;contig=0
    for a,b,cid,na,nb,tier,placed,seed in candidates:
        if (a,b) in used:continue
        if not placed:contig+=1;continue
        x,y=int(rank[a]),int(rank[b]);key=chroms[int(chrom[a])],chroms[int(chrom[b])]
        choice,tie=find_position_support(x,y,key,bp,spatial,lift)
        ambiguous+=tie
        if choice:supported.append((choice[1],a,b,cid,tier,choice[0]))
    lift_seconds=time.monotonic()-lift_start
    write_start=time.monotonic()
    candidate_filename='candidates.tsv' if export_all else 'placed_candidates.tsv'
    table(dest/candidate_filename,['reference_gene','target_gene','pre_cluster','reference_copies','target_copies','tier','main_chromosomes','seed'],
          ((gid(a),gid(b),cid,na,nb,tier,int(placed),int(seed)) for a,b,cid,na,nb,tier,placed,seed in candidates))
    table(dest/'blocks.tsv',['block_id','reference_chrom','target_chrom','orientation','anchors','one_to_one_anchors','reference_start0','reference_end0','target_start0','target_end0','chromosome_relation'],blocks)
    table(dest/'anchors.tsv',['block_id','reference_gene','target_gene','reference_chrom','target_chrom','reference_rank','target_rank','pre_cluster','tier'],
          ((bid,gid(a),gid(b),ac,bc,x,y,cid,tier) for bid,a,b,ac,bc,x,y,cid,tier in anchors))
    table(dest/'position_supported_pairs.tsv',['block_id','reference_gene','target_gene','pre_cluster','tier','rank_manhattan_distance'],
          ((bid,gid(a),gid(b),cid,tier,d) for bid,a,b,cid,tier,d in supported))
    c=Counter(b[-1] for b in blocks)
    s=dict(status='completed',pair_index=args.index,reference=left,variety=right,pair=pp['name'],
           candidate_pairs=n_candidates,seed_pairs=seed_count,seed_components=components_count,
           blocks=len(blocks),anchors=len(anchors),one_to_one_anchors=sum(a[-1]=='one_to_one' for a in anchors),
           reference_anchor_genes=len({a[1] for a in anchors}),target_anchor_genes=len({a[2] for a in anchors}),
           corresponding_chromosome_blocks=c['corresponding_chromosome'],same_homoeologous_group_blocks=c['same_homoeologous_group'],
           other_chromosome_blocks=c['other_chromosome'],reverse_blocks=sum(b[3]=='-' for b in blocks),
           position_supported_pairs=len(supported),high_copy_supported_pairs=sum(s[-2]=='high_copy' for s in supported),
           ambiguous_position_pairs=ambiguous,unplaced_contig_candidates=contig_count,
           candidate_tiers=all_counts,candidate_table=candidate_filename,candidate_table_rows=len(candidates),
           compiled_kernels=compiled_components is not None,
           timings=dict(read_and_intersect_seconds=read_seconds,candidate_seconds=candidate_seconds,
                        chaining_seconds=chain_seconds,position_support_seconds=lift_seconds,write_seconds=time.monotonic()-write_start),
           **measured(start))
    assert len(used)==len(anchors) and all(b[4]>=minimum for b in blocks)
    save(dest/'summary.json',s);print(json.dumps(s),flush=True)


def collect(args):
    root=args.output.resolve();m=read_manifest(root);rows=[]
    for pp in m['pairs']:
        d=root/'pairs'/pp['name'];s=json.loads((d/'summary.json').read_text())
        if s['status']!='completed' or s['pair_index']!=pp['index']:
            raise ValueError(f'Incomplete pair: {pp["name"]}')
        for filename,count in [(s.get('candidate_table','candidates.tsv'),s.get('candidate_table_rows',s['candidate_pairs'])),('blocks.tsv',s['blocks']),
                               ('anchors.tsv',s['anchors']),('position_supported_pairs.tsv',s['position_supported_pairs'])]:
            with (d/filename).open('rb') as f:actual=sum(1 for _ in f)-1
            if actual!=count:raise ValueError(f'Row count mismatch: {d/filename}')
        if (d/'candidate_pairs.npy').exists():
            if len(np.load(d/'candidate_pairs.npy',mmap_mode='r'))!=s['candidate_pairs']:
                raise ValueError('Binary candidate row count mismatch')
        rows.append(s)
    fields=[k for k in rows[0] if k not in ['timings','candidate_tiers']]
    table(root/'summary.tsv',fields,([s.get(k,'') for k in fields] for s in rows))
    result=dict(status='completed',comparisons=len(rows),nodes=m['nodes'],pre_clusters=m['clusters'],
                total_blocks=sum(s['blocks'] for s in rows),total_anchor_pairs=sum(s['anchors'] for s in rows),
                total_position_supported_pairs=sum(s['position_supported_pairs'] for s in rows),
                sum_pair_wall_seconds=sum(s['wall_seconds'] for s in rows),
                candidate_mode=m['candidate_mode'],parameters=m['parameters'],pairs=rows)
    save(root/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k!='pairs'}),flush=True)
    return result


def export(args):
    """Stream full candidate relations to text only when that view is needed."""
    root=args.output.resolve();m,names,meta,copies=load_index(root)
    pp=m['pairs'][args.index];d=root/'pairs'/pp['name']
    target=args.destination or d/'candidates.export.tsv'
    if target.exists():raise ValueError('Export destination already exists')
    if not target.parent.exists():raise ValueError('Export parent directory must exist')
    keys=np.load(d/'candidate_pairs.npy',mmap_mode='r',allow_pickle=False)
    chroms=m['chromosomes'];limit=m['parameters']['max_copy']
    def rows():
        for key in keys:
            a,b=int(key)>>32,int(key)&0xffffffff;cid=int(meta[a]['cluster'])
            na,nb=int(copies[cid,pp['left_index']]),int(copies[cid,pp['right_index']])
            placed=chroms[int(meta[a]['chrom'])] in m['main_chroms'] and chroms[int(meta[b]['chrom'])] in m['main_chroms']
            low=max(na,nb)<=limit;tier='one_to_one' if na==nb==1 else 'low_copy' if low else 'high_copy'
            yield names[a].decode(),names[b].decode(),cid,na,nb,tier,int(placed),int(placed and low)
    table(target,['reference_gene','target_gene','pre_cluster','reference_copies','target_copies','tier','main_chromosomes','seed'],rows())
    print(json.dumps(dict(pair=pp['name'],rows=len(keys),destination=str(target),order='packed integer gene ID; table order differs from lexicographic pilot tables')))


def configure_parser(p):
    p.add_argument('--prepared',type=Path,help='Schema-2 prepared directory')
    p.add_argument('--results',type=Path,help='construct/append results; uses pre.cluster and graph.edges.tsv')
    p.add_argument('--pre-cluster',type=Path,help='Explicit pre.cluster (alternative to --results)')
    p.add_argument('--edges',type=Path,help='Complete directed graph edges (alternative to --results)')
    p.add_argument('--main-chroms',type=Path,help='Main chromosome list; default wheat Chr1A through Chr7D')
    p.add_argument('-o','--output',type=Path,required=True,help='New output directory; reuse for later stages')
    p.add_argument('--max-copy',type=int,default=3,help='Maximum copies per variety in a pre cluster for seed eligibility')
    p.add_argument('--distance',type=int,default=20,help='Maximum adjacent anchor gene-rank distance on both axes')
    p.add_argument('--min-anchors',type=int,default=5,help='Minimum chain anchors')
    p.add_argument('--lift-distance',type=int,default=10,help='Maximum per-axis rank distance for additional position support')
    p.add_argument('--stage',choices=['run','index','shard','pair','collect','export'],default='run')
    p.add_argument('--pair-index',dest='index',type=int,help='Zero-based pair index for pair/export stages')
    p.add_argument('--export-candidates',action='store_true',help='Include unplaced candidates in text output (always retained in binary)')
    p.add_argument('--destination',type=Path,help='New candidate TSV destination for export stage')
    return p


def validate_args(args, parser):
    if args.stage in ('run','index'):
        if not args.prepared:parser.error('--prepared is required for run/index')
        if args.results:
            if args.pre_cluster or args.edges:parser.error('--results conflicts with --pre-cluster/--edges')
            args.pre_cluster=args.results/'pre.cluster';args.edges=args.results/'graph.edges.tsv'
        if not args.pre_cluster or not args.edges:parser.error('Provide --results or both --pre-cluster and --edges')
    if args.stage in ('pair','export') and (args.index is None or args.index<0):
        parser.error('--pair-index must be a nonnegative integer for pair/export')
    if min(args.max_copy,args.distance)<1 or args.min_anchors<2 or args.lift_distance<0:
        parser.error('max-copy/distance >=1, min-anchors >=2, lift-distance >=0 required')
    if args.destination and args.stage!='export':parser.error('--destination is for export stage only')


def run(args):
    root=args.output.resolve()
    if args.stage in ('run','index'):
        if root.exists():raise ValueError('Output directory already exists; use a new directory')
        args.mode='reciprocal'
        index(args)
    else:
        if not (root/'inputs.json').is_file():raise ValueError('Synteny index is missing; run index first')
    if args.stage in ('pair','export'):
        if args.index>=len(read_manifest(root)['pairs']):raise ValueError('Pair index outside available comparisons')
    if args.stage=='run':
        shard(args)
        for pp in read_manifest(root)['pairs']:
            args.index=pp['index'];pair(args)
        return collect(args)
    if args.stage in ('shard','pair','collect','export'):
        return globals()[args.stage](args)


def main(argv=None):
    p=configure_parser(argparse.ArgumentParser(description=__doc__))
    args=p.parse_args(argv);validate_args(args,p);run(args)
    return 0
