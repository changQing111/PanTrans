"""Frozen-raw-BAM five-variety threshold/structural sweep; outputs are experimental."""
import argparse,collections,csv,importlib.util,json,sys,time
from pathlib import Path
import pysam
from pantrans.cds_filter_experiment import compensated_indels,compatible_junctions

SCENARIOS=['strict95','strict90','indel12_cov90','junction3_cov90','combined_cov90']

def details(read,qlen,blocks,p):
    qp=0;rp=read.reference_start;frame=True
    for op,n in read.cigartuples or []:
        if op in (0,7,8):qp+=n;rp+=n
        elif op==1:
            if n%3 and max(qp,0)<min(qp+n,qlen):frame=False
            qp+=n
        elif op==2:
            if n%3 and any(max(rp,a)<min(rp+n,b) for a,b in blocks):frame=False
            rp+=n
        elif op==3:rp+=n
        elif op in (4,5):qp+=n
    pieces,junctions,_,_,den,_=p.layout(read)
    projected=[(ra+max(qa,0)-qa,ra+min(qb,qlen)-qa) for qa,qb,ra,rb in pieces if max(qa,0)<min(qb,qlen)]
    if not den or not projected:return None
    shared=p.overlap(projected,blocks);tlen=sum(b-a for a,b in blocks)
    lo=min(a for a,b in projected);hi=max(b for a,b in projected)
    targets={(b,c) for (a,b),(c,d) in zip(blocks,blocks[1:])}
    observed={j for j in junctions if lo<j[0] and j[1]<hi}
    splice=observed<=targets and all(pair in observed for pair in targets if lo<pair[0] and pair[1]<hi)
    return frame,splice,1-read.get_tag('NM')/den,shared/qlen,shared/tlen,targets,observed,lo,hi


def main():
    a=argparse.ArgumentParser();a.add_argument('--prepared',type=Path,required=True);a.add_argument('--history',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args()
    args.output.mkdir(exist_ok=False);started=time.monotonic()
    from pantrans import cds_construct as p
    manifest=json.loads(next(args.prepared.glob('*_preprocess.json')).read_text())['files']
    with (args.prepared/manifest['genes']).open() as f:regions={r['gene']:(int(r['body_start0']),r['variety']) for r in csv.DictReader(f,delimiter='\t')}
    models={};targets=collections.defaultdict(dict)
    with (args.prepared/manifest['models']).open() as f:
        for line in f:
            m=json.loads(line);models[m['transcript']]=(m['gene'],m['cds_len'])
            if m['valid']:
                offset=regions[m['gene']][0];blocks=tuple((a-offset,b-offset) for a,b in m['cds'])
                targets[m['gene']][blocks]=m['transcript']
    counts=collections.Counter();record_counts=collections.Counter();deltas={s:set() for s in SCENARIOS[1:]};examples={s:[] for s in SCENARIOS[2:]}
    with pysam.AlignmentFile(args.history/'cds.raw.bam') as bam:
        for read in bam:
            counts['raw_records']+=1
            if read.is_unmapped or read.is_supplementary or read.is_reverse or not read.has_tag('NM'):continue
            qgene,qlen=models[read.query_name];target=read.reference_name;passes=set();hit_examples={}
            for blocks,tid in targets[target].items():
                d=details(read,qlen,blocks,p)
                if d is None:continue
                frame,splice,pid,qcov,tcov,introns,observed,lo,hi=d
                baseline=frame and splice and min(pid,qcov,tcov)>=.95
                # Sample oracle checks in addition to full baseline retained count.
                if counts['raw_records']%997==0:
                    old=p.coding_support(read,[(0,qlen)],blocks)
                    expected=bool(old and min(old['pid'],old['query_cds_coverage'],old['target_cds_coverage'])>=.95)
                    assert baseline==expected
                    counts['oracle_checks']+=1
                if baseline:passes.add('strict95')
                if pid<.95 or qcov<.90 or tcov<.90:continue
                f2=frame or compensated_indels(read,[(0,qlen)],blocks,12)
                s2=splice or compatible_junctions(observed,introns,lo,hi,3)
                if frame and splice:passes.add('strict90')
                for name,accepted in [('indel12_cov90',f2 and splice),('junction3_cov90',frame and s2),('combined_cov90',f2 and s2)]:
                    if accepted:
                        passes.add(name)
                        if qgene!=target:hit_examples[name]=dict(query=read.query_name,query_gene=qgene,target=target,target_transcript=tid,pid=pid,query_coverage=qcov,target_coverage=tcov,cigar=read.cigarstring)
                if len(passes)==len(SCENARIOS):break
            record_counts.update(passes)
            if 'strict95' not in passes:
                for name in passes:deltas[name].add((qgene,target))
            if 'strict90' not in passes:
                for name in hit_examples:
                    if len(examples[name])<30:examples[name].append(hit_examples[name])
            if counts['raw_records']%1000000==0:print(counts['raw_records'],flush=True)
    expected=json.loads((args.history/'graph.summary.json').read_text())['retained']
    assert record_counts['strict95']==expected,(record_counts,expected)
    for name,edges in deltas.items():
        with (args.output/(name+'.candidate_edges.tsv')).open('w') as f:
            for x,y in sorted(edges):f.write(x+'\t'+y+'\n')
    report=dict(counts=dict(counts),retained_records=dict(record_counts),candidate_edges={s:len(e) for s,e in deltas.items()},examples=examples,seconds=time.monotonic()-started,
                policy=dict(pid=.95,coverage_baseline=.95,coverage_test=.90,indel_window_bp=12,junction_equal_shift_tolerance_bp=3),baseline_retained_exact=True)
    (args.output/'filter.summary.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='examples'},indent=2),flush=True)

if __name__=='__main__':main()
