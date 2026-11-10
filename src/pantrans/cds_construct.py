"""CDS homology graph construction and UTR-aware transcript annotation."""
import argparse, collections, csv, hashlib, json, os, shutil, subprocess, sys, tempfile, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from pantrans.annotation_runtime import FastaHandles, ordered_results, write_progress, local_annotation_inputs
from pantrans.annotation_evidence import EvidenceCache, cached_cluster, evidence_policy
from pantrans.stage_timing import StageTimer
from pantrans.cds_package import POLICY

import networkx as nx
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from pantrans.align_filter import minimap2_map
from pantrans.construct_di_graph import assign_sccs, get_conn_comp
from pantrans.transcript_redundancy import deduplicate_models
from pantrans.prepared_bundle import prepared_file, manifest, load_prepared_bundle


def overlap(blocks, intervals):
    return sum(max(0, min(b, y)-max(a, x)) for a,b in blocks for x,y in intervals)


def layout(read):
    q=0; r=read.reference_start; pieces=[]; junctions=[]; exons=[]; start=r; aligned=0; nmden=0; clips=0
    for op,n in read.cigartuples or []:
        if op in (0,7,8): pieces.append((q,q+n,r,r+n));q+=n;r+=n;aligned+=n;nmden+=n
        elif op==1:q+=n;aligned+=n;nmden+=n
        elif op==2:r+=n;nmden+=n
        elif op==3:exons.append((start,r));junctions.append((r,r+n));r+=n;start=r
        elif op in (4,5):q+=n;clips+=n
    if r>start:exons.append((start,r))
    return pieces,junctions,exons,aligned,nmden,clips


def coding_support(read, query_cds_blocks, target_cds_blocks):
    """Project query CDS query intervals to one target isoform; no strand synthesis."""
    if read.is_unmapped or read.is_supplementary or read.is_reverse or not read.has_tag('NM'):return None
    # A coding indel must not change the reading frame.
    qp=0;rp=read.reference_start
    for op,n in read.cigartuples or []:
        if op in (0,7,8):qp+=n;rp+=n
        elif op==1:
            if n%3 and any(max(qp,a)<min(qp+n,b) for a,b in query_cds_blocks):return None
            qp+=n
        elif op==2:
            if n%3 and any(max(rp,a)<min(rp+n,b) for a,b in target_cds_blocks):return None
            rp+=n
        elif op==3:rp+=n
        elif op in (4,5):qp+=n
    pieces,junctions,exons,aligned,den,clips=layout(read)
    if not den:return None
    projected=[];coding_segments=[]
    for qa,qb,ra,rb in pieces:
        for a,b in query_cds_blocks:
            lo,hi=max(qa,a),min(qb,b)
            if lo<hi:projected.append((ra+lo-qa,ra+hi-qa))
    qlen=sum(b-a for a,b in query_cds_blocks);tlen=sum(b-a for a,b in target_cds_blocks)
    if not qlen or not tlen:return None
    shared=overlap(projected,target_cds_blocks)
    # Only introns within coding spans must match. UTR junctions remain unrestricted.
    target_introns={(b,c) for (a,b),(c,d) in zip(target_cds_blocks,target_cds_blocks[1:])}
    qstart,qend=min(a for a,b in projected) if projected else 0,max(b for a,b in projected) if projected else 0
    coding_introns={j for j in junctions if j[0]>qstart and j[1]<qend}
    if any(j not in target_introns for j in coding_introns):return None
    # Reject a missing target coding intron spanned by a contiguous alignment block.
    if any(qstart<a and b<qend and (a,b) not in coding_introns for a,b in target_introns):return None
    for x,y in target_cds_blocks:
        segs=[(max(a,x),min(b,y)) for a,b in projected if max(a,x)<min(b,y)]
        if segs:coding_segments.extend(segs)
    return dict(query_cds_coverage=shared/qlen,target_cds_coverage=shared/tlen,
                pid=1-read.get_tag('NM')/den,junctions=junctions,exons=exons,
                projected=projected,coding_segments=coding_segments,aligned=aligned,clips=clips)


def graph_coding_support(read, query_length, target_cds_blocks):
    """Homology metrics only: splice differences and coding indels are allowed.

    Coverage retains the existing CDS-overlap definition; query bases mapping
    outside annotated target CDS do not count toward either CDS coverage.
    """
    if read.is_unmapped or read.is_supplementary or read.is_reverse or not read.has_tag('NM'):
        return None
    pieces,_,_,_,den,_=layout(read)
    target_length=sum(b-a for a,b in target_cds_blocks)
    if not den or query_length<=0 or target_length<=0:return None
    projected=[]
    for qa,qb,ra,rb in pieces:
        lo,hi=max(qa,0),min(qb,query_length)
        if lo<hi:projected.append((ra+lo-qa,ra+hi-qa))
    shared=overlap(projected,target_cds_blocks)
    return dict(pid=1-read.get_tag('NM')/den,
                query_cds_coverage=shared/query_length,
                target_cds_coverage=shared/target_length)


def writefa(path, records):
    with open(path,'w') as f:
        for name,seq in records:
            f.write('>'+name+'\n')
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')


def load_bundle(prepared, annotations, output, selected=None, graph_fasta_mode="copy"):
    if manifest(str(Path(prepared).resolve())).get('schema_version', 1) >= 2:
        return load_prepared_bundle(prepared, output, selected, graph_fasta_mode=graph_fasta_mode)
    if not annotations:
        raise ValueError('Legacy preprocessing requires --annotations; rerun preprocess for schema 2')
    regions={r['gene']:r for r in csv.DictReader(open(prepared_file(prepared, 'genes')),delimiter='\t') if selected is None or r['gene'] in selected}
    exons=collections.defaultdict(list);tidgene={}
    for r in csv.DictReader(open(prepared_file(prepared, 'exons')),delimiter='\t'):
        if r['gene'] in regions:exons[r['transcript']].append((int(r['gdna_start0']),int(r['gdna_end0'])));tidgene[r['transcript']]=r['gene']
    originals={}
    for r in csv.DictReader(open(prepared_file(prepared, 'ids')),delimiter='\t'):
        if r['feature']=='transcript' and r['output_gene'] in regions:originals[r['variety'],r['original_id']]=r['output_id']
    cds=collections.defaultdict(list)
    for ap in annotations:
        variety=Path(ap).name.split('.')[0]
        for line in open(ap):
            if line.startswith('#'):continue
            f=line.rstrip().split('\t')
            if len(f)!=9 or f[2]!='CDS':continue
            a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
            for orig in a.get('Parent','').split(','):
                tid=originals.get((variety,orig))
                if tid is None:continue
                g=tidgene[tid];r=regions[g];left,right=int(r['extract_start0']),int(r['extract_end0']);x,y=int(f[3])-1,int(f[4])
                lo,hi=(right-y,right-x) if r['strand']=='-' else (x-left,y-left)
                cds[tid].append((lo,hi,int(f[7]) if f[7] in ['0','1','2'] else -1))
    models={};gene_tids=collections.defaultdict(list);status=collections.Counter();fasta=pysam.FastaFile(str(prepared_file(prepared, 'gdna')))
    cdna=pysam.FastaFile(str(prepared_file(prepared, 'cdna')))
    with open(output/'cds_status.tsv','w') as audit,open(output/'cds.fa','w') as cfa,open(output/'body.fa','w') as bfa:
        audit.write('transcript\tgene\tstatus\tcds_length\n')
        for g,r in regions.items():
            s=fasta.fetch(g).upper();a,b=int(r['body_start0']),int(r['body_end0']);bfa.write('>'+g+'\n'+s[a:b]+'\n')
        for tid,blocks in exons.items():
            g=tidgene[tid];r=regions[g];s=fasta.fetch(g).upper();blocks=sorted(set(blocks));c=sorted(set(cds.get(tid,[])));coding=''.join(s[a:b] for a,b,p in c);qblocks=[];offset=0
            for a,b in blocks:
                for x,y,p in c:
                    if max(a,x)<min(b,y):qblocks.append((offset+max(a,x)-a,offset+min(b,y)-a))
                offset+=b-a
            if ''.join(s[a:b] for a,b in blocks)!=cdna.fetch(tid).upper():raise ValueError('Input cDNA/exon mismatch '+tid)
            reason='valid';cumulative=0
            for a,b,p in c:
                if p!=(3-cumulative%3)%3:reason='invalid_phase'
                cumulative+=b-a
            if not c:reason='no_CDS'
            elif not coding.startswith('ATG') or coding[-3:] not in ['TAA','TAG','TGA'] or len(coding)%3:reason='partial_or_invalid_ORF'
            elif '*' in str(Seq(coding).translate())[:-1]:reason='internal_stop'
            elif sum(b-a for a,b in qblocks)!=len(coding) or qblocks[-1][1]-qblocks[0][0]!=len(coding):reason='CDS_not_contiguous_in_cDNA'
            status[reason]+=1;audit.write(f'{tid}\t{g}\t{reason}\t{len(coding)}\n')
            models[tid]=dict(gene=g,exons=blocks,cds=[(a,b) for a,b,p in c],qcds=qblocks,valid=reason=='valid',cds_len=len(coding),cdna_len=offset)
            gene_tids[g].append(tid)
            if reason=='valid':cfa.write('>'+tid+'\n'+coding+'\n')
    fasta.close();cdna.close()
    (output/'cds_status.json').write_text(json.dumps(status,indent=2))
    return regions,models,gene_tids


def filter_cds_graph(raw_bam,filtered_bam,regions,models,gene_tids,G):
    stats=collections.Counter();target_cds={}
    with pysam.AlignmentFile(raw_bam) as bam,pysam.AlignmentFile(filtered_bam,'wb',header=bam.header) as filt:
        for read in bam:
            stats['records']+=1
            if read.is_unmapped or read.is_supplementary:stats['unmapped_supplementary']+=1;continue
            if read.is_reverse:stats['reverse_rejected']+=1;continue
            qm=models[read.query_name];target=read.reference_name;passed=False
            if target not in target_cds:
                shift=int(regions[target]['body_start0'])
                target_cds[target]=tuple(dict.fromkeys(
                    tuple((a-shift,b-shift) for a,b in models[t]['cds'])
                    for t in gene_tids[target] if models[t]['valid']))
            for blocks in target_cds[target]:
                hit=graph_coding_support(read,qm['cds_len'],blocks)
                if hit and hit['pid']>=POLICY['pid'] and hit['query_cds_coverage']>=POLICY['query_cds_coverage'] and hit['target_cds_coverage']>=POLICY['target_cds_coverage']:
                    passed=True;break
            if passed:filt.write(read);G.add_edge(qm['gene'],target);stats['retained']+=1
            else:stats['coding_support_rejected']+=1
    return stats


def cluster_cds_graph(output,regions,G,stats):
    with open(output/'graph.edges.tsv','w') as f:
        for a,b in sorted(G.edges()):f.write(a+'\t'+b+'\n')
    bed={g:[r['chrom'],int(r['gene_start0']),int(r['gene_end0'])] for g,r in regions.items()};lengths={g:int(r['body_end0'])-int(r['body_start0']) for g,r in regions.items()};varieties={g:r['variety'] for g,r in regions.items()}
    pre,last=assign_sccs(get_conn_comp(G),G,lengths,bed,list(dict.fromkeys(varieties.values())),'CS',main_chroms=['Chr'+str(i)+s for i in range(1,8) for s in 'ABD'],gene_variety_dic=varieties)
    for stage,clusters in [('pre',pre),('last',last)]:
        with open(output/(stage+'.cluster'),'w') as f:
            for c in clusters:f.write('\t'.join(c)+'\n')
    stats.update(nodes=G.number_of_nodes(),edges=G.number_of_edges(),pre_clusters=len(pre),last_clusters=len(last));(output/'graph.summary.json').write_text(json.dumps(stats,indent=2))
    return pre,last


def graph_stage(output,regions,models,gene_tids,threads,timer=None):
    timer=timer or StageTimer(output)
    with timer.stage('cds_mapping'):
        minimap2_map(str(output/'cds.fa'),str(output/'body.fa'),threads,str(output/'cds.raw.bam'))
    G=nx.DiGraph();G.add_nodes_from(regions)
    with timer.stage('cds_filter'):
        stats=filter_cds_graph(output/'cds.raw.bam',output/'cds.filtered.bam',regions,models,gene_tids,G)
    with timer.stage('clustering'):
        return cluster_cds_graph(output,regions,G,stats)


def merge_intervals(blocks):
    out=[]
    for a,b in sorted(blocks):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    return out


def annotation_cluster(cluster,prepared,models,gene_tids,work, fasta_handles=None, query_ids=None):
    hub=cluster[0];candidates=[];audit=[];accepted_sam=[]
    # Original representative models are deterministic seeds, explicitly labeled.
    for tid in gene_tids[hub]:
        m=models[tid];candidates.append(dict(source=tid,exons=m['exons'],cds=m['cds'] if m['valid'] else [],complete=True,origin='original_annotation',pid=1.0))
    queries=[tid for g in cluster[1:] for tid in gene_tids[g] if models[tid]['valid']]
    if query_ids is not None:
        requested=set(query_ids);queries=[t for t in queries if t in requested]
    if not queries:return hub,candidates,audit,accepted_sam
    with tempfile.TemporaryDirectory(prefix='cdna-',dir=work) as temp:
        temp=Path(temp)
        with (nullcontext(fasta_handles.get(prepared_file(prepared, 'cdna'))) if fasta_handles else pysam.FastaFile(str(prepared_file(prepared, 'cdna')))) as fa:writefa(temp/'query.fa',((t,fa.fetch(t)) for t in queries))
        with (nullcontext(fasta_handles.get(prepared_file(prepared, 'gdna'))) if fasta_handles else pysam.FastaFile(str(prepared_file(prepared, 'gdna')))) as fa:target_sequence=fa.fetch(hub).upper();writefa(temp/'target.fa',[(hub,target_sequence)])
        # One representative per invocation: no cross-cluster target competition.
        cmd=['minimap2','-ax','splice:hq','-uf','--secondary=yes','-N','100','-t','1',str(temp/'target.fa'),str(temp/'query.fa')]
        with open(temp/'raw.sam','w') as sam,open(temp/'minimap.log','w') as log:subprocess.run(cmd,stdout=sam,stderr=log,check=True)
        best={}
        with pysam.AlignmentFile(temp/'raw.sam') as bam:
            for read in bam:
                if read.is_unmapped or read.is_supplementary or read.is_reverse:continue
                qm=models[read.query_name]
                for target_tid in gene_tids[hub]:
                    tm=models[target_tid]
                    if not tm['valid']:continue
                    hit=coding_support(read,qm['qcds'],tm['cds'])
                    # Full query coding sequence must map; UTR may be partial and explicitly marked.
                    if not hit or hit['pid']<.95 or hit['query_cds_coverage']<.999999 or hit['target_cds_coverage']<.95:continue
                    projected=merge_intervals(hit['projected'])
                    if sum(b-a for a,b in projected)!=qm['cds_len']:continue
                    projected_sequence=''.join(target_sequence[a:b] for a,b in projected)
                    if len(projected_sequence)%3 or not projected_sequence.startswith('ATG') or projected_sequence[-3:] not in ['TAA','TAG','TGA'] or '*' in str(Seq(projected_sequence).translate())[:-1]:continue
                    model=dict(source=read.query_name,exons=hit['exons'],cds=projected,complete=hit['clips']==0,origin='cross_variety_alignment',pid=hit['pid'])
                    # Accepted records are forward on the gene-oriented target.
                    cigar=read.cigartuples or []
                    model['missing_5prime']=bool(cigar and cigar[0][0] in (4,5))
                    model['missing_3prime']=bool(cigar and cigar[-1][0] in (4,5))
                    rank=(model['complete'],hit['aligned'],hit['pid'])
                    if read.query_name not in best or rank>best[read.query_name][0]:best[read.query_name]=(rank,model,read.to_string())
        candidates.extend(best[t][1] for t in queries if t in best);accepted_sam.extend(best[t][2] for t in queries if t in best)
        for tid in queries:audit.append([hub,tid,'complete' if tid in best and best[tid][1]['complete'] else 'partial' if tid in best else 'no_accepted_model'])
    return hub,candidates,audit,accepted_sam


def deduplicate(candidates):
    return deduplicate_models(candidates)


def emit_model(out,hub,idx,m):
    tid=f'{hub}.{idx}';exons=m['exons'];attrs=f'ID={tid};Parent={hub};completeness={"complete" if m["complete"] else "partial"};model_origin={m["origin"]};source_transcript={m["source"]}'
    attrs+=f';utr_evidence={m.get("utr_evidence","unknown")};redundancy_status={m.get("redundancy_status","unknown")}'
    out.write(f'{hub}\tPanTransCDS\tmRNA\t{exons[0][0]+1}\t{exons[-1][1]}\t.\t+\t.\t{attrs}\n')
    for i,(a,b) in enumerate(exons,1):out.write(f'{hub}\tPanTransCDS\texon\t{a+1}\t{b}\t.\t+\t.\tID={tid}.exon{i};Parent={tid}\n')
    used=0
    for a,b in m['cds']:
        phase=(3-used%3)%3;out.write(f'{hub}\tPanTransCDS\tCDS\t{a+1}\t{b}\t.\t+\t{phase}\tParent={tid}\n');used+=b-a
    if m['cds']:
        cs,ce=m['cds'][0][0],m['cds'][-1][1]
        for a,b in exons:
            for label,x,y in [('five_prime_UTR',a,min(b,cs)),('three_prime_UTR',max(a,ce),b)]:
                if x<y:out.write(f'{hub}\tPanTransCDS\t{label}\t{x+1}\t{y}\t.\t+\t.\tParent={tid}\n')


def annotation_stage(stage,clusters,prepared,models,gene_tids,output,workers,cache_sources=None,local_fasta=False):
    context=local_annotation_inputs(prepared,output,stage) if local_fasta else nullcontext(prepared)
    with context as annotation_prepared:
        return _annotation_stage(stage,clusters,annotation_prepared,models,gene_tids,output,workers,cache_sources)


def _annotation_stage(stage,clusters,prepared,models,gene_tids,output,workers,cache_sources=None):
    stats=collections.Counter();work=Path(os.environ.get('PANTRANS_LOCAL_SCRATCH', tempfile.gettempdir()));started=time.monotonic()
    write_progress(output,stage,0,len(clusters),started)
    policy=evidence_policy([annotation_cluster,coding_support,layout,merge_intervals,writefa]) if cache_sources is not None else None
    cache_context=EvidenceCache(output/(stage+'.evidence.sqlite'),cache_sources) if cache_sources is not None else nullcontext(None)
    with cache_context as evidence, open(output/(stage+'.gff3'),'w') as gff,open(output/(stage+'.cdna.fa'),'w') as cdna,open(output/(stage+'.model_sources.tsv'),'w') as sources,open(output/(stage+'.alignment_status.tsv'),'w') as audit,pysam.FastaFile(str(prepared_file(prepared, 'gdna'))) as fa:
        gff.write('##gff-version 3\n');sources.write('gene\ttranscript\tcompleteness\tsources\tutr_evidence\tredundancy_status\tcompatible_sources\tsource_relations_json\n');audit.write('target\tquery\tstatus\n')
        bam_header=pysam.AlignmentHeader.from_dict({'HD':{'VN':'1.6'},'SQ':[{'SN':c[0],'LN':fa.get_reference_length(c[0])} for c in clusters]})
        with FastaHandles() as fasta_handles, pysam.AlignmentFile(output/(stage+'.annotation.bam'),'wb',header=bam_header) as accepted_bam, ThreadPoolExecutor(max_workers=workers) as pool:
            if evidence is None:
                fn=lambda c: annotation_cluster(c,prepared,models,gene_tids,work,fasta_handles=fasta_handles)
            else:
                fn=lambda c: cached_cluster(c,prepared,models,gene_tids,work,fasta_handles,evidence,policy,annotation_cluster)
            for index,result in enumerate(ordered_results(pool,fn,clusters,workers*4),1):
                hub,candidates,checks,sam_records=result;unique=deduplicate(candidates);stats['genes']+=bool(unique);stats['transcripts']+=len(unique)
                for sam in sam_records:accepted_bam.write(pysam.AlignedSegment.fromstring(sam,bam_header))
                for row in checks:audit.write('\t'.join(row)+'\n');stats['alignment_'+row[-1]]+=1
                if not unique:
                    if index % 1000 == 0:write_progress(output,stage,index,len(clusters),started)
                    continue
                first=min(m['exons'][0][0] for m in unique);last=max(m['exons'][-1][1] for m in unique);gff.write(f'{hub}\tPanTransCDS\tgene\t{first+1}\t{last}\t.\t+\t.\tID={hub}\n');seq=fa.fetch(hub)
                for i,m in enumerate(unique,1):
                    emit_model(gff,hub,i,m);state='complete' if m['complete'] else 'partial';stats[state]+=1
                    sources.write(f'{hub}\t{hub}.{i}\t{state}\t'+','.join(m['sources'])+'\t'+m['utr_evidence']+'\t'+m['redundancy_status']+'\t'+','.join(m['compatible_sources'])+'\t'+json.dumps(m['source_relations'],sort_keys=True)+'\n');cdna.write(f'>{hub}.{i}\n'+''.join(seq[a:b] for a,b in m['exons'])+'\n')
                if index % 1000 == 0:
                    for handle in (gff,cdna,sources,audit):handle.flush()
                    write_progress(output,stage,index,len(clusters),started)
    if evidence is not None:
        (output/(stage+'.cache_stats.json')).write_text(json.dumps(dict(evidence.counts),indent=2))
    write_progress(output,stage,len(clusters),len(clusters),started)
    (output/(stage+'.summary.json')).write_text(json.dumps(stats,indent=2));return dict(stats)


def reuse_graph(source, output, regions):
    checks={};clusters=[]
    for stage in ('pre','last'):
        path=source/(stage+'.cluster')
        rows=[line.split() for line in path.read_text().splitlines() if line.strip()]
        members=[g for row in rows for g in row]
        if len(members)!=len(set(members)) or set(members)!=set(regions):
            raise ValueError('Reused '+stage+' clusters must partition the prepared genes exactly')
        clusters.append(rows)
    for name in ('cds.raw.bam','cds.filtered.bam','graph.edges.tsv','pre.cluster','last.cluster','graph.summary.json'):
        path=(source/name).resolve()
        if not path.is_file():raise ValueError('Missing graph input: '+str(path))
        st=path.stat();checks[name]=dict(path=str(path),size=st.st_size,mtime_ns=st.st_mtime_ns)
    for name in ('pre.cluster','last.cluster','graph.summary.json'):
        shutil.copyfile(source/name,output/name)
    (output/'reused_graph.json').write_text(json.dumps(checks,indent=2))
    return clusters


def configure_parser(p):
    p.add_argument('--prepared',type=Path,required=True,help='preprocess output directory (schema 2 recommended)')
    p.add_argument('--annotations',type=Path,help='GFF3 list required only for legacy schema-1 preprocessing')
    p.add_argument('-o','--output',type=Path,required=True)
    p.add_argument('--genes',type=Path)
    p.add_argument('-t','--threads',type=int,default=32)
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--annotate-pre',action='store_true')
    p.add_argument('--reuse-graph',type=Path,help='Annotation-only restart from completed graph/cluster files')
    p.add_argument('--history-evidence-cache',type=Path,nargs='+',default=[])
    p.add_argument('--local-fasta',action='store_true',help='Stage annotation FASTAs on node-local scratch')
    p.add_argument('--no-evidence-cache',action='store_true')
    return p


def validate_args(args, parser):
    if args.workers<1 or args.threads<2:parser.error('workers must be >=1 and threads >=2')
    if args.reuse_graph and args.genes:parser.error('--reuse-graph cannot be combined with --genes')
    if args.no_evidence_cache and args.history_evidence_cache:parser.error('Cache inputs conflict with --no-evidence-cache')


def run(args):
    args.output.mkdir(exist_ok=False);started=time.monotonic();timer=StageTimer(args.output)
    selected=set(args.genes.read_text().split()) if args.genes else None
    annotations=args.annotations.read_text().splitlines() if args.annotations else []
    with timer.stage('load_bundle'):
        regions,models,gt=load_bundle(args.prepared,annotations,args.output,selected,graph_fasta_mode='link')
    reused_policy = None
    if args.reuse_graph:
        if (args.reuse_graph/'cds.graph.json').is_file():
            from pantrans.cds_package import load_package
            reused_package,_=load_package(args.reuse_graph,args.prepared)
            reused_policy=reused_package['policy']
        with timer.stage('reuse_graph'):pre,last=reuse_graph(args.reuse_graph,args.output,regions)
    else:
        pre,last=graph_stage(args.output,regions,models,gt,args.threads,timer=timer)
    if not args.genes and not args.reuse_graph and manifest(str(args.prepared.resolve())).get('schema_version')==2:
        from pantrans.cds_package import write_package
        with timer.stage('graph_package'):write_package(args.output,args.prepared,kind='construct')
    summaries={}
    for stage,clusters in ([('pre',pre)] if args.annotate_pre else [])+[('last',last)]:
        with timer.stage(stage+'_annotation'):
            summaries[stage]=annotation_stage(stage,clusters,args.prepared,models,gt,args.output,args.workers,
                cache_sources=None if args.no_evidence_cache else args.history_evidence_cache,local_fasta=args.local_fasta)
    (args.output/'summary.json').write_text(json.dumps(dict(elapsed_seconds=time.monotonic()-started,
        stages=summaries,workflow="cds-homology-v2",experimental=False,
        graph_thresholds=({k:POLICY[k] for k in ('pid','query_cds_coverage','target_cds_coverage')}
                          if not args.reuse_graph or reused_policy else None),
        graph_structure_checks=False if not args.reuse_graph or reused_policy else None,
        graph_origin='reused' if args.reuse_graph else 'fresh',
        graph_policy_verified=not args.reuse_graph or reused_policy is not None,
        original_representative_models_seeded=True),indent=2))
    print(json.dumps(summaries),flush=True)

def main(argv=None):
    parser=configure_parser(argparse.ArgumentParser(description=__doc__))
    args=parser.parse_args(argv)
    validate_args(args,parser)
    run(args)
    return 0


if __name__=='__main__':
    sys.exit(main())
