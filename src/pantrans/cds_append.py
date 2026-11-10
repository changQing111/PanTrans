"""Incremental CDS graph extension followed by fresh, last-only annotation."""
import argparse
import collections
import json
import time
from pathlib import Path

import networkx as nx
import pysam

from pantrans import cds_construct as pipeline
from pantrans.stage_timing import StageTimer
from pantrans.cds_merge import check_inputs, merge_prepared
from pantrans.cds_package import digest, load_package, write_package, mapping_policy
from pantrans.prepared_bundle import load_prepared_bundle, prepared_file


def add_history_edges(path, graph, historical_genes):
    count=0
    with Path(path).open() as stream:
        for line in stream:
            a,b=line.rstrip('\n').split('\t')
            if a not in historical_genes or b not in historical_genes:
                raise ValueError('Historical edge has an unknown endpoint')
            graph.add_edge(a,b);count+=1
    return count


def validate_legacy(history, prepared):
    """Bind legacy sequence inputs; coordinates/policy remain explicitly assumed."""
    for name,key in (('body.fa','gene'),('cds.fa','cds_valid')):
        if digest(history/name)!=digest(prepared_file(prepared,key)):
            raise ValueError('Legacy history sequence mismatch: '+name)
    raw=history/'cds.raw.bam'
    with pysam.AlignmentFile(raw) as bam, pysam.FastaFile(str(prepared_file(prepared,'gene'))) as fa:
        if dict(zip(bam.references,bam.lengths))!=dict(zip(fa.references,fa.lengths)):
            raise ValueError('Legacy raw BAM references differ from prepared gene bodies')
    return raw


def run(args):
    started=time.monotonic();timer=StageTimer(args.output)
    with timer.stage('validate_history'):
        checked=check_inputs(args.history_prepared,args.query_prepared)
        _,gene_sets=checked
        history=args.history_results.resolve()
        # Validate before creating output or starting a cross alignment.
        if args.import_legacy_history:
            if (history/'cds.graph.json').exists():
                raise ValueError('Packaged history must use normal validation, not legacy import')
            raw=validate_legacy(history,args.history_prepared)
            parent=dict(path=str(history),legacy=True,coordinates_and_old_mapping_assumed=True,
                        raw_bam_sha256=digest(raw),mapping_at_import=mapping_policy())
        else:
            package,edges=load_package(history,args.history_prepared)
            parent=dict(path=str(history/'cds.graph.json'),sha256=digest(history/'cds.graph.json'),
                        legacy=False,ancestors=package.get('parents',[]))
    args.output.mkdir(exist_ok=False)
    prepared=args.output/'merged_prepared'
    merged_hashes={}
    with timer.stage('merge_prepared'):
        merge_prepared(args.history_prepared,args.query_prepared,prepared,checked=checked,hashes=merged_hashes)
    with timer.stage('load_bundle'):
        regions,models,tids=load_prepared_bundle(prepared,args.output,graph_fasta_mode='link')
    graph=nx.DiGraph();graph.add_nodes_from(regions)
    stats=collections.Counter();stages={}
    with timer.stage('historical_edges'):
        if args.import_legacy_history:
            stages['history_refilter']=dict(pipeline.filter_cds_graph(
                raw,args.output/'history.refiltered.bam',
                {g:regions[g] for g in gene_sets[0]},
                {t:m for t,m in models.items() if m['gene'] in gene_sets[0]},tids,graph))
        else:
            add_history_edges(edges,graph,gene_sets[0])
    stats['historical_edges']=graph.number_of_edges()
    for name,query,target in (
        ('new_to_all',prepared_file(args.query_prepared,'cds_valid'),args.output/'body.fa'),
        ('history_to_new',prepared_file(args.history_prepared,'cds_valid'),prepared_file(args.query_prepared,'gene')),
    ):
        print(json.dumps(dict(stage=name,event='start')),flush=True)
        if not query.stat().st_size:
            stages[name]={'skipped_empty_cds':True};continue
        raw=args.output/(name+'.raw.bam');filtered=args.output/(name+'.filtered.bam')
        with timer.stage(name+'_mapping'):
            pipeline.minimap2_map(str(query),str(target),args.threads,str(raw))
        with timer.stage(name+'_filter'):
            counts=pipeline.filter_cds_graph(raw,filtered,regions,models,tids,graph)
        stages[name]=dict(counts);stats.update(counts)
    with timer.stage('clustering'):
        pre,last=pipeline.cluster_cds_graph(args.output,regions,graph,stats)
    with timer.stage('graph_package'):
        write_package(args.output,prepared,kind='append',parents=[parent],prepared_hashes=merged_hashes)
    summary=dict(operation='cds_append',workflow="cds-homology-v2",experimental=False,cross_alignments=stages,
                 historical_genes=len(gene_sets[0]),new_genes=len(gene_sets[1]),
                 annotation_strategy='reuse_transfer_evidence_then_rededuplicate_last',
                 direct_construct_equivalence_guaranteed=False,stages={})
    caches=list(getattr(args,'history_evidence_cache',[]) or [])
    default_cache=history/'last.evidence.sqlite'
    if default_cache.is_file() and default_cache not in caches:caches.append(default_cache)
    with timer.stage('last_annotation'):
        summary['stages']['last']=pipeline.annotation_stage('last',last,prepared,models,tids,args.output,args.workers,
            cache_sources=None if getattr(args,'no_evidence_cache',False) else caches,local_fasta=getattr(args,'local_fasta',False))
    summary['elapsed_seconds']=time.monotonic()-started
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
    return summary


def configure_parser(parser):
    parser.add_argument('--history-results',type=Path,required=True)
    parser.add_argument('--history-prepared',type=Path,required=True,
                        help='Original historical schema-2 inputs; for chained append use merged_prepared')
    parser.add_argument('--query-prepared',type=Path,required=True)
    parser.add_argument('-o','--output',type=Path,required=True)
    parser.add_argument('-t','--threads',type=int,default=32)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--import-legacy-history',action='store_true',
                        help='Explicitly import unpackaged CDS results; original coordinates/mapping policy cannot be proven')
    parser.add_argument('--history-evidence-cache',type=Path,nargs='+',default=[])
    parser.add_argument('--local-fasta',action='store_true')
    parser.add_argument('--no-evidence-cache',action='store_true')
    return parser


def validate_args(args, parser):
    if args.no_evidence_cache and args.history_evidence_cache:parser.error('Conflicting cache options')
    if args.workers<1 or args.threads<2:parser.error('workers >=1 and threads >=2 required')


def main(argv=None):
    parser=configure_parser(argparse.ArgumentParser(description=__doc__))
    args=parser.parse_args(argv)
    validate_args(args,parser)
    run(args)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
