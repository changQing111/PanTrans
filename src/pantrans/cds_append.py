"""Incremental CDS graph extension followed by fresh, last-only annotation."""
import argparse
import json
import time
from pathlib import Path

from pantrans import cds_construct as pipeline
from pantrans.stage_timing import StageTimer
from pantrans.cds_merge import check_inputs, merge_prepared
from pantrans.cds_package import digest, load_package, write_package
from pantrans.cds_target_graph import require_longest_inputs, graph_stage as per_variety_graph_stage
from pantrans.prepared_bundle import load_prepared_bundle


def run(args):
    started=time.monotonic();timer=StageTimer(args.output)
    with timer.stage('validate_history'):
        require_longest_inputs(args.history_prepared)
        require_longest_inputs(args.query_prepared)
        if args.import_legacy_history:
            raise ValueError('N2000 append does not import unattested legacy history; rebuild with current construct')
        checked=check_inputs(args.history_prepared,args.query_prepared)
        manifests,gene_sets=checked
        history=args.history_results.resolve()
        package,edges=load_package(history,args.history_prepared,longest=True)
        if not (history/'summary.json').is_file():
            raise ValueError('Historical annotation is incomplete: missing summary.json')
        parent=dict(path=str(history/'cds.graph.json'),sha256=digest(history/'cds.graph.json'),
                    legacy=False,ancestors=package.get('parents',[]))
    args.output.mkdir(exist_ok=False)
    prepared=args.output/'merged_prepared'
    merged_hashes={}
    with timer.stage('merge_prepared'):
        merge_prepared(args.history_prepared,args.query_prepared,prepared,checked=checked,hashes=merged_hashes)
    with timer.stage('load_bundle'):
        regions,models,tids=load_prepared_bundle(prepared,args.output,graph_fasta_mode='link')
    pre,last=per_variety_graph_stage(prepared,args.output,regions,models,tids,args.threads,timer,
        historical_edges=edges,historical_genes=gene_sets[0],
        historical_varieties=[v['variety'] for v in manifests[0]['varieties']])
    with timer.stage('graph_package'):
        write_package(args.output,prepared,kind='append',parents=[parent],prepared_hashes=merged_hashes,longest=True)
    mapping=json.loads((args.output/'graph.mapping.json').read_text())
    stages={t['variety']:dict(query_scope=t['query_scope'],queries=t['query_count'],
                            skipped_empty_cds=t.get('skipped_empty_cds',False),counts=t['counts'])
            for t in mapping['targets']}
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
    from pantrans.cds_outputs import export_last
    with timer.stage('export_validation'):
        summary['validation']=export_last(prepared,args.output,last)
    summary['graph_mapping']=mapping['mapping']
    summary['reused_variety_pairs']=mapping['reused_variety_pairs']
    summary['new_variety_pairs']=mapping['ordered_variety_pairs']
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
                        help='Unsupported by N2000 append; rebuild historical graph using current construct')
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
