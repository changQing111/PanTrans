"""Disk-backed comparison of CDS append and direct construct (Linux sort)."""
import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

from Bio import SeqIO
from pantrans.prepared_bundle import prepared_file


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def clusters(path, representatives):
    with path.open() as stream:
        for line in stream:
            members=line.split()
            yield encoded([members[0],sorted(members)] if representatives else sorted(members))


def models(path):
    """Canonicalize each gene's model multiset, independent of transcript numbering."""
    gene=None;transcripts={}
    with path.open() as stream:
        for line in stream:
            if line.startswith('#'):continue
            f=line.rstrip('\n').split('\t')
            if len(f)!=9:continue
            attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
            if f[2]=='gene':
                if gene is not None:yield encoded([gene,sorted(transcripts.values(),key=encoded)])
                gene=f[0];transcripts={}
            elif f[2]=='mRNA':
                transcripts[attrs['ID']]=dict(completeness=attrs.get('completeness'),
                    utr_evidence=attrs.get('utr_evidence'),redundancy_status=attrs.get('redundancy_status'),features=[])
            else:
                for parent in attrs.get('Parent','').split(','):
                    if parent in transcripts:
                        transcripts[parent]['features'].append([f[2],int(f[3]),int(f[4]),f[6],f[7]])
    if gene is not None:yield encoded([gene,sorted(transcripts.values(),key=encoded)])


def sequences(path, by_gene=False):
    with path.open() as stream:
        for record in SeqIO.parse(stream,'fasta'):
            name=record.id.rsplit('.',1)[0] if by_gene else record.id
            yield encoded([name,hashlib.sha256(str(record.seq).upper().encode()).hexdigest()])


def model_inputs(path):
    with path.open() as stream:
        for line in stream:yield encoded(json.loads(line))


def lines(path, header=False):
    with path.open() as stream:
        if header:next(stream,None)
        for line in stream:yield line.rstrip('\n')


def compare(left, right, output, name):
    """Multiset difference: retain duplicate sequences/models rather than collapsing."""
    counts=dict(append_records=0,construct_records=0,shared=0,append_only=0,construct_only=0)
    with tempfile.TemporaryDirectory(dir=output) as td:
        paths=[]
        for side,rows in [('append',left),('construct',right)]:
            raw=Path(td)/(side+'.raw');ordered=Path(td)/(side+'.sorted')
            with raw.open('w') as stream:
                for row in rows:
                    stream.write(row+'\n');counts[side+'_records']+=1
            subprocess.run(['sort','-S','256M','-T',td,'-o',str(ordered),str(raw)],
                           check=True,env=dict(os.environ,LC_ALL='C'))
            paths.append(ordered)
        with paths[0].open() as a,paths[1].open() as b, (output/(name+'.append_only.txt')).open('w') as ao,(output/(name+'.construct_only.txt')).open('w') as bo:
            x=next(a,None);y=next(b,None)
            while x is not None or y is not None:
                if y is None or (x is not None and x<y):
                    counts['append_only']+=1;ao.write(x);x=next(a,None)
                elif x is None or y<x:
                    counts['construct_only']+=1;bo.write(y);y=next(b,None)
                else:
                    counts['shared']+=1;x=next(a,None);y=next(b,None)
    counts['equal']=not(counts['append_only'] or counts['construct_only'])
    return counts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--append',type=Path,required=True);p.add_argument('--construct',type=Path,required=True)
    p.add_argument('--construct-prepared',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    # Success summaries are mandatory; never compare partial outputs.
    summaries={key:json.loads((directory/'summary.json').read_text()) for key,directory in [('append',a.append),('construct',a.construct)]}
    for key,summary in summaries.items():
        if 'last' not in summary['stages']:raise ValueError('Missing completed last annotation: '+key)
    a.output.mkdir(exist_ok=False)
    report=dict(paths={key:str(value.resolve()) for key,value in vars(a).items()},runs=summaries,comparisons={})
    jobs=[('graph_edges',lines(a.append/'graph.edges.tsv'),lines(a.construct/'graph.edges.tsv'))]
    for stage in ('pre','last'):
        for reps in (False,True):
            jobs.append((stage+('_representatives_members' if reps else '_partition'),clusters(a.append/(stage+'.cluster'),reps),clusters(a.construct/(stage+'.cluster'),reps)))
    jobs.extend([
        ('last_models',models(a.append/'last.gff3'),models(a.construct/'last.gff3')),
        ('last_cdna_by_hub',sequences(a.append/'last.cdna.fa',True),sequences(a.construct/'last.cdna.fa',True)),
        ('last_alignment_status',lines(a.append/'last.alignment_status.tsv',True),lines(a.construct/'last.alignment_status.tsv',True)),
        ('last_model_sources_exact_ids',lines(a.append/'last.model_sources.tsv',True),lines(a.construct/'last.model_sources.tsv',True)),
    ])
    for key in ('gene','cds_valid','cdna','gdna'):
        jobs.append(('input_'+key,sequences(prepared_file(a.append/'merged_prepared',key)),sequences(prepared_file(a.construct_prepared,key))))
    jobs.append(('input_models',model_inputs(prepared_file(a.append/'merged_prepared','models')),model_inputs(prepared_file(a.construct_prepared,'models'))))
    for name,left,right in jobs:
        report['comparisons'][name]=compare(left,right,a.output,name)
        print(name,report['comparisons'][name],flush=True)
        (a.output/'comparison.progress.json').write_text(json.dumps(report,indent=2))
    report['all_equal']=all(x['equal'] for x in report['comparisons'].values())
    (a.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
