"""Merge disjoint schema-2 prepared datasets without losing original transcripts."""
import csv
import json
import shutil
import hashlib
from pathlib import Path

import pysam
from .prepared_bundle import manifest, prepared_file

FASTA_KEYS = ('gene', 'cds_valid', 'cds', 'gdna', 'cdna')
TABLE_KEYS = ('genes', 'exons', 'cds_regions', 'cds_status', 'bed', 'ids')


def check_inputs(history, query):
    manifests = [manifest(str(Path(p).resolve())) for p in (history, query)]
    for data in manifests:
        if data.get('schema_version') != 2:
            raise ValueError('Both prepared datasets must use schema 2')
    for key in ('flank', 'namespace_ids', 'cds_validation', 'gdna_orientation',
                'coordinate_system', 'graph_query', 'graph_target'):
        if manifests[0].get(key) != manifests[1].get(key):
            raise ValueError('Incompatible preprocessing policy: '+key)
    names=[]
    for data in manifests:
        values=[v['variety'] for v in data['varieties']]
        if not values or len(values)!=len(set(values)):
            raise ValueError('Missing or duplicate variety names')
        names.append(set(values))
    if names[0] & names[1]:
        raise ValueError('New varieties overlap historical varieties')
    gene_sets=[];transcript_sets=[]
    for directory, allowed in zip((history,query),names):
        genes=set();transcripts=set()
        with prepared_file(directory,'genes').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['gene'] in genes or row['variety'] not in allowed:
                    raise ValueError('Duplicate gene or inconsistent variety in gene regions')
                genes.add(row['gene'])
        with prepared_file(directory,'models').open() as stream:
            for line in stream:
                model=json.loads(line)
                if model['transcript'] in transcripts or model['gene'] not in genes:
                    raise ValueError('Duplicate transcript or unknown model gene')
                transcripts.add(model['transcript'])
        if not genes:
            raise ValueError('Empty prepared gene set')
        gene_sets.append(genes);transcript_sets.append(transcripts)
    if gene_sets[0] & gene_sets[1] or transcript_sets[0] & transcript_sets[1]:
        raise ValueError('Historical and query gene/transcript IDs overlap')
    return manifests,gene_sets


def merge_prepared(history, query, output, checked=None, hashes=None):
    manifests,gene_sets=checked if checked is not None else check_inputs(history,query)
    output=Path(output);output.mkdir(exist_ok=False)
    files={}
    # Preserve the optional per-gene query artifacts when both inputs have them.
    # Older schema-2 bundles remain valid, without advertising absent outputs.
    with_longest = all(all(key in data['files'] for key in ('cds_longest', 'cds_longest_selection'))
                       for data in manifests)
    fasta_keys = FASTA_KEYS + (('cds_longest',) if with_longest else ())
    table_keys = TABLE_KEYS + (('cds_longest_selection',) if with_longest else ())
    transcript_ids=set();valid_ids=set();coding_ids=set()
    longest = {}
    for directory in (history,query):
        with prepared_file(directory,'models').open() as stream:
            for line in stream:
                model=json.loads(line);transcript_ids.add(model['transcript'])
                if model['valid']:
                    valid_ids.add(model['transcript'])
                    if with_longest:
                        candidate = (-model['cds_len'], model['transcript'])
                        gene = model['gene']
                        if gene not in longest or candidate < longest[gene]:
                            longest[gene] = candidate
                if model['cds_len']:coding_ids.add(model['transcript'])
    expected={'gene':gene_sets[0]|gene_sets[1], 'gdna':gene_sets[0]|gene_sets[1],
              'cdna':transcript_ids, 'cds_valid':valid_ids, 'cds':coding_ids}
    if with_longest:
        expected['cds_longest'] = {candidate[1] for candidate in longest.values()}
    for key in fasta_keys:
        name='merged_'+key+'.fa';files[key]=name
        seen=set();hasher=hashlib.sha256()
        with (output/name).open('wb') as dst:
            for directory in (history,query):
                last=b''
                with prepared_file(directory,key).open('rb') as src:
                    while True:
                        batch=src.readlines(8*1024*1024)
                        if not batch:break
                        block=b''.join(batch);dst.write(block);hasher.update(block);last=batch[-1]
                        for line in batch:
                            if line.startswith(b'>'):
                                identifier=line[1:].split()[0].decode()
                                if identifier in seen:
                                    raise ValueError('Duplicate FASTA identifier: '+identifier)
                                seen.add(identifier)
                if last and not last.endswith(b'\n'):
                    dst.write(b'\n');hasher.update(b'\n')
        if hashes is not None:hashes[key]=hasher.hexdigest()
        if seen!=expected[key]:
            raise ValueError('FASTA identifiers do not match prepared models: '+key)
        if (output/name).stat().st_size:
            pysam.faidx(str(output/name))
    for key in table_keys+('models',):
        name='merged_'+key+('.jsonl' if key=='models' else '.tsv');files[key]=name
        hasher=hashlib.sha256()
        with (output/name).open('wb') as dst:
            header=None
            for directory in (history,query):
                with prepared_file(directory,key).open('rb') as src:
                    if key not in ('bed','models'):
                        current=src.readline()
                        if header is None:
                            header=current;dst.write(current);hasher.update(current)
                        elif current!=header:
                            raise ValueError('Incompatible table header: '+key)
                    last=b''
                    for block in iter(lambda:src.read(8*1024*1024),b''):
                        dst.write(block);hasher.update(block);last=block
                    if last and not last.endswith(b'\n'):
                        dst.write(b'\n');hasher.update(b'\n')
        if hashes is not None:hashes[key]=hasher.hexdigest()
    data=dict(manifests[0]);data['files']=files
    data['varieties']=manifests[0]['varieties']+manifests[1]['varieties']
    data['empty_fastas']=[key for key in fasta_keys if not (output/files[key]).stat().st_size]
    if with_longest:
        data['longest_cds_selection'] = dict(manifests[0]['longest_cds_selection'])
        for key in ('genes_with_valid_cds', 'genes_without_valid_cds', 'genes_with_max_cds_length_ties'):
            data['longest_cds_selection'][key] = sum(v[key] for v in data['varieties'])
    else:
        data.pop('longest_cds_selection', None)
    data['files']['manifest']='merged_preprocess.json'
    (output/'merged_preprocess.json').write_text(json.dumps(data,indent=2)+'\n')
    return gene_sets
