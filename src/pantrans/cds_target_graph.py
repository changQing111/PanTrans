"""Longest-valid-CDS queries against complete per-variety gene-body targets."""
import collections
import csv
import hashlib
import json
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser
import networkx as nx
import pysam

from .align_filter import minimap2_map
from .cds_package import CONSTRUCT_OPTIONS, mapping_policy
from .prepared_bundle import manifest, prepared_file


def save(path, data):
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def require_longest_inputs(prepared):
    data = manifest(str(Path(prepared).resolve()))
    keys = ('cds_longest', 'cds_longest_selection')
    if data.get('schema_version') != 2 or any(k not in data.get('files', {}) for k in keys):
        raise ValueError('CDS graph construction/append requires schema-2 inputs with longest valid CDS; rerun pantrans preprocess into a new directory')
    if any(not prepared_file(prepared, k).is_file() for k in keys):
        raise ValueError('Missing longest valid CDS inputs; rerun pantrans preprocess')
    return data


def prepare_queries(prepared, output, regions, models, gene_tids, subset=False, allow_empty=False):
    """Verify deterministic selection and sequence identity without loading FASTAs."""
    chosen = {}
    for g in regions:
        valid = [t for t in gene_tids.get(g, ()) if models[t]['valid']]
        if valid:
            chosen[g] = min(valid, key=lambda t: (-models[t]['cds_len'], t))
    queries = {t: g for g, t in chosen.items()}
    if not queries and not allow_empty:
        raise ValueError('No valid CDS for graph mapping')
    seen = set()
    with prepared_file(prepared, 'cds_longest_selection').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            g = row['gene']
            if g not in regions:
                if subset:continue
                raise ValueError('Unknown gene in longest CDS selection: ' + g)
            if g in seen:raise ValueError('Duplicate gene in longest CDS selection: ' + g)
            seen.add(g)
            if g in chosen:
                tid = chosen[g]
                if row['status'] != 'selected' or row['transcript'] != tid or int(row['CDS_length']) != models[tid]['cds_len']:
                    raise ValueError('Longest CDS selection disagrees with models: ' + g)
            elif row['transcript'] or int(row['CDS_length']) or row['status'] not in ('no_valid_CDS', 'no_transcripts'):
                raise ValueError('Invalid longest CDS selection status: ' + g)
    if seen != set(regions):raise ValueError('Longest CDS selection must cover every gene')
    expected = {}
    with prepared_file(prepared, 'cds_valid').open() as f:
        for title, seq in SimpleFastaParser(f):
            tid = title.split()[0]
            if tid in queries:
                if tid in expected:raise ValueError('Duplicate selected CDS: ' + tid)
                expected[tid] = hashlib.sha256(seq.upper().encode()).digest()
    if set(expected) != set(queries):raise ValueError('Selected CDS missing from all-valid CDS FASTA')
    seen = set()
    # Write a verified, self-contained query snapshot; preserve transcript IDs.
    with prepared_file(prepared, 'cds_longest').open() as f, (output/'longest.cds.fa').open('w') as out:
        for title, seq in SimpleFastaParser(f):
            tid = title.split()[0]
            if subset and tid not in queries:continue
            if tid not in queries or tid in seen:raise ValueError('Unexpected/duplicate longest CDS: ' + tid)
            if len(seq) != models[tid]['cds_len'] or hashlib.sha256(seq.upper().encode()).digest() != expected[tid]:
                raise ValueError('Longest CDS sequence mismatch: ' + tid)
            seen.add(tid)
            out.write('>' + tid + '\n' + seq + '\n')
    if seen != set(queries):raise ValueError('Longest CDS FASTA is missing selected transcripts')
    if queries:pysam.faidx(str(output/'longest.cds.fa'))
    with (output/'query_metadata.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['transcript', 'gene', 'CDS_length', 'variety'])
        for tid, g in queries.items():w.writerow([tid, g, models[tid]['cds_len'], regions[g]['variety']])
    return queries


def prepare_targets(output, regions):
    varieties = list(dict.fromkeys(r['variety'] for r in regions.values()))
    targets = []
    handles = {}
    try:
        for index, variety in enumerate(varieties):
            directory = output/'targets'/str(index + 1).zfill(4)
            directory.mkdir(parents=True)
            handles[variety] = (directory/'gene.fa').open('w')
            targets.append(dict(variety=variety, directory=str(directory.relative_to(output)), genes=0, bases=0))
        stats = {t['variety']: t for t in targets}
        seen = set()
        with (output/'body.fa').open() as f:
            for title, seq in SimpleFastaParser(f):
                g = title.split()[0]
                if g not in regions or g in seen:raise ValueError('Unexpected/duplicate gene-body sequence: ' + g)
                seen.add(g)
                if len(seq) != int(regions[g]['body_end0']) - int(regions[g]['body_start0']):
                    raise ValueError('Gene-body length mismatch: ' + g)
                v = regions[g]['variety']
                handles[v].write('>' + g + '\n' + seq + '\n')
                stats[v]['genes'] += 1
                stats[v]['bases'] += len(seq)
        if seen != set(regions):raise ValueError('Gene-body FASTA must cover every gene')
    finally:
        for handle in handles.values():handle.close()
    return targets


def filter_target(directory, variety, queries, regions, models, gene_tids, graph):
    from .cds_construct import filter_cds_graph
    primary = set()
    mapped = set()
    secondary = collections.Counter()
    by_variety = {v: collections.Counter(queries=0, edges=0) for v in dict.fromkeys(r['variety'] for r in regions.values())}
    for g in queries.values():by_variety[regions[g]['variety']]['queries'] += 1
    expected_targets = {g for g in regions if regions[g]['variety'] == variety}
    with pysam.AlignmentFile(directory/'raw.bam') as bam:
        lengths = {g: int(regions[g]['body_end0']) - int(regions[g]['body_start0']) for g in expected_targets}
        if dict(zip(bam.references, bam.lengths)) != lengths:raise ValueError('Target BAM header differs from gene bodies')

    def observe(read):
        q = read.query_name
        if q not in queries:raise ValueError('Unexpected query in target BAM: ' + q)
        if not read.is_supplementary and not read.is_secondary:
            if q in primary:raise ValueError('Duplicate primary query in target BAM: ' + q)
            primary.add(q)
        if not read.is_unmapped and not read.is_supplementary:
            mapped.add(q)
            secondary[q] += int(read.is_secondary)

    counts = filter_cds_graph(directory/'raw.bam', directory/'filtered.bam', regions, models, gene_tids, graph, record_observer=observe)
    if primary != set(queries):raise ValueError('Target BAM primary queries do not cover the longest CDS set')
    saturated = set()
    for q in mapped:by_variety[regions[queries[q]]['variety']]['mapped_queries'] += 1
    with (directory/'saturated_CDS.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['query', 'gene', 'query_variety', 'secondary'])
        for q, n in secondary.items():
            c = by_variety[regions[queries[q]]['variety']]
            for cutoff in (500, 1000, 2000):
                if n >= cutoff:
                    counts['secondary_ge_' + str(cutoff)] += 1
                    c['secondary_ge_' + str(cutoff)] += 1
            if n >= 2000:
                g = queries[q]
                saturated.add(g)
                w.writerow([q, g, regions[g]['variety'], n])
    edges = 0
    sources = set(queries.values())
    with (directory/'edges.tsv').open('w') as f:
        for target in sorted(expected_targets):
            for source in sorted(graph.predecessors(target)):
                if source not in sources:continue
                f.write(source + '\t' + target + '\n')
                edges += 1
                by_variety[regions[source]['variety']]['edges'] += 1
    counts.update(queries=len(primary), mapped_queries=len(mapped), unmapped_queries=len(primary)-len(mapped),
                  edges=edges, saturated_CDS=len(saturated), max_secondary=max(secondary.values(), default=0))
    result = dict(status='completed', variety=variety, secondary_limit=2000,
                  counts=dict(counts), query_variety_counts={v: dict(c) for v, c in by_variety.items()})
    save(directory/'summary.json', result)
    return result, saturated


def graph_stage(prepared, output, regions, models, gene_tids, threads, timer, subset=False,
                historical_edges=None, historical_genes=None, historical_varieties=None):
    from .cds_construct import cluster_cds_graph
    with timer.stage('prepare_longest_cds'):
        queries = prepare_queries(prepared, output, regions, models, gene_tids, subset,
                                  allow_empty=historical_edges is not None)
        targets = prepare_targets(output, regions)
    graph = nx.DiGraph()
    graph.add_nodes_from(regions)
    historical_varieties = set(historical_varieties or ())
    new_queries = {t:g for t,g in queries.items() if g not in (historical_genes or ())}
    historical_count = 0
    if historical_edges is not None:
        with timer.stage('historical_edges'), Path(historical_edges).open() as f:
            for line in f:
                a,b = line.split()
                if a not in historical_genes or b not in historical_genes:
                    raise ValueError('Historical edge has an unknown endpoint')
                graph.add_edge(a,b)
                historical_count += 1
        if historical_count != graph.number_of_edges():raise ValueError('Duplicate historical graph edges')
    with (output/'new.longest.cds.fa').open('w') as out, (output/'longest.cds.fa').open() as f:
        for title,seq in SimpleFastaParser(f):
            tid = title.split()[0]
            if tid in new_queries:out.write('>'+tid+'\n'+seq+'\n')
    mapping = dict(status='running', mapping=mapping_policy(longest=True), targets=targets,
                   query_count=len(queries), target_gene_count=len(regions),
                   historical_edges=historical_count, historical_varieties=sorted(historical_varieties),
                   target_execution='sequential; each task uses the shared threads budget')
    save(output/'graph.mapping.json', mapping)
    total = collections.Counter()
    saturated = set()
    matrix = []
    try:
        for target in targets:
            d = output/target['directory']
            old_target = target['variety'] in historical_varieties
            task_queries = new_queries if old_target else queries
            query_path = output/('new.longest.cds.fa' if old_target else 'longest.cds.fa')
            target.update(status='mapping', raw_bam=str((d/'raw.bam').relative_to(output)),
                          filtered_bam=str((d/'filtered.bam').relative_to(output)),
                          edges=str((d/'edges.tsv').relative_to(output)))
            target.update(query_scope='new' if old_target else 'all', query_count=len(task_queries))
            save(output/'graph.mapping.json', mapping)
            if task_queries:
                with timer.stage('mapping_' + d.name):
                    command = minimap2_map(query_path, d/'gene.fa', threads, d/'raw.bam', options=CONSTRUCT_OPTIONS)
                    save(d/'command.json', command)
                with timer.stage('filter_' + d.name):
                    result, sat = filter_target(d, target['variety'], task_queries, regions, models, gene_tids, graph)
            else:
                sat=set()
                (d/'edges.tsv').write_text('')
                (d/'saturated_CDS.tsv').write_text('query\tgene\tquery_variety\tsecondary\n')
                result=dict(status='completed',variety=target['variety'],secondary_limit=2000,
                            skipped_empty_cds=True, counts=dict(queries=0,edges=0),
                            query_variety_counts={v:dict(queries=0,edges=0) for v in dict.fromkeys(r['variety'] for r in regions.values())})
                save(d/'summary.json',result)
                target['skipped_empty_cds']=True
                target.pop('raw_bam');target.pop('filtered_bam')
            saturated.update(sat)
            total.update({k: v for k, v in result['counts'].items() if k != 'max_secondary'})
            target.update(status='completed', counts=result['counts'])
            matrix.extend(dict(query_variety=v, target_variety=target['variety'], **c)
                          for v, c in result['query_variety_counts'].items()
                          if not(old_target and v in historical_varieties))
            save(output/'graph.mapping.json', mapping)
        fields = list(dict.fromkeys(k for row in matrix for k in row))
        with (output/'query_target_variety_counts.tsv').open('w') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', restval=0)
            w.writeheader();w.writerows(matrix)
        (output/'saturated_genes_any_target.txt').write_text(''.join(g+'\n' for g in sorted(saturated)))
        mapping.update(status='completed', ordered_variety_pairs=len(matrix), totals=dict(total),
                       reused_variety_pairs=len(historical_varieties)**2,
                       total_variety_pairs=len(targets)**2,
                       audit_scope='Current generation mappings only; historical saturation and alignment statistics remain in parent provenance',
                       saturated_genes_any_target=len(saturated),
                       limitation='N2000 is a finite reporting cap; below-cap counts do not establish exhaustive alignment.')
        save(output/'graph.mapping.json', mapping)
        with timer.stage('clustering'):
            pre, last = cluster_cds_graph(output, regions, graph, total)
        for rows in (pre, last):
            members = [g for row in rows for g in row]
            if len(members) != len(set(members)) or set(members) != set(regions):
                raise ValueError('Clusters do not partition all prepared genes')
        main = {'Chr'+str(i)+s for i in range(1, 8) for s in 'ABD'}
        for row in last:
            if any(g != row[0] and not (graph.has_edge(g, row[0]) and graph.has_edge(row[0], g)) for g in row):
                raise ValueError('Last cluster member lacks reciprocal hub support')
            c = collections.Counter(regions[g]['variety'] for g in row if regions[g]['chrom'] in main)
            if any(n > 1 for n in c.values()):raise ValueError('Duplicate main-chromosome variety occupancy')
        save(output/'graph.validation.json', dict(status='completed', partition_verified=True,
             reciprocal_hub_verified=True, variety_occupancy_verified=True))
        return pre, last
    except BaseException as exc:
        mapping.update(status='failed', error=repr(exc))
        save(output/'graph.mapping.json', mapping)
        raise
