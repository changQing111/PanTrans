"""Export the exact representative DNA and verify every reconstructed transcript."""
import collections
import json
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pysam

from .prepared_bundle import prepared_file
from .cds_target_graph import save


def export_last(prepared, output, clusters):
    wanted = {row[0] for row in clusters}
    sequences = {}
    dna = output/'last.gdna.fa'
    with prepared_file(prepared, 'gdna').open() as src, Path(str(dna)+'.partial').open('w') as out:
        for title, seq in SimpleFastaParser(src):
            g = title.split()[0]
            if g in wanted:
                if g in sequences:raise ValueError('Duplicate representative DNA: ' + g)
                sequences[g] = seq
                out.write('>' + g + '\n' + seq + '\n')
    if set(sequences) != wanted:raise ValueError('Missing representative DNA')
    features = collections.Counter()
    genes, transcripts, ids = set(), set(), set()
    counts = collections.Counter()
    complete_genes = set()
    model = None
    exons = []
    complete = output/'last.complete.cdna.fa'
    selected = output/'last.complete.transcripts.txt'
    with (output/'last.cdna.fa').open() as original, Path(str(complete)+'.partial').open('w') as out, Path(str(selected)+'.partial').open('w') as chosen:
        records = iter(SimpleFastaParser(original))

        def check():
            if model is None:return
            tid, parent, ref, status = model
            if not exons or exons != sorted(exons) or any(b > c for (a,b),(c,d) in zip(exons,exons[1:])):
                raise ValueError('Invalid exon structure: ' + tid)
            seq = ''.join(sequences[ref][a:b] for a,b in exons)
            record = next(records, None)
            if record is None or record[0].split()[0] != tid or record[1].upper() != seq.upper():
                raise ValueError('Reconstructed cDNA mismatch: ' + tid)
            counts['transcripts'] += 1
            counts[status] += 1
            if status == 'complete':
                complete_genes.add(parent)
                out.write('>' + tid + '\n' + seq + '\n')
                chosen.write(tid + '\n')

        with (output/'last.gff3').open() as f:
            for line in f:
                if line.startswith('#'):continue
                x = line.rstrip().split('\t')
                if len(x) != 9 or x[0] not in sequences or x[6] != '+':raise ValueError('Invalid GFF reference/strand')
                a, b = int(x[3])-1, int(x[4])
                if not 0 <= a < b <= len(sequences[x[0]]):raise ValueError('GFF feature outside representative DNA')
                attrs = dict(z.split('=',1) for z in x[8].split(';') if '=' in z)
                features[x[2]] += 1
                ident = attrs.get('ID')
                if ident:
                    if ident in ids:raise ValueError('Duplicate feature ID: ' + ident)
                    ids.add(ident)
                if x[2] in ('gene','mRNA'):
                    check();model=None;exons=[]
                if x[2] == 'gene':
                    if attrs['ID'] != x[0]:raise ValueError('Gene ID/reference mismatch')
                    genes.add(attrs['ID'])
                elif x[2] == 'mRNA':
                    if attrs['Parent'] not in genes or attrs['Parent'] != x[0]:raise ValueError('Invalid mRNA parent')
                    status = attrs['completeness']
                    if status not in ('complete','partial'):raise ValueError('Invalid completeness')
                    transcripts.add(attrs['ID'])
                    model = (attrs['ID'],attrs['Parent'],x[0],status)
                else:
                    if model is None or attrs['Parent'] != model[0] or x[0] != model[2]:raise ValueError('Invalid child parent')
                    if x[2] == 'exon':exons.append((a,b))
        check()
        if next(records,None) is not None:raise ValueError('Extra output cDNA records')
    expected = json.loads((output/'last.summary.json').read_text())
    if len(genes) != expected['genes'] or counts['transcripts'] != expected['transcripts'] or counts['complete'] != expected.get('complete',0):
        raise ValueError('Annotation counts differ from reconstructed outputs')
    missing = wanted-genes
    # Genes lacking any original transcript remain graph nodes, with no invented model.
    with (output/'last.unannotated_representatives.tsv').open('w') as f:
        f.write('gene\treason\n')
        for g in sorted(missing):f.write(g+'\tno_transcript_model\n')
    for path in (dna, complete, selected):Path(str(path)+'.partial').replace(path)
    for path in (dna,complete):
        if path.stat().st_size:pysam.faidx(str(path))
    report = dict(status='completed', representatives=len(wanted), annotated_genes=len(genes),
                  unannotated_representatives=len(missing), transcripts=counts['transcripts'],
                  complete_mRNAs=counts['complete'], partial_mRNAs=counts['partial'],
                  genes_with_complete_mRNA=len(complete_genes), features=dict(features),
                  all_cdna_reconstructed_and_verified=True, gff_bounds_verified=True,
                  id_parent_verified=True, coordinate_system='Gene-oriented representative DNA with configured flanks; local coordinates')
    save(output/'validation.json', report)
    return report
