"""Versioned CDS graph provenance; intentionally separate from legacy graph packages."""
import hashlib
import json
import subprocess
from pathlib import Path

from .align_filter import MINIMAP2_OPTIONS
from .prepared_bundle import manifest, prepared_file

KEYS = ('gene', 'cds_valid', 'cdna', 'gdna', 'genes', 'models')
POLICY = dict(filter='cds-homology-v2', pid=0.90, query_cds_coverage=0.90,
              target_cds_coverage=0.90, structure_checks=False, reference='CS',
              main_chroms=['Chr'+str(i)+s for i in range(1,8) for s in 'ABD'])


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


CONSTRUCT_OPTIONS = "-ax splice:hq -uf -I 2G -K 50M --secondary=yes -f 20000 -N 2000"


def prepared_identity(prepared, hashes=None, longest=False):
    data = manifest(str(Path(prepared).resolve()))
    if data.get('schema_version') != 2:
        raise ValueError('CDS append requires schema-2 preprocessing')
    keys = KEYS + (('cds_longest', 'cds_longest_selection') if longest else ())
    return dict(files={key: (hashes[key] if hashes is not None else digest(prepared_file(prepared, key))) for key in keys},
                varieties=data['varieties'], cds_validation=data['cds_validation'])


def mapping_policy(longest=False):
    data = dict(options=CONSTRUCT_OPTIONS if longest else MINIMAP2_OPTIONS,
                minimap2_version=subprocess.check_output(['minimap2','--version'],text=True).strip())
    if longest:
        data.update(query='longest-valid-CDS-per-gene', targets='per-variety-complete-gene-bodies',
                    within_variety=True, tie_break='output-transcript-ID-ascending')
    return data


def write_package(output, prepared, kind, parents=None, prepared_hashes=None, longest=False):
    data = dict(schema_version=1, kind='pantrans-cds-graph', operation=kind,
                policy=POLICY, mapping=mapping_policy(longest), prepared=prepared_identity(prepared, prepared_hashes, longest),
                prepared_path=str(Path(prepared).resolve()),
                edges=dict(path='graph.edges.tsv',sha256=digest(output/'graph.edges.tsv')),
                parents=parents or [])
    if longest:
        data['target_manifest'] = dict(path='graph.mapping.json', sha256=digest(output/'graph.mapping.json'))
    (output/'cds.graph.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


def load_package(history, prepared, longest=False):
    history = Path(history)
    path = history/'cds.graph.json'
    if not path.is_file():
        raise ValueError('Missing cds.graph.json; legacy CDS results need explicit --import-legacy-history')
    data=json.loads(path.read_text())
    if data.get('schema_version') != 1 or data.get('kind') != 'pantrans-cds-graph':
        raise ValueError('Not a supported CDS graph package')
    if data['policy'] != POLICY or data['mapping'] != mapping_policy(longest):
        raise ValueError('Historical CDS filtering or mapping policy differs; use compatible longest-CDS per-variety N2000 history for current construct/append')
    if data['prepared'] != prepared_identity(prepared, longest=longest):
        raise ValueError('Historical prepared inputs changed or do not match the graph')
    edges=history/data['edges']['path']
    if digest(edges) != data['edges']['sha256']:
        raise ValueError('Historical graph edges changed')
    if longest and digest(history/data['target_manifest']['path']) != data['target_manifest']['sha256']:
        raise ValueError('Historical target mapping manifest changed')
    return data,edges
