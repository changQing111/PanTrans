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


def prepared_identity(prepared, hashes=None):
    data = manifest(str(Path(prepared).resolve()))
    if data.get('schema_version') != 2:
        raise ValueError('CDS append requires schema-2 preprocessing')
    return dict(files={key: (hashes[key] if hashes is not None else digest(prepared_file(prepared, key))) for key in KEYS},
                varieties=data['varieties'], cds_validation=data['cds_validation'])


def mapping_policy():
    return dict(options=MINIMAP2_OPTIONS,
                minimap2_version=subprocess.check_output(['minimap2','--version'],text=True).strip())


def write_package(output, prepared, kind, parents=None, prepared_hashes=None):
    data = dict(schema_version=1, kind='pantrans-cds-graph', operation=kind,
                policy=POLICY, mapping=mapping_policy(), prepared=prepared_identity(prepared, prepared_hashes),
                prepared_path=str(Path(prepared).resolve()),
                edges=dict(path='graph.edges.tsv',sha256=digest(output/'graph.edges.tsv')),
                parents=parents or [])
    (output/'cds.graph.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


def load_package(history, prepared):
    history = Path(history)
    path = history/'cds.graph.json'
    if not path.is_file():
        raise ValueError('Missing cds.graph.json; legacy CDS results need explicit --import-legacy-history')
    data=json.loads(path.read_text())
    if data.get('schema_version') != 1 or data.get('kind') != 'pantrans-cds-graph':
        raise ValueError('Not a supported CDS graph package')
    if data['policy'] != POLICY or data['mapping'] != mapping_policy():
        raise ValueError('Historical CDS filtering or mapping policy differs')
    if data['prepared'] != prepared_identity(prepared):
        raise ValueError('Historical prepared inputs changed or do not match the graph')
    edges=history/data['edges']['path']
    if digest(edges) != data['edges']['sha256']:
        raise ValueError('Historical graph edges changed')
    return data,edges
