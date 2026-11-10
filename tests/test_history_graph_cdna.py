import tempfile
from pathlib import Path
from unittest.mock import patch
from pantrans import pipeline
from pantrans.main import read_parameters


def test_cli_separates_annotation_and_graph_queries():
    args = read_parameters(['legacy-append', '--cdna', 'final.fa', '--history-gff3', 'final.gff3',
                            '--history-graph-cdna', 'original.fa', '--query-cdna', 'new.fa',
                            '--gdna', 'new.gdna.fa', '--bed', 'all.bed',
                            '--history-graph', 'old.json', '--output', 'out'])
    assert args.cdna == 'final.fa'
    assert args.history_graph_cdna == 'original.fa'


def test_graph_alignment_does_not_overwrite_annotation_bam():
    with tempfile.TemporaryDirectory() as folder:
        annotation = Path(folder) / 'test_history_to_query.filtered.bam'
        annotation.write_bytes(b'annotation sentinel')
        with patch.object(pipeline, 'minimap2_map') as align, patch.object(pipeline, 'filter_bam') as filt:
            result = pipeline._align_history_graph_cdna('original.fa', 'new.fa', 32, folder, 'test', {})
        assert result != str(annotation)
        assert annotation.read_bytes() == b'annotation sentinel'
        assert align.call_args.args[:3] == ('original.fa', 'new.fa', 32)
        assert filt.call_args.args[:2] == (align.call_args.args[3], result)
        assert filt.call_args.kwargs == {'collect_edges': False}
