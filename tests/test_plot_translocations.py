"""Standalone plotting contracts: coordinate windows, filtering and validation."""
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'plot_translocations.py'
spec = importlib.util.spec_from_file_location('plot_translocations', SCRIPT)
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.block = dict(block_id='TR000001', reference_variety='R', query_variety='Q',
            reference_chrom='Chr1', query_chrom='Chr2', reference_start=100, reference_end=400,
            query_start=500, query_end=800, orientation='-', anchor_count=3,
            pre_cluster_count=3, status='supported_block')
        self.write('blocks.tsv', [self.block])
        rows = []
        bed = []
        for i in range(3):
            row = dict(self.block)
            row.update(reference_gene='r'+str(i), query_gene='q'+str(i),
                       reference_start=100+i*100, reference_end=200+i*100,
                       query_start=700-i*100, query_end=800-i*100)
            rows.append(row)
            bed += ['Chr1\t{}\t{}\tr{}\tR'.format(100+i*100, 200+i*100, i),
                    'Chr2\t{}\t{}\tq{}\tQ'.format(700-i*100, 800-i*100, i)]
        self.write('genes.tsv', rows)
        self.write('sizes.tsv', [dict(variety='R', chrom='Chr1', length=1000),
                                 dict(variety='Q', chrom='Chr2', length=900)])
        (self.root/'genes.bed').write_text('\n'.join(bed)+'\n', encoding='utf-8')
        self.args = ['--blocks', str(self.root/'blocks.tsv'), '--genes', str(self.root/'genes.tsv'),
            '--chrom-sizes', str(self.root/'sizes.tsv'), '--bed', str(self.root/'genes.bed'),
            '-o', str(self.root/'result'), '--formats', 'png', '--dpi', '50']

    def write(self, name, rows):
        with (self.root/name).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def run_plot(self, extra):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(plot.main(self.args+extra), 0)
        return json.loads((self.root/'result.plot_manifest.json').read_text())

    def test_flanks_clip_to_both_chromosome_ends(self):
        manifest = self.run_plot(['--view','local','--flank','1000'])
        self.assertEqual(manifest['outputs'][0]['windows'], {'reference':[0,1000], 'query':[0,900]})
        self.assertTrue((self.root/'result.Q.TR000001.local.png').is_file())

    def test_zero_flank_preserves_reverse_coordinates(self):
        manifest = self.run_plot(['--view','local','--flank','0'])
        self.assertEqual(manifest['outputs'][0]['windows'], {'reference':[100,400], 'query':[500,800]})

    def test_overview_and_empty_selection(self):
        self.assertEqual(self.run_plot([])['selected_blocks'], 1)
        manifest = self.run_plot(['--status','low_support_block'])
        self.assertEqual(manifest['outputs'], [])

    def test_inconsistent_reports_fail(self):
        self.block['anchor_count'] = 4
        self.write('blocks.tsv', [self.block])
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            plot.main(self.args+['--view','local'])

    def test_negative_flank_fails(self):
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            plot.main(self.args+['--flank','-1'])

    def test_bezier_paths_preserve_endpoints_and_reverse_order(self):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt, patches
        from matplotlib.path import Path as MplPath
        fig, ax = plt.subplots()
        self.addCleanup(plt.close, fig)
        plot.connection(ax, patches, .2, .8, edgecolor='black')
        path = ax.patches[0].get_path()
        self.assertEqual(list(path.codes), [MplPath.MOVETO] + [MplPath.CURVE4]*3)
        self.assertEqual(tuple(path.vertices[0]), (.2, 1))
        self.assertEqual(tuple(path.vertices[-1]), (.8, 0))
        plot.ribbon(ax, patches, .1, .3, .6, .9, reverse=True, facecolor='gray')
        path = ax.patches[1].get_path()
        self.assertEqual(tuple(path.vertices[4]), (.6, 0))
        self.assertEqual(tuple(path.vertices[5]), (.9, 0))
        self.assertEqual(sum(path.codes == MplPath.CURVE4), 6)

    def test_additional_edges_require_reciprocal_support(self):
        with (self.root/'genes.bed').open('a') as handle:
            handle.write('Chr1\t10\t30\ta\tR\nChr2\t810\t830\tb\tQ\n')
        edges = self.root/'edges.tsv'
        edges.write_text('a\tb\n')
        extra = ['--view','local','--edges',str(edges)]
        self.assertEqual(self.run_plot(extra)['outputs'][0]['additional_matches'], 0)
        edges.write_text('a\tb\nb\ta\na\tb\n')
        self.assertEqual(self.run_plot(extra)['outputs'][0]['additional_matches'], 1)


if __name__ == '__main__':
    unittest.main()
