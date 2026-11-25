#!/usr/bin/env python3
"""Plot PanTrans candidate blocks without importing or changing pantrans.

Requires Python >= 3.6 and matplotlib >= 3.3. See plot_translocations.md.
All genomic intervals are 0-based, half-open. No rearrangements are inferred.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys


STATUSES = ('supported_block', 'low_support_block')
TEAL = '#247f78'
ORANGE = '#bc803e'
GRAY = '#b5bdc5'


def read_table(path, required):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError('{}: missing columns: {}'.format(path, ', '.join(sorted(missing))))
        rows = []
        for number, row in enumerate(reader, 2):
            if None in row or any(row.get(k) in (None, '') for k in required):
                raise ValueError('{}:{}: incomplete or malformed row'.format(path, number))
            rows.append(row)
        return rows


def interval(row, side):
    start, end = int(row[side + '_start']), int(row[side + '_end'])
    if start < 0 or end <= start:
        raise ValueError('Invalid {} interval: {}-{}'.format(side, start, end))
    return start, end


def natural(value):
    return tuple((0, int(x)) if x.isdigit() else (1, x.lower())
                 for x in re.split(r'(\d+)', value) if x)


def token(value):
    # Append a digest if sanitizing, preventing distinct IDs from colliding.
    clean = re.sub(r'[^A-Za-z0-9_.-]', '_', value).strip('.') or 'item'
    if clean != value:
        clean += '-' + hashlib.sha256(value.encode()).hexdigest()[:8]
    return clean


def load_sizes(path):
    sizes = {}
    if path:
        for row in read_table(path, ('variety', 'chrom', 'length')):
            key = row['variety'], row['chrom']
            length = int(row['length'])
            if length <= 0 or key in sizes:
                raise ValueError('Chromosome sizes must be positive and keys unique: {}'.format(key))
            sizes[key] = length
    return sizes


def load_bed(path):
    positions = {}
    if path:
        with open(path, encoding='utf-8-sig') as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) < 5:
                    raise ValueError('{}:{}: BED needs at least five columns'.format(path, number))
                chrom, start, end, gene, variety = fields[:5]
                value = (variety, chrom, int(start), int(end))
                if not gene or not variety or not chrom or value[2] < 0 or value[3] <= value[2]:
                    raise ValueError('{}:{}: invalid BED record'.format(path, number))
                if gene in positions and positions[gene] != value:
                    raise ValueError('Conflicting BED records for ' + gene)
                positions[gene] = value
    return positions


def load_inputs(args):
    required = ('block_id', 'reference_variety', 'query_variety', 'reference_chrom',
                'query_chrom', 'reference_start', 'reference_end', 'query_start',
                'query_end', 'orientation', 'anchor_count', 'pre_cluster_count', 'status')
    blocks = read_table(args.blocks, required)
    ids = set()
    for block in blocks:
        if block['block_id'] in ids:
            raise ValueError('Duplicate block ID: ' + block['block_id'])
        ids.add(block['block_id'])
        if block['orientation'] not in ('+', '-') or block['status'] not in STATUSES:
            raise ValueError('Invalid block orientation/status: ' + block['block_id'])
        interval(block, 'reference')
        interval(block, 'query')
        if not 1 <= int(block['pre_cluster_count']) <= int(block['anchor_count']):
            raise ValueError('Invalid anchor counts: ' + block['block_id'])
    if len({b['reference_variety'] for b in blocks}) > 1:
        raise ValueError('Use a report with a single reference variety')
    if args.query and set(args.query) - {b['query_variety'] for b in blocks}:
        raise ValueError('Unknown query variety in --query')
    if args.block_id and set(args.block_id) - ids:
        raise ValueError('Unknown block ID in --block-id')
    selected = [b for b in blocks if b['status'] in args.status
                and (not args.query or b['query_variety'] in args.query)
                and (not args.block_id or b['block_id'] in args.block_id)]
    sizes, positions = load_sizes(args.chrom_sizes), load_bed(args.bed)
    for b in selected:
        for side in ('reference', 'query'):
            key = b[side + '_variety'], b[side + '_chrom']
            if sizes and (key not in sizes or interval(b, side)[1] > sizes[key]):
                raise ValueError('Missing chromosome size or block exceeds chromosome: {}'.format(key))
    genes = []
    if args.genes:
        genes = read_table(args.genes, ('block_id', 'reference_gene', 'query_gene',
            'reference_variety', 'query_variety', 'reference_chrom', 'query_chrom',
            'reference_start', 'reference_end', 'query_start', 'query_end', 'status'))
        for row in genes:
            for side in ('reference', 'query'):
                start, end = interval(row, side)
                gene = row[side + '_gene']
                pos = row[side + '_variety'], row[side + '_chrom'], start, end
                if gene in positions and positions[gene] != pos:
                    raise ValueError('Gene report and BED disagree: ' + gene)
                positions[gene] = pos
    anchors = {}
    for b in selected:
        members = [r for r in genes if r['block_id'] == b['block_id']]
        if args.view == 'local' and len(members) != int(b['anchor_count']):
            raise ValueError('Anchor count differs between reports: ' + b['block_id'])
        if len({(r['reference_gene'], r['query_gene']) for r in members}) != len(members):
            raise ValueError('Duplicate anchor pair: ' + b['block_id'])
        for r in members:
            for side in ('reference', 'query'):
                if any(r[side + '_' + k] != b[side + '_' + k] for k in ('variety', 'chrom')):
                    raise ValueError('Anchor identity differs from block: ' + b['block_id'])
                start, end = interval(r, side)
                bs, be = interval(b, side)
                if start < bs or end > be:
                    raise ValueError('Anchor falls outside block: ' + b['block_id'])
        anchors[b['block_id']] = members
    edges = set()
    if args.edges:
        directed = set()
        with open(args.edges, encoding='utf-8-sig') as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 2 or not all(fields):
                    raise ValueError('{}:{}: expected two gene IDs'.format(args.edges, number))
                a, b = fields
                if a not in positions or b not in positions:
                    raise ValueError('Graph endpoint absent from BED: {} / {}'.format(a, b))
                if a != b:
                    directed.add((a, b))
        edges = {(a, b) for a, b in directed if (b, a) in directed}
    return blocks, selected, sizes, positions, genes, anchors, edges


def ordered_chroms(chroms, path):
    order = []
    if path:
        with open(path, encoding='utf-8-sig') as handle:
            order = [s.strip() for s in handle if s.strip() and not s.lstrip().startswith('#')]
        if len(order) != len(set(order)):
            raise ValueError('Duplicate entries in --chrom-order')
    return [c for c in order if c in chroms] + sorted(set(chroms) - set(order), key=natural)


def connection(ax, patches, top, bottom, **style):
    """Cubic Bezier with vertical tangents where it meets the gene tracks."""
    from matplotlib.path import Path as MplPath
    path = MplPath([(top, 1), (top, .62), (bottom, .38), (bottom, 0)],
                   [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    ax.add_patch(patches.PathPatch(path, facecolor='none', **style))


def ribbon(ax, patches, rs, re_, qs, qe, reverse=False, **style):
    """Bezier boundaries preserve interval ends, including reverse-order blocks."""
    from matplotlib.path import Path as MplPath
    q_right, q_left = (qs, qe) if reverse else (qe, qs)
    vertices = [(rs, 1), (re_, 1), (re_, .62), (q_right, .38), (q_right, 0),
                (q_left, 0), (q_left, .38), (rs, .62), (rs, 1), (rs, 1)]
    codes = [MplPath.MOVETO, MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CURVE4, MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(MplPath(vertices, codes), **style))


def overview(plt, patches, blocks, sizes, args):
    ref, query = blocks[0]['reference_variety'], blocks[0]['query_variety']
    chroms = ordered_chroms({c for v, c in sizes if v in (ref, query)}, args.chrom_order)
    longest = max(sizes[v, c] for v, c in sizes if v in (ref, query))
    gap = longest * .08
    offsets, cursor = {}, 0
    for c in chroms:
        offsets[c] = cursor
        cursor += max(sizes.get((ref, c), 0), sizes.get((query, c), 0)) + gap
    total = cursor - gap
    fig, ax = plt.subplots(figsize=(args.width or max(10, len(chroms)*.55), args.height or 5))
    fig.subplots_adjust(left=.07, right=.98, bottom=.23, top=.78)
    for b in blocks:
        rs, re_ = interval(b, 'reference')
        qs, qe = interval(b, 'query')
        rs, re_ = rs + offsets[b['reference_chrom']], re_ + offsets[b['reference_chrom']]
        qs, qe = qs + offsets[b['query_chrom']], qe + offsets[b['query_chrom']]
        color = TEAL if b['orientation'] == '+' else ORANGE
        supported = b['status'] == 'supported_block'
        ribbon(ax, patches, rs, re_, qs, qe, reverse=b['orientation']=='-',
               facecolor=color, edgecolor=color, alpha=.32 if supported else .16,
               linewidth=.6, linestyle='-' if supported else '--')
        # A center line keeps very short spans visible at chromosome scale.
        connection(ax, patches, (rs+re_)/2, (qs+qe)/2, edgecolor=color,
                   lw=.7, alpha=.8, ls='-' if supported else '--', zorder=2)
        for start, end, y in [(rs, re_, 1), (qs, qe, 0)]:
            ax.plot([start, end], [y, y], color=color, lw=5, zorder=4)
        if args.label_blocks:
            ax.text((rs+re_)/2, .88, b['block_id'], fontsize=6, rotation=55, rotation_mode='anchor', ha='right')
    for variety, y in [(ref, 1), (query, 0)]:
        for chrom in chroms:
            if (variety, chrom) not in sizes:
                continue
            start, end = offsets[chrom], offsets[chrom]+sizes[variety, chrom]
            ax.plot([start, end], [y, y], color=GRAY, lw=7, solid_capstyle='butt', zorder=3)
            ax.text((start+end)/2, y + (.15 if y else -.15), chrom,
                    ha='center', va='center', fontsize=7, rotation=45 if len(chroms)>12 else 0, rotation_mode='anchor')
    ax.set_xlim(-total*.015, total*1.015)
    ax.set_ylim(-.32, 1.32)
    ax.axis('off')
    fig.suptitle('Candidate interchromosomal blocks: {} vs {}'.format(query, ref), fontsize=12, y=.96)
    fig.text(.07, .87, 'Top: {}     Bottom: {}     Blocks: {}'.format(ref, query, len(blocks)), fontsize=9)
    fig.text(.07, .12, 'Teal: + order   |   Ochre: - order   |   Solid: supported   |   Dashed: low support', fontsize=8)
    scale = 10 ** math.floor(math.log10(longest / 3))
    ax.plot([0, scale], [-.3, -.3], color='#263442', lw=1.5)
    fig.text(.07, .065, 'Scale bar: {:g} Mb. Shared bp scale; chromosome slots include display gaps.'.format(scale/1e6), fontsize=8)
    return fig, {'reference': ref, 'query': query, 'block_ids': [b['block_id'] for b in blocks]}


def local(plt, patches, block, members, sizes, positions, genes, edges, args):
    windows = {}
    for side in ('reference', 'query'):
        start, end = interval(block, side)
        key = block[side+'_variety'], block[side+'_chrom']
        right = min(end+args.flank, sizes[key]) if sizes else end+args.flank
        windows[side] = (max(0, start-args.flank), right)

    def xcoord(side, value):
        start, end = windows[side]
        return (value-start)/(end-start)

    def visible(gene, side):
        variety, chrom, start, end = positions[gene]
        lo, hi = windows[side]
        return (variety == block[side+'_variety'] and chrom == block[side+'_chrom']
                and start < hi and end > lo)

    fig, ax = plt.subplots(figsize=(args.width or 11, args.height or 5))
    fig.subplots_adjust(left=.12, right=.96, bottom=.24, top=.73)
    anchor_pairs = {(r['reference_gene'], r['query_gene']) for r in members}
    background = {(a, b) for a, b in edges if visible(a, 'reference') and visible(b, 'query')}
    if args.show_ambiguous:
        background.update((r['reference_gene'], r['query_gene']) for r in genes
                          if r['status'] == 'ambiguous_copy'
                          and visible(r['reference_gene'], 'reference') and visible(r['query_gene'], 'query'))
    background -= anchor_pairs
    for a, b in sorted(background):
        connection(ax, patches, xcoord('reference', sum(positions[a][2:])/2),
                   xcoord('query', sum(positions[b][2:])/2),
                   edgecolor='#9ca7b3', lw=.7, ls=':', alpha=.65, zorder=1)
    for side, y in [('reference', 1), ('query', 0)]:
        start, end = interval(block, side)
        xs, xe = xcoord(side, start), xcoord(side, end)
        ax.add_patch(patches.Rectangle((xs, y-.08), xe-xs, .16, facecolor='#e4f1ee', edgecolor='none', zorder=1))
        ax.plot([0, 1], [y, y], color=GRAY, lw=1, zorder=2)
        for boundary in [xs, xe]:
            ax.plot([boundary, boundary], [y-.12, y+.12], color='#608d86', ls='--', lw=.8)
        anchor_ids = {r[side+'_gene'] for r in members}
        for gene, pos in positions.items():
            if not visible(gene, side):
                continue
            a, b = max(0, xcoord(side, pos[2])), min(1, xcoord(side, pos[3]))
            color = TEAL if gene in anchor_ids else GRAY
            ax.plot([a, b], [y, y], color=color, lw=6, solid_capstyle='butt', zorder=4)
            if gene in anchor_ids:
                ax.plot([(a+b)/2], [y], marker='|', color=color, markersize=7, zorder=4)
        label = '{} / {}'.format(block[side+'_variety'], block[side+'_chrom'])
        ax.text(0, y+(.25 if y else -.25), label, ha='left', fontsize=9)
        lo, hi = windows[side]
        for fraction in (0, .5, 1):
            value = lo + fraction*(hi-lo)
            ax.text(fraction, y+(.13 if y else -.16), '{:,.0f}'.format(value),
                    ha={0:'left', .5:'center', 1:'right'}[fraction], fontsize=7, color='#596675')
    color = TEAL if block['orientation'] == '+' else ORANGE
    for row in members:
        xs = [xcoord(side, sum(interval(row, side))/2) for side in ('reference', 'query')]
        connection(ax, patches, xs[0], xs[1], edgecolor=color, lw=1.1,
                   ls='-' if block['status']=='supported_block' else '--', alpha=.8, zorder=3)
        if args.label_genes:
            for side, x, y in [('reference', xs[0], .90), ('query', xs[1], .10)]:
                ax.text(x, y, row[side+'_gene'], fontsize=6, rotation=60, rotation_mode='anchor',
                        ha='right' if side=='reference' else 'left', va='center', zorder=5)
    ax.set_xlim(-.015, 1.015)
    ax.set_ylim(-.32, 1.32)
    ax.axis('off')
    fig.suptitle('Candidate block {} | {} vs {}'.format(block['block_id'], block['query_variety'], block['reference_variety']), fontsize=12, y=.96)
    fig.text(.12, .87, '{} | order {} | {} anchors / {} pre-clusters | flank {:,} bp'.format(
        block['status'], block['orientation'], block['anchor_count'], block['pre_cluster_count'], args.flank), fontsize=9)
    fig.text(.12, .15, 'Each track has its own increasing bp scale. Shading marks anchor span, not breakpoints.', fontsize=8)
    fig.text(.12, .105, 'Gray: background genes. Dotted: additional matches within these two windows only.', fontsize=8)
    fig.text(.12, .06, ('BED context provided.' if args.bed else 'No BED: only genes present in the candidate report are shown.') +
             (' Upper flank clipped to chromosome length.' if sizes else 'No chromosome sizes: upper flank may extend beyond assembly.'), fontsize=7)
    return fig, {'block_id': block['block_id'], 'windows': windows,
                 'anchors': len(members), 'additional_matches': len(background)}


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--blocks', required=True, type=Path)
    p.add_argument('--genes', type=Path, help='Gene report; required for local view')
    p.add_argument('--chrom-sizes', type=Path, help='Headered TSV: variety, chrom, length; required for overview')
    p.add_argument('--bed', type=Path, help='Complete BED: chrom, start, end, gene, variety')
    p.add_argument('--edges', type=Path, help='Optional directed edge TSV; requires BED; local view only')
    p.add_argument('-o', '--output-prefix', required=True, type=Path)
    p.add_argument('--view', choices=('overview', 'local'), default='overview')
    p.add_argument('--query', nargs='+')
    p.add_argument('--block-id', nargs='+')
    p.add_argument('--status', nargs='+', choices=STATUSES, default=list(STATUSES))
    p.add_argument('--chrom-order', type=Path, help='One chromosome per line; unlisted names follow naturally')
    p.add_argument('--flank', type=int, default=100000, help='Local extension in bp on each side (default: 100000)')
    p.add_argument('--show-ambiguous', action='store_true', help='Overlay ambiguous report pairs within local windows')
    p.add_argument('--label-blocks', action='store_true', help='Label overview links (may crowd dense plots)')
    p.add_argument('--label-genes', action='store_true', help='Label local anchor genes (may crowd dense plots)')
    p.add_argument('--formats', nargs='+', choices=('pdf', 'png', 'svg'), default=['pdf', 'png'])
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--width', type=float, help='Figure width in inches')
    p.add_argument('--height', type=float, help='Figure height in inches')
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.flank < 0 or args.dpi <= 0 or any(x is not None and (not math.isfinite(x) or x <= 0) for x in (args.width, args.height)):
        p.error('Require flank >= 0 and positive finite dimensions / dpi')
    if args.view == 'overview' and not args.chrom_sizes:
        p.error('--chrom-sizes is required for overview')
    if args.view == 'local' and not args.genes:
        p.error('--genes is required for local view')
    if args.edges and (not args.bed or args.view != 'local'):
        p.error('--edges requires --bed and --view local')
    if args.show_ambiguous and args.view != 'local':
        p.error('--show-ambiguous requires --view local')
    try:
        all_blocks, blocks, sizes, positions, genes, anchors, edges = load_inputs(args)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import patches
        plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans'],
                             'font.size':8, 'pdf.fonttype':42, 'svg.fonttype':'none'})
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        outputs = []
        inputs = {str(value.resolve()) for key, value in vars(args).items()
                  if isinstance(value, Path) and key != 'output_prefix'}
        def save(fig, suffix, details):
            paths = [Path(str(args.output_prefix)+'.'+suffix+'.'+fmt) for fmt in dict.fromkeys(args.formats)]
            if any(str(path.resolve()) in inputs for path in paths):
                raise ValueError('Output would overwrite an input')
            for path in paths:
                if path.suffix == '.pdf':
                    fig.savefig(path, format='pdf', facecolor='white')
                elif path.suffix == '.svg':
                    fig.savefig(path, format='svg', facecolor='white')
                else:
                    fig.savefig(path, format='png', dpi=args.dpi, facecolor='white')
                print(path)
            plt.close(fig)
            outputs.append(dict(details, files=[str(path) for path in paths]))
        if args.view == 'overview':
            for query in sorted({b['query_variety'] for b in blocks}, key=natural):
                subset = [b for b in blocks if b['query_variety'] == query]
                fig, details = overview(plt, patches, subset, sizes, args)
                save(fig, token(query)+'.overview', details)
        else:
            for b in blocks:
                fig, details = local(plt, patches, b, anchors[b['block_id']], sizes, positions, genes, edges, args)
                save(fig, token(b['query_variety'])+'.'+token(b['block_id'])+'.local', details)
        manifest = Path(str(args.output_prefix)+'.plot_manifest.json')
        if str(manifest.resolve()) in inputs:
            raise ValueError('Manifest would overwrite an input')
        manifest.write_text(json.dumps({'parameters': {k:str(v) if isinstance(v, Path) else v for k,v in vars(args).items()},
            'input_blocks':len(all_blocks), 'selected_blocks':len(blocks), 'outputs':outputs,
            'interpretation':'Candidate anchor spans, not validated rearrangements or breakpoints.'}, indent=2)+'\n', encoding='utf-8')
        print('Selected {} / {} blocks; manifest: {}'.format(len(blocks), len(all_blocks), manifest))
        return 0
    except (ValueError, OSError, ImportError) as exc:
        p.exit(1, 'Error: {}\n'.format(exc))


if __name__ == '__main__':
    sys.exit(main())
