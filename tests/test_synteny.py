"""Installed-command-friendly synteny tests with real graph/coordinate fixtures."""
import csv
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from pantrans import synteny
from pantrans.main import read_parameters


def reference_chains(points, distance, minimum):
    """Independent full predecessor scan for the bounded production DP."""
    remaining=set(points);result=[]
    while len(remaining)>=minimum:
        ordered=sorted(remaining);best=[];best_direction=1
        for direction in (1,-1):
            paths=[]
            for i,(x,y,_) in enumerate(ordered):
                path=[(x,y,1)]
                for j in range(i-1,-1,-1):
                    a,b,_=ordered[j]
                    if 0<x-a<=distance and 0<direction*(y-b)<=distance and len(paths[j])+1>len(path):
                        path=paths[j]+[(x,y,1)]
                paths.append(path)
            longest=max(paths,key=len)
            if len(longest)>len(best):best=longest;best_direction=direction
        if len(best)<minimum:break
        result.append((best_direction,best));remaining.difference_update(best)
    return result


def reference_support(x,y,key,blocks,spatial,lift):
    """Inspect every block, independent of the production spatial index."""
    from bisect import bisect_left
    choices=[]
    for bid,direction,points in blocks[key]:
        if not points[0][0]-lift<=x<=points[-1][0]+lift:continue
        pos=bisect_left(points,(x,-1))
        near=[abs(x-a)+abs(y-b) for a,b in points[max(0,pos-1):pos+1]
              if abs(x-a)<=lift and abs(y-b)<=lift]
        if not near:continue
        if 0<pos<len(points):
            if not min(points[pos-1][1],points[pos][1])<=y<=max(points[pos-1][1],points[pos][1]):continue
        elif pos==0:
            if direction*(y-points[0][1])>0:continue
        elif direction*(y-points[-1][1])<0:continue
        choices.append((min(near),bid))
    choices.sort()
    if not choices:return None,False
    if len(choices)>1 and choices[0][0]==choices[1][0]:return None,True
    return choices[0],False


class SyntenyTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.prepared=self.root/'prepared';self.prepared.mkdir()
        self.results=self.root/'results';self.results.mkdir()
        varieties=[dict(variety='A'),dict(variety='B')]
        (self.prepared/'Pan_preprocess.json').write_text(json.dumps(dict(schema_version=2,varieties=varieties,files=dict(bed='genes.bed'))))
        rows=[];clusters=[];edges=[]
        for i in range(5):
            a,b=f'a{i}',f'b{i}'
            rows += [f'Chr1A\t{i*100}\t{i*100+50}\t{a}\tA\t+',
                     f'Chr1A\t{i*100}\t{i*100+50}\t{b}\tB\t+']
            clusters.append(a+'\t'+b);edges += [a+'\t'+b,b+'\t'+a]
        # Reverse chromosome block, different chromosome labels permitted.
        for i in range(5):
            a,b=f'r{i}',f's{i}'
            rows += [f'Chr2A\t{i*100}\t{i*100+50}\t{a}\tA\t+',
                     f'Chr3D\t{(4-i)*100}\t{(4-i)*100+50}\t{b}\tB\t-']
            clusters.append(a+'\t'+b);edges += [a+'\t'+b,b+'\t'+a]
        rows += ['ctg\t0\t50\tu\tA\t+','ctg\t0\t50\tv\tB\t+']
        clusters.append('u\tv');edges += ['u\tv','v\tu']
        # A one-way same-pre relationship must not become a synteny candidate.
        rows += ['Chr1A\t1000\t1050\tx\tA\t+','Chr1A\t1000\t1050\ty\tB\t+']
        clusters.append('x\ty');edges.append('x\ty')
        self.bed=self.prepared/'genes.bed';self.bed.write_text('\n'.join(rows)+'\n')
        (self.results/'pre.cluster').write_text('\n'.join(clusters)+'\n')
        (self.results/'graph.edges.tsv').write_text('\n'.join(edges)+'\n')
        self.out=self.root/'synteny'

    def args(self,*extra):
        return read_parameters(['synteny','--prepared',str(self.prepared),'--results',str(self.results),
                                '-o',str(self.out),*map(str,extra)])

    def rows(self,path):
        with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

    def test_complete_run_forward_reverse_and_unplaced_candidates(self):
        before=self.bed.read_bytes()
        synteny.run(self.args())
        summary=json.loads((self.out/'summary.json').read_text())
        self.assertEqual(summary['comparisons'],1)
        pair=self.out/'pairs/pair_0001'
        details=json.loads((pair/'summary.json').read_text())
        self.assertEqual(details['candidate_pairs'],11)
        self.assertEqual(details['unplaced_contig_candidates'],1)
        blocks=self.rows(pair/'blocks.tsv')
        self.assertEqual({x['orientation'] for x in blocks},{'+','-'})
        self.assertEqual(len(blocks),2)
        self.assertEqual({int(x['anchors']) for x in blocks},{5})
        self.assertEqual(len(self.rows(pair/'anchors.tsv')),10)
        self.assertEqual(self.bed.read_bytes(),before)
        synteny.run(read_parameters(['synteny','--stage','export','-o',str(self.out),'--pair-index','0']))
        self.assertEqual(len(self.rows(pair/'candidates.export.tsv')),11)
        self.assertNotIn('x\ty',(pair/'candidates.export.tsv').read_text())
        with self.assertRaisesRegex(ValueError,'already exists'):synteny.run(self.args())

    def test_stages_and_cli_no_external_algorithm_dependency(self):
        base=[sys.executable,'-m','pantrans','synteny','-o',str(self.out)]
        env=dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
        for options in [ ['--stage','index','--prepared',str(self.prepared),'--results',str(self.results)],
                         ['--stage','shard'], ['--stage','pair','--pair-index','0'], ['--stage','collect'] ]:
            done=subprocess.run(base+options,env=env,capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stderr)
        self.assertTrue((self.out/'summary.json').is_file())

    def test_empty_graph_and_invalid_graph(self):
        (self.results/'graph.edges.tsv').write_text('')
        synteny.run(self.args())
        self.assertEqual(json.loads((self.out/'summary.json').read_text())['total_blocks'],0)
        (self.results/'graph.edges.tsv').write_text('unknown\ta0\n')
        self.out=self.root/'bad_edges'
        with self.assertRaisesRegex(ValueError,'Unknown gene'):synteny.run(self.args())
        self.assertFalse((self.out/'summary.json').exists())

    def test_input_identity_and_partition_validation(self):
        synteny.run(self.args('--stage','index'))
        self.bed.write_text(self.bed.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed after indexing'):
            synteny.run(read_parameters(['synteny','-o',str(self.out),'--stage','shard']))
        self.out=self.root/'bad_partition'
        self.bed.write_text(self.bed.read_text().rstrip()+'\n')
        (self.results/'pre.cluster').write_text('a0\tb0\na0\n')
        with self.assertRaisesRegex(ValueError,'duplicate membership'):synteny.run(self.args())

    def test_invalid_parameters_and_removed_command(self):
        for extra in [ ['--max-copy','0'],['--distance','0'],['--min-anchors','1'],['--lift-distance','-1'],
                       ['--edges','x'],['--stage','pair'],['--destination','x'] ]:
            with self.subTest(extra=extra),self.assertRaises(SystemExit):self.args(*extra)
        with self.assertRaises(SystemExit):read_parameters(['translocation'])

    def test_chaining_matches_exhaustive_reference(self):
        rng=random.Random(768)
        for _ in range(25):
            points=list({(rng.randrange(40),rng.randrange(40),1) for i in range(50)})
            self.assertEqual(synteny.monotonic_chains(points,5,3),reference_chains(points,5,3))
        # Component connectivity must include repeated-axis neighbors just as the pilot kernel does.
        self.assertEqual(synteny.components([(0,0,1),(0,2,1),(2,3,1),(20,20,1)],3,2),
                         [[(0,0,1),(0,2,1),(2,3,1)]])

    def test_position_support_and_ambiguity_match_exhaustive_reference(self):
        from collections import defaultdict
        rng=random.Random(716)
        for lift in [0,1,3,10]:
            key=('Chr1A','Chr2D');blocks=[];spatial=defaultdict(set)
            for bi in range(15):
                direction=rng.choice([-1,1]);x=rng.randrange(30);y=rng.randrange(40,80);pts=[]
                for _ in range(6):
                    x+=rng.randrange(1,10);y+=direction*rng.randrange(1,10);pts.append((x,y))
                blocks.append((str(bi),direction,pts))
                for x,y in pts:spatial[key,x//(lift+1)].add(bi)
            for _ in range(500):
                x,y=rng.randrange(100),rng.randrange(100)
                self.assertEqual(synteny.find_position_support(x,y,key,{key:blocks},spatial,lift),
                                 reference_support(x,y,key,{key:blocks},spatial,lift))

    def test_three_varieties_and_custom_chromosome_list(self):
        mf=self.prepared/'Pan_preprocess.json';data=json.loads(mf.read_text())
        data['varieties'].append(dict(variety='C'));mf.write_text(json.dumps(data))
        with self.bed.open('a') as f:f.write('Chr1A\t0\t50\tc0\tC\t+\n')
        pre=self.results/'pre.cluster';text=pre.read_text().replace('a0\tb0','a0\tb0\tc0');pre.write_text(text)
        edges=self.results/'graph.edges.tsv'
        with edges.open('a') as f:f.write('a0\tc0\nc0\ta0\nb0\tc0\nc0\tb0\n')
        main=self.root/'main.txt';main.write_text('Chr1A\n')
        from pantrans.prepared_bundle import manifest
        manifest.cache_clear()
        synteny.run(self.args('--main-chroms',main))
        summary=json.loads((self.out/'summary.json').read_text())
        self.assertEqual(summary['comparisons'],3)
        self.assertEqual(summary['total_blocks'],1)


if __name__=='__main__':unittest.main()
