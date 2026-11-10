"""Run with a wheel-installed Python, outside the checkout and without PYTHONPATH."""
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

import pantrans
import pysam


def main():
    expected_prefix=Path(sys.prefix).resolve()
    package_path=Path(pantrans.__file__).resolve()
    assert expected_prefix in package_path.parents,(expected_prefix,package_path)
    executable=Path(sys.executable).parent/'pantrans'
    assert executable.is_file()
    env=dict(os.environ)
    env.pop('PYTHONPATH',None)
    with tempfile.TemporaryDirectory(prefix='pantrans-wheel-smoke-') as directory:
        root=Path(directory)
        def run(*args,ok=True):
            result=subprocess.run([str(executable),*map(str,args)],cwd=root,env=env,capture_output=True,text=True)
            assert (result.returncode==0)==ok,result.stdout+'\n'+result.stderr
            return result
        assert '0.8.0' in run('--version').stdout
        assert '--prepared' in run('construct','--help').stdout
        assert '--query-prepared' in run('append','--help').stdout
        assert '--cdna' in run('legacy-construct','--help').stdout
        assert '--history-gff3' in run('legacy-append','--help').stdout
        run('construct','-c','old.fa','-o','bad',ok=False)
        rng=random.Random(772)
        codons=['AAA','AAC','ACC','AGC','ATC','CAA','CAC','CCA','CGC','CTC','GAA','GAC','GCA','GGC','GTC','TAC','TCC','TGC','TTC']
        cds='ATG'+''.join(rng.choices(codons,k=298))+'TAA'
        genomes=[];annotations=[]
        for name in ('X','CS'):
            seq='C'*40+cds+'G'*40;gene=name+'_gene'
            genome=root/(name+'.fa');genome.write_text('>Chr1A\n'+seq+'\n')
            annotation=root/(name+'.gff3')
            annotation.write_text('##gff-version 3\n'+'\n'.join([
                f'Chr1A\tx\tgene\t1\t{len(seq)}\t.\t+\t.\tID={gene}',
                f'Chr1A\tx\tmRNA\t1\t{len(seq)}\t.\t+\t.\tID={gene}.1;Parent={gene}',
                f'Chr1A\tx\texon\t1\t{len(seq)}\t.\t+\t.\tParent={gene}.1',
                f'Chr1A\tx\tCDS\t41\t{40+len(cds)}\t.\t+\t0\tParent={gene}.1'])+'\n')
            genomes.append(genome);annotations.append(annotation)
            run('preprocess','--genomes',genome,'--annotations',annotation,'--varieties',name,'-o',root/(name+'_prepared'))
        run('construct','--prepared',root/'X_prepared','-o',root/'history','-t',2,'--workers',2)
        run('append','--history-results',root/'history','--history-prepared',root/'X_prepared',
            '--query-prepared',root/'CS_prepared','-o',root/'appended','-t',2,'--workers',2)
        run('preprocess','--genomes',*genomes,'--annotations',*annotations,'--varieties','X','CS','-o',root/'all_prepared')
        run('construct','--prepared',root/'all_prepared','-o',root/'direct','-t',2,'--workers',2)
        for name in ('graph.edges.tsv','pre.cluster','last.cluster','last.gff3','last.cdna.fa','last.model_sources.tsv','last.alignment_status.tsv','last.summary.json'):
            assert (root/'appended'/name).read_bytes()==(root/'direct'/name).read_bytes(),name
        with pysam.AlignmentFile(root/'appended/last.annotation.bam') as a,pysam.AlignmentFile(root/'direct/last.annotation.bam') as b:
            assert [r.to_string() for r in a]==[r.to_string() for r in b]
        assert not (root/'direct/pre.gff3').exists()
        assert (root/'direct/last.evidence.sqlite').exists()
        policy=json.loads((root/'direct/cds.graph.json').read_text())['policy']
        assert policy['filter']=='cds-homology-v2' and policy['pid']==0.9 and not policy['structure_checks']
        summary={'wheel_installed':True,'package_path':str(package_path),'console':str(executable),
                 'version':pantrans.__version__,'real_construct_append_equivalent':True,
                 'last_only':True,'graph_policy':policy,'legacy_help_available':True}
        print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
