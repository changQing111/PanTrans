import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('comparison',Path(__file__).resolve().parents[1]/'scripts/compare_cds_construct_append.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class CompareCDSTests(unittest.TestCase):
    def test_multiset_and_differences(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            result=m.compare(iter(['b','a','a']),iter(['a','b','c']),root,'x')
            self.assertEqual((result['shared'],result['append_only'],result['construct_only']),(2,1,1))
            self.assertEqual((root/'x.append_only.txt').read_text(),'a\n')
            self.assertEqual((root/'x.construct_only.txt').read_text(),'c\n')

    def test_model_identity_independent_of_transcript_number(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);a=root/'a';b=root/'b'
            text='g\tx\tgene\t1\t10\t.\t+\t.\tID=g\ng\tx\tmRNA\t1\t10\t.\t+\t.\tID=g.1;Parent=g;completeness=complete\ng\tx\texon\t1\t10\t.\t+\t.\tParent=g.1\n'
            a.write_text(text);b.write_text(text.replace('g.1','g.7'))
            self.assertEqual(list(m.models(a)),list(m.models(b)))
            b.write_text(text.replace('exon\t1\t10','exon\t2\t10'))
            self.assertNotEqual(list(m.models(a)),list(m.models(b)))
