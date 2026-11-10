import csv,json,hashlib,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'server'))
from app import Backend
class PortableTests(unittest.TestCase):
    def test_bundled_data(self):
        rows=Backend(ROOT).list_reference_data()['references']
        self.assertEqual(len(rows),7)
        self.assertEqual(sum(x['annual_points'] for x in rows),494)
    def test_frozen_example(self):
        (ROOT/'docs').mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT/'docs') as temp:
            b=Backend(ROOT);b.output_root=Path(temp)
            r=b.run_wiggle_match(str(ROOT/'server/examples/example_samples.csv'),'demo_ce400',350,431,0,50,'portable')
            self.assertEqual(r['summary']['map_start_astro_year'],391)
            self.assertEqual(r['scientific_acceptance'],'not_automatically_evaluated')
if __name__=='__main__':unittest.main()
