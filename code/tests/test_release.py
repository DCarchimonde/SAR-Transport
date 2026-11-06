import json,tempfile,unittest
from pathlib import Path
from sar_transport_resource.core import *

class ReleaseTests(unittest.TestCase):
    def test_initial_four_are_paid_and_invalid_access_never_calls_provider(self):
        ids=['p'+str(i) for i in range(10)]+['q1','q2']
        v=TaskView(ACQUISITION_TRACK+'::task',ACQUISITION_TRACK,'task','g','s','t','Ki','pKi','target',
            tuple(Compound(x,'C',1.) for x in ids),tuple(ids[:10]),(),('q1','q2'),degree_weighted_edges([('q1','q2')]))
        calls=[]
        def provider(cid):calls.append(cid);return 1.
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'trace.jsonl';o=BudgetOracle(v,provider,8,path)
            for x in ids[:4]:o.reveal(x)
            self.assertEqual(o.spent,4);self.assertEqual(o.remaining,4)
            for x in ids[:4]:o.reveal(x)
            self.assertEqual(len(calls),4)
            for x in ids[4:8]:o.reveal(x)
            for x,exc in [('p8',RuntimeError),('q1',ValueError),('unknown',KeyError)]:
                with self.assertRaises(exc):o.reveal(x)
            self.assertEqual(len(calls),8);self.assertEqual(o.spent,8)
            self.assertEqual(o.current_view().query_signature,v.query_signature)
            trace=[json.loads(x) for x in path.read_text().splitlines()]
            self.assertEqual(sum(x.get('cost',0) for x in trace),8)
            self.assertEqual(sum(x['event']=='rejected' for x in trace),3)

    def test_extra_truth_fields_rejected(self):
        p=Path(__file__).resolve().parents[1]/'examples/fixed_c8_task.json'
        raw=json.loads(p.read_text());raw['query_truth']={}
        with self.assertRaises(ValueError):TaskView.from_dict(raw)

if __name__=='__main__':unittest.main()
