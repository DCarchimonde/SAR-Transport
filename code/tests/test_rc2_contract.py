import copy,hashlib,json,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from sar_transport_resource import (ACQUISITION_TRACK,Compound,TaskView,BudgetOracle,
    degree_weighted_edges,AcquisitionContract,FormalAcquisitionRunner,audit_acquisition_trace)
from sar_transport_resource.acquisition import PROTOCOL_ID,object_signature

def fixture():
    pool=tuple('p'+str(i) for i in range(10));q=('q1','q2','q3')
    view=TaskView(ACQUISITION_TRACK+'::fixture',ACQUISITION_TRACK,'fixture','group','a','b','Ki','pKi','target',
        tuple(Compound(x,'C',1.) for x in pool+q),pool,(),q,degree_weighted_edges([('q1','q2'),('q2','q3')]))
    c={'protocol_id':PROTOCOL_ID,'task_id':view.task_id,'endpoint':'Ki',
       'initial_reference_ids':['p3','p0','p5','p2'],'total_budget':8,'additional_unique_references':4,
       'query_signature':view.query_signature,'task_definition_sha256':'0'*64,
       'source_file':'synthetic_fixture','source_sha256':'0'*64}
    c['protocol_signature']=object_signature(c)
    return view,AcquisitionContract.from_dict(c)

class RC2Tests(unittest.TestCase):
    def test_common_four_paid_before_callback_and_four_remaining(self):
        v,c=fixture();calls=[]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'trace';r=FormalAcquisitionRunner(v,c,lambda x:(calls.append(x) or 1.),p)
            self.assertEqual(calls,list(c.initial_reference_ids))
            def callback(view,o):
                self.assertEqual((o.spent,o.remaining),(4,4))
                self.assertEqual([x for x,_ in view.reference_labels],list(c.initial_reference_ids))
                for cid in c.initial_reference_ids:o.reveal(cid)
                self.assertEqual(len(calls),4)
                for cid in ('p1','p4','p6','p7'):o.reveal(cid)
                self.assertEqual(o.spent,8)
                for cid,exc in [('p8',RuntimeError),('q1',ValueError),('unknown',KeyError)]:
                    with self.assertRaises(exc):o.reveal(cid)
                self.assertEqual(len(calls),8)
                self.assertEqual(o.current_view().query_edges,v.query_edges)
                return {x:1. for x in v.query_ids}
            out=r.run(callback);self.assertTrue(out['protocol_compliant'])
            self.assertEqual(out['query_signature'],v.query_signature)
            self.assertEqual(out['protocol_signature'],c.protocol_signature)
            self.assertEqual(out['additional_unique_references'],4)
            with self.assertRaises(RuntimeError):r.run(callback)

    def test_generic_eight_cannot_be_called_formal(self):
        v,c=fixture()
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'generic';o=BudgetOracle(v,lambda _:1.,8,p)
            for cid in v.reference_pool_ids[:8]:o.reveal(cid)
            report=audit_acquisition_trace(v,c,[json.loads(x) for x in p.read_text().splitlines()])
            self.assertFalse(report['protocol_compliant'])
            self.assertIn('missing_or_mismatched_formal_header',report['failure_reasons'])

    def test_under_spending_is_not_a_formal_comparison(self):
        v,c=fixture()
        with tempfile.TemporaryDirectory() as td:
            out=FormalAcquisitionRunner(v,c,lambda _:1.,Path(td)/'t').run(lambda view,o:{x:1. for x in view.query_ids})
            self.assertFalse(out['protocol_compliant']);self.assertEqual(out['cost_used'],4)

    def test_changed_initials_fail_signature_before_provider_access(self):
        v,c=fixture();raw=c.to_dict();raw['initial_reference_ids'][0]='p9'
        with self.assertRaises(ValueError):AcquisitionContract.from_dict(raw)
        with self.assertRaises(ValueError):replace(c,total_budget=9)
        with self.assertRaises(ValueError):replace(c,additional_unique_references=5)

    def test_wrong_query_or_endpoint_is_rejected_before_access(self):
        v,c=fixture();calls=[]
        with tempfile.TemporaryDirectory() as td:
            wrong=replace(v,endpoint='IC50',potency_unit='pIC50')
            with self.assertRaises(ValueError):FormalAcquisitionRunner(wrong,c,lambda x:calls.append(x),Path(td)/'t')
            wrong=replace(v,query_edges=degree_weighted_edges([('q1','q2'),('q1','q3')]))
            with self.assertRaises(ValueError):FormalAcquisitionRunner(wrong,c,lambda x:calls.append(x),Path(td)/'t')
            self.assertEqual(calls,[])

    def test_bad_prediction_keys_and_nonfinite_values_rejected(self):
        v,c=fixture()
        with tempfile.TemporaryDirectory() as td:
            for n,pred in enumerate([{'q1':1.},{x:float('nan') for x in v.query_ids}]):
                r=FormalAcquisitionRunner(v,c,lambda _:1.,Path(td)/str(n))
                with self.assertRaises(ValueError):r.run(lambda view,o:pred)

    def test_tampered_billing_is_detected(self):
        v,c=fixture()
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'t';r=FormalAcquisitionRunner(v,c,lambda _:1.,p)
            def callback(view,o):
                for x in ('p1','p4','p6','p7'):o.reveal(x)
                return {x:1. for x in view.query_ids}
            r.run(callback)
            e=[json.loads(x) for x in p.read_text().splitlines()];e[2]['cost']=0
            self.assertFalse(audit_acquisition_trace(v,c,e)['protocol_compliant'])

if __name__=='__main__':unittest.main()
