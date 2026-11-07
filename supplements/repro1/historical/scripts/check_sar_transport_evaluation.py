"""Synthetic checks for the decision metric; no source outcomes or fitting."""
import ast
from decimal import Decimal
import datetime as dt
from pathlib import Path
import numpy as np
from sar_transport_methods import per_task_metrics, aggregate, prediction_view
from sar_transport_eligibility import OUT, dump, sha, ROOT

def main():
    t={'task_id':'synthetic','group_id':'g1','target_id':'T','source_unknown':False,
       'ids':['B','A','C'],'q':np.array([0,1,2]),'c8':np.array([0]),
       'yt':np.array([4.,6.,5.]),'xs':np.array([4.,5.,6.]),
       'edge_a':np.array([0,1]),'edge_b':np.array([1,2]),'weights':np.array([.25,.75])}
    tied=per_task_metrics(t,np.array([3.,3.,3.]))
    assert tied['L']==0., tied  # ID A is preferred on both edges.
    wrong=per_task_metrics(t,np.array([7.,4.,6.]))
    assert np.isclose(wrong['L'],1.25)
    assert np.isclose(wrong['delta_MSE'],13.)
    near_tie=per_task_metrics(t,np.array([3.+.5e-8,3.,3.]))
    assert near_tie['L']==0.
    # Two tasks in g1 must not give g1 twice the weight of the single g2 task.
    rr=[dict(wrong,group_id='g1',L=2.,MAE=2.,delta_MSE=2.),
        dict(wrong,group_id='g1',L=4.,MAE=4.,delta_MSE=4.),
        dict(wrong,group_id='g2',L=9.,MAE=9.,delta_MSE=9.)]
    assert aggregate(rr)=={'L':6.,'MAE':6.,'delta_MSE':6.}
    v=prediction_view(t)
    assert 'yt' not in v and list(v['yc'])==[4.]
    evalfile=ROOT/'scripts/evaluate_sar_transport.py'
    tree=ast.parse(evalfile.read_text())
    pred=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='predict')
    calls=[n.func.id for n in ast.walk(pred) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
    assert 'load_eval_truth' not in calls and 'per_task_metrics' not in calls
    assert Decimal('3.3')-Decimal('3.0')>=Decimal('.3')
    assert Decimal('3.299999999')-Decimal('3.0')<Decimal('.3')
    result={'status':'PASS','synthetic_tests':['canonical-ID tie','1e-8 tie tolerance',
      'weighted selection loss','weighted gap MSE','task-then-group weights','query-label-free view',
      'prediction stage does not call truth reader or scorer','exact Decimal margin boundary'],
      'source_outcomes_read':0,'fits':0,'optimizer_steps':0,
      'evaluator_sha256':sha(evalfile),'checked_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    dump(OUT/'EVALUATION_UNIT_CHECKS.json',result)
    print(result)

if __name__=='__main__':main()
