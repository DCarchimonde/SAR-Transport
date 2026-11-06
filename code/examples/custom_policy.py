"""NON-FORMAL integration demo. Use formal_saved_policy.py for common-C4 replay.

The trivial policy measures eight lexicographically first pool IDs and returns
the source ordering. It illustrates integration; it is NOT an added benchmark
strategy, and is never included in the Ki results.
"""
import argparse,json
from pathlib import Path
from dataclasses import asdict
from sar_transport_resource import TaskView,BudgetOracle,score_predictions

def policy(view,oracle):
    for cid in sorted(view.reference_pool_ids)[:8]:oracle.reveal(cid)
    xs={c.compound_id:c.source_potency for c in view.compounds}
    return {cid:xs[cid] for cid in view.query_ids}

def main():
    p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--pool',required=True);p.add_argument('--truth',required=True);p.add_argument('--trace',required=True);a=p.parse_args()
    view=TaskView.from_dict(json.loads(Path(a.task).read_text()))
    pool=json.loads(Path(a.pool).read_text());assert pool['task_id']==view.task_id
    oracle=BudgetOracle(view,pool['pool_truth'].__getitem__,8,a.trace)
    predictions=policy(view,oracle)
    # The example driver/evaluator, never policy(), opens query truth after prediction.
    truth=json.loads(Path(a.truth).read_text());assert truth['task_id']==view.task_id
    print(json.dumps({'demonstration_only':True,'protocol_compliant':False,'benchmark_comparable':False,
       'reason':'custom starting references, not the formal common-C4 track','new_fits':0,'cost_used':oracle.spent,
       'score':asdict(score_predictions(view,predictions,truth['query_truth']))},indent=2))
if __name__=='__main__':main()
