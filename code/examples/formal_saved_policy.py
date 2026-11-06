"""One saved formal IC50 acquisition trajectory: integration, not a new policy.

The callback only chooses the four already-recorded additional references and
returns frozen saved predictions. No predictor is fitted or inferred. The
evaluator opens the query truth after the callback has finished.
"""
import argparse,json,csv
from pathlib import Path
from dataclasses import asdict
from sar_transport_resource import FormalAcquisitionRunner,load_task_entry,score_predictions

DEFAULT_TASK='CHEMBL1013099__CHEMBL1013112'

def run(data_root,trace_path,legacy_task_id=DEFAULT_TASK):
    root=Path(data_root)
    entry=load_task_entry(root,'IC50','acquisition_sentinel_v1::'+legacy_task_id)
    view=entry.view
    saved=next(r for r in json.loads((root/'IC50/P0_development_evidence/PREDICTIONS.json').read_text())
               if r['task_id']==legacy_task_id and r['policy']=='maxmin')
    assert saved['context8_ids'][:4]==list(entry.initial_reference_contract.initial_reference_ids)
    def callback(paid_view,oracle):
        assert oracle.spent==4 and oracle.remaining==4
        assert [x for x,_ in paid_view.reference_labels]==saved['context8_ids'][:4]
        for cid in saved['context8_ids'][4:]:oracle.reveal(cid)
        return dict(zip(saved['query_ids'],saved['predicted']))
    pool=json.loads((root/'IC50/acquisition_sentinel_v1/pool_truth'/f'{legacy_task_id}.json').read_text())
    assert pool['task_id']==view.task_id
    runner=FormalAcquisitionRunner(view,entry.initial_reference_contract,pool['pool_truth'].__getitem__,trace_path)
    result=runner.run(callback)
    truth=json.loads((root/'IC50/acquisition_sentinel_v1/truth'/f'{legacy_task_id}.json').read_text())
    assert truth['task_id']==view.task_id
    score=asdict(score_predictions(view,result['predictions'],truth['query_truth']))
    with (root/'IC50/P0_development_evidence/PER_TASK_METRICS.tsv').open() as f:
        expected=next(r for r in csv.DictReader(f,delimiter='\t') if r['task_id']==legacy_task_id and r['policy']=='maxmin')
    errors={k:abs(score[k]-float(expected[k])) for k in ('L','MAE')}
    assert max(errors.values())<1e-12 and result['protocol_compliant']
    return {'scope':'one saved-trajectory software regression; not a new scientific comparison',
      'exposure_status':entry.exposure_status,'run':result,'score':score,'original_score_absolute_difference':errors,
      'scientific_fits':0,'old_model_inference':0,'matrix_scores':0}

def main():
    p=argparse.ArgumentParser();p.add_argument('--data-root',required=True);p.add_argument('--trace',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();Path(a.out).write_text(json.dumps(run(a.data_root,a.trace),indent=2)+'\n')
if __name__=='__main__':main()
