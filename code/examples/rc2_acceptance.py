"""Exactly one saved fixed task and one saved acquisition task; no matrix run."""
import argparse,importlib.util,json,platform,sys,time
from pathlib import Path
from dataclasses import asdict
from sar_transport_resource import load_task_entry,score_predictions
from sar_transport_resource.rebuild import rebuild_exported_task

def main():
    p=argparse.ArgumentParser();p.add_argument('--data-root',required=True);p.add_argument('--out-dir',required=True);a=p.parse_args()
    start=time.process_time();root=Path(a.data_root).resolve();out=Path(a.out_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    code=Path(__file__).resolve().parents[1]
    v=load_task_entry(root,'IC50','fixed_c8_v1::CHEMBL1000620__CHEMBL1000621').view
    saved=json.loads((code/'examples/saved_predictions.json').read_text())
    truth=json.loads((root/'IC50/fixed_c8_v1/truth/CHEMBL1000620__CHEMBL1000621.json').read_text())
    expected=json.loads((code/'private_fixtures/expected_scores.json').read_text())['expected_from_original_PER_TASK_METRICS']
    fixed={}
    for name in ('source','tanimoto'):
        score=asdict(score_predictions(v,saved['methods'][name],truth['query_truth']))
        diffs={k:abs(score[k]-expected[name][k]) for k in ('L','MAE','delta_MSE')}
        assert max(diffs.values())<1e-12
        fixed[name]={'score':score,'absolute_difference':diffs}
    rebuilt=rebuild_exported_task(root,'IC50',v.task_id)
    assert rebuilt.to_dict()==v.to_dict()
    spec=importlib.util.spec_from_file_location('formal_demo',code/'examples/formal_saved_policy.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    acq=module.run(root,out/'formal_example_trace.jsonl')
    av=load_task_entry(root,'IC50',acq['run']['task_id']).view
    rebuilt_a=rebuild_exported_task(root,'IC50',av.task_id)
    assert rebuilt_a.to_dict()==av.to_dict()
    result={'status':'PASS','scope':'software-only saved-prediction regression, one real task per track',
       'real_fixed_tasks':1,'real_acquisition_tasks':1,'fixed':fixed,'acquisition':acq,
       'portable_table_reconstruction':{'two_task_semantic_identity':True},
       'python':sys.version,'platform':platform.platform(),'installed_module':__import__('sar_transport_resource').__file__,
       'Windows_tested':False,'macOS_tested':False,'cpu_seconds':time.process_time()-start,
       'new_scientific_fits':0,'old_model_inference':0,'scientific_matrix_scores':0,'resampling':0}
    (out/'REAL_TASK_REGRESSION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','real_fixed_tasks','real_acquisition_tasks','cpu_seconds')}))
if __name__=='__main__':main()
