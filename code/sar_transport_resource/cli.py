"""Standard-library task access and saved-prediction evaluation."""
import argparse,csv,json
from pathlib import Path
from dataclasses import asdict
from .core import TaskView,score_predictions,aggregate_scores
from .catalog import list_tasks,load_task_entry

def load_task(path):
    return TaskView.from_dict(json.loads(Path(path).read_text()))

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    ls=s.add_parser('list');ls.add_argument('--manifest');ls.add_argument('--data-root');ls.add_argument('--split');ls.add_argument('--endpoint')
    ins=s.add_parser('inspect');ins.add_argument('--data-root',required=True);ins.add_argument('--endpoint',required=True);ins.add_argument('--task-id',required=True)
    one=s.add_parser('predict-source');one.add_argument('--task',required=True);one.add_argument('--out',required=True)
    sc=s.add_parser('score');sc.add_argument('--task',required=True);sc.add_argument('--predictions',required=True);sc.add_argument('--truth',required=True)
    batch=s.add_parser('score-saved');batch.add_argument('--data-root',required=True);batch.add_argument('--manifest',required=True);batch.add_argument('--predictions-tsv',required=True);batch.add_argument('--split',default='evaluation');batch.add_argument('--out',required=True)
    a=p.parse_args()
    if a.command=='inspect':
        e=load_task_entry(a.data_root,a.endpoint,a.task_id)
        print(json.dumps({'exposure_status':e.exposure_status,'split':e.split,'task':e.view.to_dict(),
            'initial_reference_contract':e.initial_reference_contract.to_dict() if e.initial_reference_contract else None},indent=2));return
    if a.command=='list':
        if bool(a.manifest)==bool(a.data_root):p.error('Use --data-root for the current catalogue, or an explicit historical --manifest')
        if a.data_root:items=list_tasks(a.data_root)
        else:
            with open(a.manifest,newline='') as f:items=list(csv.DictReader(f,delimiter='\t'))
        items=[r for r in items if (not a.split or r['split']==a.split) and (not a.endpoint or r['endpoint']==a.endpoint)]
        print(json.dumps(items,indent=2));return
    if a.command=='score-saved':
        bytask={}
        with open(a.predictions_tsv,newline='') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                d=bytask.setdefault(r['task_id'],{});cid=r['canonical_compound_id']
                if cid in d:raise ValueError('duplicate prediction')
                d[cid]=float(r['prediction'])
        with open(a.manifest,newline='') as f:manifest=[r for r in csv.DictReader(f,delimiter='\t') if r['split']==a.split and r['track_id']=='fixed_c8_v1']
        if set(bytask)!={r['legacy_task_id'] for r in manifest}:raise ValueError('Prediction task coverage differs from requested split')
        records=[]
        for r in manifest:
            v=load_task(Path(a.data_root)/r['task_file']);t=json.loads((Path(a.data_root)/r['truth_file']).read_text())
            if t['task_id']!=v.task_id:raise ValueError('Truth task mismatch')
            records.append(score_predictions(v,bytask[r['legacy_task_id']],t['query_truth']))
        Path(a.out).write_text(json.dumps({'aggregate':aggregate_scores(records),'tasks':[asdict(x) for x in records]},indent=2)+'\n');return
    v=load_task(a.task)
    if a.command=='predict-source':
        xs={c.compound_id:c.source_potency for c in v.compounds}
        Path(a.out).write_text(json.dumps({'task_id':v.task_id,'query_signature':v.query_signature,'predictions':{x:xs[x] for x in v.query_ids}},indent=2)+'\n');return
    pr=json.loads(Path(a.predictions).read_text());tr=json.loads(Path(a.truth).read_text())
    if pr['task_id']!=v.task_id or tr['task_id']!=v.task_id:raise ValueError('Task identity mismatch')
    if pr.get('query_signature',v.query_signature)!=v.query_signature:raise ValueError('Query signature mismatch')
    print(json.dumps(asdict(score_predictions(v,pr['predictions'],tr['query_truth'])),indent=2))
if __name__=='__main__':main()
