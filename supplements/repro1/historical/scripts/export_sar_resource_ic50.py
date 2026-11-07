"""Export existing IC50 tasks/results only: zero fits and zero new inference."""
from pathlib import Path
import csv,json,gzip,shutil,time,sys,collections
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'resources/sar_transport_v1'))
from sar_transport_resource.core import *
from run_sar_resource_ki import rows,table,dump,make_view,write_task,RES,OUT,digest

def main():
    start=time.process_time();R=ROOT/'reports/sar_transport_v1'
    obs={(r['assay_id'],r['canonical_compound_id']):r['pchembl'] for r in rows(ROOT/'data/sar_transport/restricted/UNIQUE_OBSERVATIONS.tsv.gz')}
    identity=list(rows(R/'COMPOUND_IDENTITY.tsv'))
    smiles={r['canonical_compound_id']:r['canonical_isomeric_smiles'] for r in identity if r['status']=='eligible'}
    contexts=collections.defaultdict(list);edges=collections.defaultdict(list)
    for r in rows(R/'C8.tsv'):contexts[r['task_id']].append((int(r['context_order']),r['canonical_compound_id']))
    for r in rows(R/'QUERY_EDGES.tsv.gz'):edges[r['task_id']].append((r['compound_a'],r['compound_b']))
    tasks={r['task_id']:r for r in rows(R/'TASKS.tsv') if r['status']=='eligible'}
    manifest=[]
    for tid,t in tasks.items():
        c=[x for _,x in sorted(contexts[tid])]
        t={**t,'target_id':json.loads(t['target_ids'])[0]}
        v=make_view(t,'IC50',obs,smiles,c,c,edges[tid],FIXED_TRACK)
        manifest.append(write_task(v,obs,t['split'],'Historical IC50 data and evaluation outputs already exposed; resource export does not restore blindness'))
    P=ROOT/'data/sar_resource_v1/planning_inputs/reference_policy_p0_20261104'
    for design in json.loads((P/'TASK_DESIGNS.json').read_text()):
        t={**design,'target_id':json.loads(tasks[design['task_id']]['target_ids'])[0]}
        v=make_view(t,'IC50',obs,smiles,t['pool_ids'],[],t['edges'],ACQUISITION_TRACK)
        manifest.append(write_task(v,obs,'development_P0','Completed exploratory P0; no rerun and no confirmatory interval'))
    table(RES/'IC50/TASK_MANIFEST.tsv',manifest)
    for filename in ['TASKS.tsv','C8.tsv','QUERY_POINTS.tsv','QUERY_EDGES.tsv.gz','ALL_ELIGIBLE_EDGES.tsv.gz','COMPOUND_IDENTITY.tsv','ASSAY_SOURCE_MAPPING.tsv','OBSERVATION_ELIGIBILITY.tsv']:
        dest=RES/'IC50/provenance'/filename;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(R/filename,dest)
    shutil.copytree(ROOT/'data/sar_transport/metadata',RES/'IC50/provenance/metadata',dirs_exist_ok=True)
    dest=RES/'IC50/saved_predictions';dest.mkdir(parents=True,exist_ok=True)
    for f in (ROOT/'runs/sar_transport_v1/evaluation_predictions').glob('*.tsv'):shutil.copy2(f,dest/f.name)
    resultdest=RES/'IC50/reference_results';resultdest.mkdir(parents=True,exist_ok=True)
    for filename in ['ALL_METRICS.json','PRIMARY_CONTRASTS.json','PER_SEED_METRICS.json','VALIDATION_SELECTION.json','RULES.json','RUN_AND_BUDGET_LEDGER.json']:
        shutil.copy2(R/filename,resultdest/filename)
    shutil.copytree(P,RES/'IC50/P0_development_evidence',dirs_exist_ok=True)
    counts={}
    for tr in [FIXED_TRACK,ACQUISITION_TRACK]:
        m=[r for r in manifest if r['track_id']==tr]
        counts[tr]={'tasks':len(m),'groups':len({r['group_id'] for r in m}),'targets':len({r['target_id'] for r in m}),
           'query_edges':sum(r['query_edges'] for r in m),'query_task_points':sum(r['query_points'] for r in m),
           'splits':dict(collections.Counter(r['split'] for r in m))}
    dump(OUT/'IC50_EXPORT.json',{'counts':counts,'CPU_seconds':time.process_time()-start,'new_fits':0,'new_model_inference':0,
           'original_query_identity_preserved':True,'original_exposure_history_preserved':True})
    print(json.dumps(counts,indent=2))
if __name__=='__main__':main()
