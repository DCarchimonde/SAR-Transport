"""Export saved C4 identities and coverage, without prediction or scoring."""
from pathlib import Path
import collections, csv, hashlib, json, sys, time
ROOT=Path(__file__).resolve().parents[1]
PUB=ROOT/'PUBLIC_RELEASE_CANDIDATE_RC2';DATA=PUB/'resource_data';OUT=ROOT/'reports/sar_resource_rc2'
sys.path.insert(0,str(PUB/'code'))
from sar_transport_resource.acquisition import object_signature,PROTOCOL_ID

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2)+'\n')
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def main():
    start=time.process_time();manifest=rows(DATA/'CURRENT_TASK_MANIFEST.tsv')
    icp=DATA/'IC50/P0_development_evidence/TASK_DESIGNS.json'
    ik={d['task_id']:d for d in json.loads(icp.read_text())}
    kp=DATA/'Ki/reference_results/KI_ACQUISITION_TRACES.json'
    kt=collections.defaultdict(list)
    for r in json.loads(kp.read_text()):kt[r['task_id']].append(r)
    old_ic_tr=json.loads((DATA/'IC50/P0_development_evidence/ACQUISITION_TRACES.json').read_text())
    old_ic_pred=json.loads((DATA/'IC50/P0_development_evidence/PREDICTIONS.json').read_text())
    checks=[]
    for m in manifest:
        if m['track_id']!='acquisition_sentinel_v1':continue
        p=DATA/m['task_file'];raw=json.loads(p.read_text());cid=m['legacy_task_id']
        if m['endpoint']=='IC50':
            initial=ik[cid]['initial4_ids'];source=icp
            matching=[x for x in old_ic_pred if x['task_id']==cid]
            trajectories=[x for x in old_ic_tr if x['task_id']==cid]
            assert len(matching)==len(trajectories)==9
            for pred in matching:
                assert pred['context8_ids'][:4]==initial
                trajectory=next(x for x in trajectories if x['policy']==pred['policy'])
                assert initial+[x['selected_id'] for x in trajectory['trace']]==pred['context8_ids']
            count=len(trajectories)
        else:
            matching=kt[m['task_id']];assert len(matching)==9
            initial=matching[0]['initial4'];source=kp
            assert all(x['initial4']==initial and x['query_signature']==m['query_signature'] for x in matching)
            assert all(len(x['steps'])==4 and len(set(initial+[s['selected_id'] for s in x['steps']]))==8 for x in matching)
            count=len(matching)
        assert len(initial)==len(set(initial))==4 and set(initial)<=set(raw['reference_pool_ids'])
        contract={'protocol_id':PROTOCOL_ID,'task_id':m['task_id'],'endpoint':m['endpoint'],
          'initial_reference_ids':initial,'total_budget':8,'additional_unique_references':4,
          'query_signature':m['query_signature'],'task_definition_sha256':sha(p),
          'source_file':str(source.relative_to(DATA)),'source_sha256':sha(source)}
        contract['protocol_signature']=object_signature(contract)
        dump(DATA/m['initial_reference_contract'],contract)
        checks.append({'endpoint':m['endpoint'],'task_id':m['task_id'],'initial_reference_ids':initial,
                       'formal_trajectories_checked':count,'identical_initial_ids':True,
                       'protocol_signature':contract['protocol_signature'],'query_signature':m['query_signature']})
    dump(OUT/'INITIAL_REFERENCE_CONTRACT.json',{'protocol_id':PROTOCOL_ID,'task_contracts':len(checks),
       'tasks_by_endpoint':dict(collections.Counter(x['endpoint'] for x in checks)),
       'saved_trajectories_checked':sum(x['formal_trajectories_checked'] for x in checks),
       'ids_reselected':False,'scientific_predictions_or_scores_changed':False,
       'formal_runner':'sar_transport_resource.FormalAcquisitionRunner',
       'schema':'code/schemas/acquisition_contract.schema.json','tasks':checks})
    schema={'$schema':'https://json-schema.org/draft/2020-12/schema',
        '$id':'urn:sar-transport:acquisition-contract:rc2','title':'Saved common-C4 acquisition sidecar',
        'type':'object','additionalProperties':False,'required':list(contract),
        'properties':{k:{'type':'string'} for k in contract}}
    schema['properties'].update({'initial_reference_ids':{'type':'array','items':{'type':'string'},'minItems':4,'maxItems':4,'uniqueItems':True},
        'total_budget':{'const':8},'additional_unique_references':{'const':4},'protocol_id':{'const':PROTOCOL_ID}})
    dump(PUB/'code/schemas/acquisition_contract.schema.json',schema)
    coverage=[]
    for (endpoint,track,split),rr in sorted(_groups(manifest,lambda x:(x['endpoint'],x['track_id'],x['split'])).items()):
        coverage.append({'asset':'task_definitions','endpoint':endpoint,'track':track,'split':split,'method_or_policy':'all',
            'tables':1,'tasks':len(rr),'groups':len({r['group_id'] for r in rr}),'targets':len({r['target_id'] for r in rr}),
            'query_rows':sum(int(r['query_points']) for r in rr),'query_edges':sum(int(r['query_edges']) for r in rr),
            'status':'exposed; task definitions do not imply saved predictions'})
    evalm={r['legacy_task_id']:r for r in manifest if r['endpoint']=='IC50' and r['track_id']=='fixed_c8_v1' and r['split']=='evaluation'}
    for p in sorted((DATA/'IC50/saved_predictions').glob('*.tsv')):
        rr=rows(p);tt={r['task_id'] for r in rr}
        if not rr or 'prediction' not in rr[0]:
            assert p.name=='TASK_ADAPTATIONS.tsv' and len(rr)==252 and tt==set(evalm)
            coverage.append({'asset':'saved_support_adaptation_log','endpoint':'IC50','track':'fixed_c8_v1','split':'evaluation',
              'method_or_policy':'ActFound seed17/43','tables':1,'tasks':126,'groups':92,'targets':75,
              'query_rows':'NA (252 task/seed adaptation records)','query_edges':25568,
              'status':'adaptation metadata, not a query-prediction table'})
            continue
        assert tt==set(evalm) and len(rr)==2898
        coverage.append({'asset':'saved_prediction_table','endpoint':'IC50','track':'fixed_c8_v1','split':'evaluation',
            'method_or_policy':p.stem,'tables':1,'tasks':len(tt),'groups':92,'targets':75,'query_rows':len(rr),
            'query_edges':25568,'status':'saved reference predictions; no development/validation predictions supplied'})
    for endpoint,trs in [('IC50',old_ic_tr),('Ki',json.loads(kp.read_text()))]:
        for policy,rr in sorted(_groups(trs,lambda x:x['policy']).items()):
            mm=[r for r in manifest if r['endpoint']==endpoint and r['track_id']=='acquisition_sentinel_v1']
            coverage.append({'asset':'saved_acquisition_trajectory','endpoint':endpoint,'track':'acquisition_sentinel_v1',
              'split':mm[0]['split'],'method_or_policy':policy,'tables':1,'tasks':len(rr),'groups':len({x['group_id'] for x in mm}),
              'targets':len({x['target_id'] for x in mm}),'query_rows':sum(int(x['query_points']) for x in mm),
              'query_edges':sum(int(x['query_edges']) for x in mm),'status':'original saved common4+4; not rerun'})
    # Ki reference-model prediction coverage is recorded without reading/scoring new values.
    for method in ('source','affine','rbf','tanimoto'):
        coverage.append({'asset':'saved_prediction_method','endpoint':'Ki','track':'fixed_c8_v1','split':'frozen_Ki_application',
            'method_or_policy':method,'tables':'combined_JSON','tasks':64,'groups':50,'targets':39,
            'query_rows':1111,'query_edges':7873,'status':'original saved reference predictions'})
    for p in (OUT/'COVERAGE.tsv',PUB/'COVERAGE.tsv'):
        with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(coverage[0]),delimiter='\t');w.writeheader();w.writerows(coverage)
    dump(OUT/'CONTRACT_EXPORT_BUDGET.json',{'cpu_seconds':time.process_time()-start,'scientific_fits':0,
        'scientific_scores':0,'resampling':0,'old_model_inference':0})
    print(json.dumps({'contracts':len(checks),'saved_trajectories_checked':sum(x['formal_trajectories_checked'] for x in checks),
                      'prediction_tables_IC50':sum(x['asset']=='saved_prediction_table' for x in coverage)}))

def _groups(items,key):
    out=collections.defaultdict(list)
    for x in items:out[key(x)].append(x)
    return out
if __name__=='__main__':main()
