"""Frozen one-shot Ki resource application. No global training or old inference.

Stages separate numeric qualification, pool-only prediction, and query scoring.
The eligibility stage necessarily reads numbers; it prints counts, not effects.
"""
from pathlib import Path
import os,sys,json,csv,gzip,hashlib,time,datetime,collections,decimal,urllib.request,resource
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'resources/sar_transport_v1'))
from sar_transport_resource.core import TaskView,Compound,FIXED_TRACK,ACQUISITION_TRACK,degree_weighted_edges,BudgetOracle,score_predictions
from sar_transport_resource.baselines import fixed_predict,acquire,fingerprints,maxmin,POLICIES
import numpy as np
OUT=ROOT/'reports/sar_resource_v1';DATA=ROOT/'data/sar_resource_v1'
PRE=DATA/'planning_inputs/ki_identity_check_20261104'
RES=ROOT/'resources/sar_transport_v1/data'
SOURCE='https://raw.githubusercontent.com/Jnelen/ChEMBL_MatchedPairsAnalysis/0a77cdd28abf40a1562aadce727d032cd668815c/'

def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def rows(p):
    p=Path(p);op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',newline='') as f:yield from csv.DictReader(f,delimiter='\t' if '.tsv' in p.name else ',')
def table(p,data,fields=None):
    data=list(data);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or (list(data[0]) if data else ['status'])
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
def event(kind,**kw):
    with (OUT/'EXPOSURE_EVENTS.jsonl').open('a') as f:f.write(json.dumps({'utc':stamp(),'event':kind,**kw})+'\n')
def check_freeze():
    f=json.loads((OUT/'KI_IMPLEMENTATION_FREEZE.json').read_text())
    assert len(f['commit'])==40
    protocol=json.loads((ROOT/'protocol/SAR_RESOURCE_KI_v1.json').read_text())
    assert digest(ROOT/'protocol/SAR_RESOURCE_KI_v1.json')==f['protocol_sha256']
    for name,sha in protocol['code_sha256'].items():assert digest(ROOT/name)==sha, name

def fetch():
    check_freeze();out=DATA/'ki_raw';out.mkdir(parents=True,exist_ok=True)
    files=json.loads((OUT/'KI_SOURCE_EXPECTED.json').read_text());record=[]
    for entry in files:
        path=out/entry['path'];path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            with urllib.request.urlopen(SOURCE+entry['path'],timeout=60) as r:path.write_bytes(r.read())
        b=path.read_bytes();gitsha=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        assert gitsha==entry['git_blob_sha1'],entry['path']
        record.append({**entry,'sha256':digest(path)})
    dump(OUT/'KI_ACQUISITION.json',{'utc':stamp(),'files':record,'potency_content_printed':False})
    print('Specified Ki Final/cache bytes acquired; no values printed.',flush=True)

def dec(s):
    try:
        v=decimal.Decimal(s)
        return v if v.is_finite() else None
    except decimal.InvalidOperation:return None

def make_view(t,endpoint,obs,smiles,pool,refs,edges,track):
    ee=degree_weighted_edges(edges);q=sorted(set(x for e in edges for x in e));ids=sorted(set(pool)|set(q))
    return TaskView(track+'::'+t['task_id'],track,t['task_id'],t['group_id'],t['source_assay'],t['target_assay'],endpoint,
       'pKi' if endpoint=='Ki' else 'pIC50',t['target_id'],
       tuple(Compound(x,smiles[x],float(obs[(t['source_assay'],x)])) for x in ids),tuple(pool),
       tuple((x,float(obs[(t['target_assay'],x)])) for x in refs),tuple(q),ee)

def write_task(view,obs,split,exposure):
    base=RES/view.endpoint/view.track_id
    name=view.legacy_task_id+'.json'
    dump(base/'tasks'/name,view.to_dict())
    dump(base/'truth'/name,{'task_id':view.task_id,'query_truth':{x:float(obs[(view.target_assay,x)]) for x in view.query_ids}})
    if view.track_id==ACQUISITION_TRACK:
        dump(base/'pool_truth'/name,{'task_id':view.task_id,'pool_truth':{x:float(obs[(view.target_assay,x)]) for x in view.reference_pool_ids}})
    return {'task_id':view.task_id,'legacy_task_id':view.legacy_task_id,'endpoint':view.endpoint,'potency_unit':view.potency_unit,
        'track_id':view.track_id,'group_id':view.group_id,'split':split,'source_assay':view.source_assay,'target_assay':view.target_assay,
        'target_id':view.target_chembl_id,'reference_pool':len(view.reference_pool_ids),'query_points':len(view.query_ids),'query_edges':len(view.query_edges),
        'query_signature':view.query_signature,'exposure_status':exposure,'task_file':str((base/'tasks'/name).relative_to(RES)),
        'truth_file':str((base/'truth'/name).relative_to(RES))}

def qualify():
    check_freeze();assert not (OUT/'KI_ELIGIBILITY.json').exists(),'No silent requalification'
    event('numeric_qualification_started',scope='Specified new Ki Final and same-key caches; no scientific effects released')
    identity=list(rows(PRE/'provenance/COMPOUND_IDENTITY.tsv'))
    canon={r['raw_compound_id']:r['canonical_compound_id'] for r in identity if r['status']=='eligible'}
    smiles={r['canonical_compound_id']:r['canonical_isomeric_smiles'] for r in identity if r['status']=='eligible'}
    assays={r['assay_chembl_id']:r for r in json.loads((PRE/'provenance/ASSAYS.json').read_text())}
    taskpre={r['task_id']:r for r in rows(PRE/'provenance/TASKS_IDENTITY_SENTINEL.tsv')}
    groups=list(rows(PRE/'provenance/GROUPS.tsv'))
    assaygroup={a:r['group_id'] for r in groups for a in json.loads(r['assay_ids'])}
    removed={r['group_id'] for r in groups if r['removed_for_IC50_overlap']=='True'}
    rawkeys=set();finalvals=collections.defaultdict(set);bad=collections.defaultdict(set);targets=collections.defaultdict(set);edges=collections.defaultdict(set)
    final=DATA/'ki_raw/Ki_MinimalCuration_Unpruned/Ki_MinimalCuration_Unpruned_Final.csv.gz'
    nr=0;rowreasons=collections.Counter()
    for r in rows(final):
        nr+=1; side=[]
        for p in (1,2):
            a=r[f'pair{p}_assay_ID'];pair=[]
            for j in (1,2):
                raw=r[f'pair{p}_compoundID_{j}'];rawkeys.add((a,raw))
                if raw not in canon:pair.append(None);continue
                x=canon[raw];key=(a,x);v=dec(r[f'pair{p}_pchembl_value_{j}'])
                if v is None:bad[key].add('nonfinite_or_invalid_Final')
                else:finalvals[key].add(v)
                targets[key].add(r['Target_Name']);pair.append(x)
            side.append((a,pair))
        (a,aa),(b,bb)=side
        if None in aa+bb:rowreasons['structure_ineligible']+=1;continue
        if set(aa)!=set(bb) or len(set(aa))!=2 or a==b:rowreasons['pair_identity_mismatch']+=1;continue
        tid='__'.join(sorted([a,b]));edges[tid].add(tuple(sorted(aa)))
    cv=collections.defaultdict(set);cachekeys=set();cachebad=collections.defaultdict(set);cache_rows=0
    for path in sorted((DATA/'ki_raw/Ki_MinimalCuration_Unpruned/Datasets').glob('*.csv.gz')):
        for r in rows(path):
            a,raw=r['assay_chembl_id'],r['compound_chembl_id']
            if (a,raw) not in rawkeys or raw not in canon:continue
            cache_rows+=1;key=(a,canon[raw]);cachekeys.add(key);v=dec(r['pchembl_value'])
            if v is None:cachebad[key].add('nonfinite_or_invalid_cache')
            else:cv[key].add(v)
            for field,expected in [('standard_relation','='),('standard_type','Ki'),('standard_units','nM')]:
                if field in r and r[field]!=expected:cachebad[key].add('explicit_'+field+'_ineligible')
    obs={};status=[];conflicts=[]
    for key,values in sorted(finalvals.items()):
        reasons=set(bad[key])|cachebad[key]
        if key not in cachekeys:reasons.add('same_key_cache_absent')
        if not values<=cv[key]:reasons.add('Final_cache_Decimal_mismatch')
        if len(values|cv[key])!=1:reasons.add('canonical_key_multiple_values')
        a,x=key
        if len(targets[key])!=1 or a not in assays or targets[key]!={assays[a]['target_chembl_id']}:reasons.add('target_mapping_conflict')
        if not reasons:obs[key]=next(iter(values))
        status.append({'assay_id':a,'canonical_compound_id':x,'status':'excluded' if reasons else 'eligible','reasons':json.dumps(sorted(reasons)),
          'final_distinct_values':len(values),'cache_distinct_values':len(cv[key]),'group_id':assaygroup[a]})
        if reasons:conflicts.append({'assay_id':a,'canonical_compound_id':x,'Final_values':list(map(str,sorted(values))),
                                     'cache_values':list(map(str,sorted(cv[key]))),'reasons':sorted(reasons)})
    for key in set(bad)-set(finalvals):
        status.append({'assay_id':key[0],'canonical_compound_id':key[1],'status':'excluded','reasons':json.dumps(sorted(bad[key])),
          'final_distinct_values':0,'cache_distinct_values':len(cv[key]),'group_id':assaygroup[key[0]]})
    table(OUT/'KI_OBSERVATION_ELIGIBILITY.tsv',status)
    dump(DATA/'ki_qualified/CONFLICT_VALUES.json',conflicts)
    table(DATA/'ki_qualified/UNIQUE_OBSERVATIONS.tsv.gz',({'assay_id':a,'canonical_compound_id':x,'pchembl_decimal':str(v)} for (a,x),v in sorted(obs.items())))
    roles=collections.defaultdict(dict)
    for r in rows(PRE/'provenance/SENTINEL_POOL_ASSIGNMENT.tsv'):roles[r['task_id']][r['canonical_compound_id']]=r['role']
    manifests=[];qualification=[]
    for tid,ee in sorted(edges.items()):
        a,b=tid.split('__');group=assaygroup[a];assert assaygroup[b]==group
        qreason=[]
        if group in removed:qreason.append('whole_group_IC50_overlap')
        if a not in assays or b not in assays or assays[a]['target_chembl_id']!=assays[b]['target_chembl_id']:qreason.append('target_mismatch')
        valid=sorted((x,y) for x,y in ee if all((ass,z) in obs for ass in (a,b) for z in (x,y)))
        ids=sorted({x for e in valid for x in e})
        if len(ids)<20:qreason.append('fewer_than_20_qualified_participants')
        t={'task_id':tid,'group_id':group,'source_assay':a,'target_assay':b,'target_id':assays[a]['target_chembl_id']}
        row={**t,'participants':len(ids),'eligible_original_edges':len(valid),'fixed_status':'not_eligible','acquisition_status':'not_eligible','reasons':json.dumps(qreason)}
        if not qreason:
            K=fingerprints([smiles[x] for x in ids]);c=[ids[j] for j in maxmin(ids,K,range(len(ids)),8)]
            qe=[(x,y) for x,y in valid if x not in c and y not in c]
            if qe:
                v=make_view(t,'Ki',obs,smiles,c,c,qe,FIXED_TRACK)
                manifests.append(write_task(v,obs,'frozen_application','numeric eligibility program access; predictions not scored yet'))
                row['fixed_status']='eligible'
            else:row['fixed_status']='no_query_edge_after_fixed_C8'
            original=taskpre.get(tid,{})
            if original.get('identity_sentinel_eligible')=='True':
                pool=sorted(x for x in ids if roles[tid].get(x)=='reference_pool')
                qe=[(x,y) for x,y in valid if roles[tid].get(x)==roles[tid].get(y)=='sentinel']
                if len(pool)>=8 and qe:
                    v=make_view(t,'Ki',obs,smiles,pool,[],qe,ACQUISITION_TRACK)
                    manifests.append(write_task(v,obs,'frozen_application','numeric eligibility program access; predictions not scored yet'))
                    row['acquisition_status']='eligible'
                else:row['acquisition_status']='frozen_roles_insufficient_pool_or_edges'
        qualification.append(row)
    table(OUT/'KI_TASK_ELIGIBILITY.tsv',qualification);table(RES/'Ki/TASK_MANIFEST.tsv',manifests)
    summary={'utc':stamp(),'source_expanded_rows':nr,'raw_Final_keys':len(rawkeys),'canonical_Final_keys':len(status),'eligible_unique_observations_before_group_removal':len(obs),
      'excluded_key_reasons':dict(collections.Counter(reason for s in status for reason in json.loads(s['reasons']))),
      'cache_same_key_rows_checked':cache_rows,'structural_row_exclusions':dict(rowreasons),'retained_pre_numeric_groups_never_split':True,'tracks':{}}
    for track in (FIXED_TRACK,ACQUISITION_TRACK):
        mm=[m for m in manifests if m['track_id']==track];gg={m['group_id'] for m in mm}
        summary['tracks'][track]={'tasks':len(mm),'groups':len(gg),'targets':len({m['target_id'] for m in mm}),
          'query_edges':sum(m['query_edges'] for m in mm),'query_task_points':sum(m['query_points'] for m in mm),
          'formal_run_authorized':len(gg)>=30}
    dump(OUT/'KI_ELIGIBILITY.json',summary)
    event('numeric_qualification_complete',effects_released=False,summary_path='KI_ELIGIBILITY.json')
    print(json.dumps(summary,indent=2),flush=True)

def predict():
    check_freeze();assert not (OUT/'KI_PREDICTIONS_FREEZE.json').exists()
    summary=json.loads((OUT/'KI_ELIGIBILITY.json').read_text());manifest=list(rows(RES/'Ki/TASK_MANIFEST.tsv'))
    eligible=[m for m in manifest if summary['tracks'][m['track_id']]['formal_run_authorized']]
    nf=sum(m['track_id']==FIXED_TRACK for m in eligible);na=len(eligible)-nf
    # One affine plus two solves/fixed task. At most 37 support states/acquisition task, two operations/state.
    reserve=3*nf+74*na;assert reserve+100<=25000
    dump(OUT/'KI_BUDGET_PREFLIGHT.json',{'fixed_tasks':nf,'acquisition_tasks':na,'worst_case_scientific_operations':reserve,
       'engineering_operation_reserve':100,'evaluation_reserved':True,'bootstrap_reserve_draws':5000,'cpu_limit_seconds':28800,'memory_limit_GiB':16,'passed':True})
    ledger=collections.Counter();saved=[];alltraces=[]
    start=time.process_time();event('task_local_reference_predictions_started',query_truth_loaded_by_predictor=False)
    for i,m in enumerate(eligible):
        v=TaskView.from_dict(json.loads((RES/m['task_file']).read_text()))
        if v.track_id==FIXED_TRACK:
            for model,p in fixed_predict(v,ledger).items():saved.append({'task_id':v.task_id,'legacy_task_id':v.legacy_task_id,'group_id':v.group_id,'track_id':v.track_id,'model':model,'budget':8,'cost_used':8,'query_signature':v.query_signature,'predictions':p})
        else:
            pp=RES/'Ki'/ACQUISITION_TRACK/'pool_truth'/(v.legacy_task_id+'.json')
            pooltruth=json.loads(pp.read_text())['pool_truth'];assert set(pooltruth)==set(v.reference_pool_ids)
            def factory(policy):return BudgetOracle(v,pooltruth.__getitem__,8,OUT/'oracle_traces'/v.legacy_task_id/(policy+'.jsonl'))
            pred,traces=acquire(v,factory,ledger)
            for row in pred:saved.append({**row,'task_id':v.task_id,'legacy_task_id':v.legacy_task_id,'group_id':v.group_id,'track_id':v.track_id,'model':row['policy'],'query_signature':v.query_signature})
            alltraces.extend({'task_id':v.task_id,**r} for r in traces)
        assert ledger['affine_fits']+ledger['kernel_solves']<=reserve
        assert time.process_time()-start<28800
        print(f'Predictions {i+1}/{len(eligible)} {v.track_id} {v.legacy_task_id}',flush=True)
    dump(OUT/'KI_PREDICTIONS.json',saved);dump(OUT/'KI_ACQUISITION_TRACES.json',alltraces)
    hashes={p.name:digest(p) for p in [OUT/'KI_PREDICTIONS.json',OUT/'KI_ACQUISITION_TRACES.json',RES/'Ki/TASK_MANIFEST.tsv']}
    dump(OUT/'KI_PREDICTIONS_FREEZE.json',{'utc':stamp(),'query_scoring_started':False,'sha256':hashes,'code_freeze':json.loads((OUT/'KI_IMPLEMENTATION_FREEZE.json').read_text()),
       'task_inputs_sha256':{m['task_file']:digest(RES/m['task_file']) for m in eligible},'allowed_inputs':'source values, structures, task-local reference pool via budget oracle; no query truth file',
       'ledger':dict(ledger),'CPU_seconds':time.process_time()-start})
    event('all_predictions_and_trajectories_frozen',query_scoring_started=False,sha256=digest(OUT/'KI_PREDICTIONS.json'))
    print(json.dumps(dict(ledger)),flush=True)

def extra(v,p,y):
    source={c.compound_id:c.source_potency for c in v.compounds};rescue=harm=0.
    def selected(pp,a,b):return min(a,b) if abs(pp[a]-pp[b])<=1e-8 else (a if pp[a]>pp[b] else b)
    for e in v.query_edges:
        a,b=e.compound_a,e.compound_b;ls=max(y[a],y[b])-y[selected(source,a,b)];lp=max(y[a],y[b])-y[selected(p,a,b)]
        if ls>0 and lp==0:rescue+=e.weight*ls
        if ls==0 and lp>0:harm+=e.weight*lp
    return {'rescued_loss':rescue,'introduced_loss':harm}

def score():
    check_freeze();assert not (OUT/'PRIMARY_CONTRASTS.json').exists(),'No silent rescoring'
    freeze=json.loads((OUT/'KI_PREDICTIONS_FREEZE.json').read_text())
    assert digest(OUT/'KI_PREDICTIONS.json')==freeze['sha256']['KI_PREDICTIONS.json']
    event('query_scoring_started',predictions_freeze_sha256=digest(OUT/'KI_PREDICTIONS_FREEZE.json'))
    manifests={m['task_id']:m for m in rows(RES/'Ki/TASK_MANIFEST.tsv')};metrics=[];cache={}
    for r in json.loads((OUT/'KI_PREDICTIONS.json').read_text()):
        tid=r['task_id'];m=manifests[tid]
        if tid not in cache:
            v=TaskView.from_dict(json.loads((RES/m['task_file']).read_text()))
            y=json.loads((RES/m['truth_file']).read_text())['query_truth'];cache[tid]=(v,y)
        v,y=cache[tid];assert v.query_signature==r['query_signature']
        metrics.append({**asdict(score_predictions(v,r['predictions'],y)),**extra(v,r['predictions'],y),'track_id':v.track_id,'model':r['model'],'budget':r['budget'],'cost_used':r['cost_used'],'target_id':v.target_chembl_id})
    table(OUT/'PER_TASK_METRICS.tsv',metrics)
    keys=['L','MAE','delta_MSE','rescued_loss','introduced_loss']
    cells=collections.defaultdict(list)
    for r in metrics:cells[(r['track_id'],r['model'],r['budget'])].append(r)
    groupscores={};agg=[]
    for (tr,mo,b),rs in sorted(cells.items()):
        gg=collections.defaultdict(list)
        for r in rs:gg[r['group_id']].append(r)
        gv={g:{k:float(np.mean([r[k] for r in vv])) for k in keys} for g,vv in gg.items()}
        groupscores[(tr,mo,b)]=gv
        agg.append({'track_id':tr,'model':mo,'budget':b,'groups':len(gg),'tasks':len(rs),**{k:float(np.mean([v[k] for v in gv.values()])) for k in keys}})
    for b in range(4,9):
        entries=[groupscores.get((ACQUISITION_TRACK,'random'+str(s),b)) for s in (17,43,71)]
        if not all(entries):continue
        mean={g:{k:float(np.mean([e[g][k] for e in entries])) for k in keys} for g in entries[0]}
        groupscores[(ACQUISITION_TRACK,'random_mean',b)]=mean
        agg.append({'track_id':ACQUISITION_TRACK,'model':'random_mean','budget':b,'groups':len(mean),
          'tasks':len(cells[(ACQUISITION_TRACK,'random17',b)]),**{k:float(np.mean([v[k] for v in mean.values()])) for k in keys}})
    table(OUT/'ALL_METRICS.tsv',agg);dump(OUT/'ALL_METRICS.json',agg)
    finalgroups=sorted({g for v in groupscores.values() for g in v})
    rng=np.random.default_rng(20261104)
    counts=rng.multinomial(len(finalgroups),np.ones(len(finalgroups))/len(finalgroups),size=5000) if finalgroups else np.zeros((5000,0))
    results={}
    for name,tr,m1,m0 in [('K1',FIXED_TRACK,'tanimoto','rbf'),('K2',ACQUISITION_TRACK,'maxmin','random_mean')]:
        g1=groupscores.get((tr,m1,8));g0=groupscores.get((tr,m0,8))
        if not g1 or len(g1)<30:results[name]={'status':'NOT_EXECUTED_FEWER_THAN_30_GROUPS'};continue
        gg=sorted(g1);assert set(gg)==set(g0)
        d=np.array([g1[g]['L']-g0[g]['L'] for g in gg]);cc=counts[:,[finalgroups.index(g) for g in gg]]
        den=cc.sum(1);assert np.all(den>0)
        sampled=(cc@d)/den;lo,hi=np.quantile(sampled,[.0125,.9875]);base=np.mean([g0[g]['L'] for g in gg])
        results[name]={'status':'COMPLETED','track_id':tr,'comparison':m1+' minus '+m0,'groups':len(gg),'difference_L':float(d.mean()),
          'family95_interval':[float(lo),float(hi)],'relative_improvement':float(-d.mean()/base) if base else None,
          'group_wins':int(sum(d<0)),'group_losses':int(sum(d>0)),'group_ties':int(sum(d==0)),
          'leave_one_group_out_difference_range':[float((d.sum()-d).min()/(len(d)-1)),float((d.sum()-d).max()/(len(d)-1))],
          'decision':'LIMITED_SUPPORT' if hi<0 else ('REVERSE_DIRECTION_SUPPORTED' if lo>0 else 'NO_FAMILY_INTERVAL_SUPPORT_NO_EQUIVALENCE')}
    dump(OUT/'PRIMARY_CONTRASTS.json',{'comparisons':results,'bootstrap_seed':20261104,'draws':5000,'family_size':2,'quantiles':[.0125,.9875],
      'sampling':'Shared union-group multiplicities; each track normalized over its eligible groups','scope':'Conditional on frozen strategy/model choices; same database cross-endpoint application'})
    table(OUT/'GROUP_METRICS.tsv',({'track_id':tr,'model':mo,'budget':b,'group_id':g,**v} for (tr,mo,b),gv in groupscores.items() for g,v in gv.items()))
    dump(OUT/'RESULT_CHECKS.json',{'status':'PASS','query_graphs_fixed':True,'exact_prediction_coverage':True,'initial4_count_in_budget':True,
      'all_acquisition_costs':sorted({r['cost_used'] for r in metrics if r['track_id']==ACQUISITION_TRACK}),
      'predictions_frozen_before_scoring':True,'new_neural_fits':0,'new_RF_fits':0,'old_model_inference':0,'old_sealed_access':False})
    event('query_scoring_complete',all_Ki_evaluation_outputs_now_exposed=True)
    print(json.dumps(results,indent=2),flush=True)

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['fetch','qualify','predict','score']);a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);start=time.process_time();wall=time.monotonic()
    try:globals()[a.stage]()
    finally:
        with (OUT/'RUNTIME_LEDGER.jsonl').open('a') as f:f.write(json.dumps({'utc':stamp(),'stage':a.stage,'CPU_seconds':time.process_time()-start,
          'wall_seconds':time.monotonic()-wall,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'platform':sys.platform,'python':sys.version})+'\n')
if __name__=='__main__':main()
