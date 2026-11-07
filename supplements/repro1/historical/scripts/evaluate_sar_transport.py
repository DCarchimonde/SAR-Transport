"""Select on validation, freeze, predict from C8, then separately score query values."""
from __future__ import annotations
import argparse,collections as co,csv,datetime as dt,gzip,json,pathlib,pickle,time,resource
import joblib
from decimal import Decimal
import numpy as np
import torch
import sar_transport_methods as methods
from sar_transport_methods import *
from sar_transport_eligibility import dump,sha,tsv,VAULT
from run_sar_transport import RUNS,GRID,check_freeze,save_predictions,duration

NAMES=['source','affine','rbf','tanimoto','rf','fscap','actfound','no_condition','conditional','isotonic']
SEEDS=[17,43]

def load_predictions(p,ts):
    values={}
    with open(p) as f:
        for r in csv.DictReader(f,delimiter='\t'):values[(r['task_id'],r['canonical_compound_id'])]=float(r['prediction'])
    return [np.array([values[(t['task_id'],t['ids'][j])] for j in t['q']]) for t in ts]

def selection(data):
    if (OUT/'EVALUATION_FREEZE.json').exists():raise RuntimeError('Evaluation already frozen; do not reselect')
    start=time.perf_counter();cpu=time.process_time();val=[t for t in data['tasks'] if t['split']=='validation']
    cal=json.loads((RUNS/'calibration/SELECTED.json').read_text());al=cal['affine_lambda'];results={};files={};kernel_fits=0
    for name in ['source','affine','rbf','tanimoto','isotonic']:
        p=RUNS/'calibration'/f'validation_{name}.tsv';pp=load_predictions(p,val)
        results[name]={'deterministic':aggregate([per_task_metrics(t,x) for t,x in zip(val,pp)])};files[str(p.relative_to(ROOT))]=sha(p)
    for name in ['rf','fscap','actfound','no_condition','conditional']:
        results[name]={}
        for seed in SEEDS:
            d=RUNS/f'{name}_seed{seed}';assert (d/'COMPLETE.json').exists(),d
            p=d/('validation_predictions.tsv' if name=='rf' else 'validation_best.tsv')
            pp=load_predictions(p,val);results[name][str(seed)]=aggregate([per_task_metrics(t,x) for t,x in zip(val,pp)])
            for f in [d/'COMPLETE.json',d/('model.joblib' if name=='rf' else 'best.pt'),p]:files[str(f.relative_to(ROOT))]=sha(f)
    grids={}
    for name in ['no_condition','conditional']:
        grids[name]={str(lam):{} for lam in GRID}
        for seed in SEEDS:
            model=build_model(name,data,seed);state=torch.load(RUNS/f'{name}_seed{seed}'/'best.pt',weights_only=False,map_location='cpu');model.load_state_dict(state['model'])
            for lam in GRID:
                if lam==1.:
                    pp=load_predictions(RUNS/f'{name}_seed{seed}'/'validation_best.tsv',val)
                else:
                    pp=[neural_predict(model,name,prediction_view(t),data,al,lam)[0] for t in val];kernel_fits+=len(val)
                met=aggregate([per_task_metrics(t,p) for t,p in zip(val,pp)]);grids[name][str(lam)][str(seed)]=met
                save_predictions(RUNS/f'{name}_seed{seed}'/f'validation_lambda{lam}.tsv',val,pp)
        lam=min(GRID,key=lambda l:(np.mean([grids[name][str(l)][str(s)]['L'] for s in SEEDS]),GRID.index(l)))
        cal[name]={'lambda':lam};results[name]=grids[name][str(lam)]
    averaged={n:{k:float(np.mean([m[k] for m in ss.values()])) for k in ['L','MAE','delta_MSE']} for n,ss in results.items()}
    no_struct=min(['source','affine','rbf'],key=lambda n:(averaged[n]['L'],NAMES.index(n)))
    comparator=min(NAMES[:8],key=lambda n:(averaged[n]['L'],NAMES.index(n)))
    # Descriptive development/validation diagnostics. Quantiles are development-only.
    diagnostics=[];desc_fits=0
    cached_val={kind:{t['task_id']:p for t,p in zip(val,load_predictions(RUNS/'calibration'/f'validation_{kind}.tsv',val))} for kind in ['rbf','tanimoto']}
    for t in data['tasks']:
        if t['split']=='evaluation':continue
        v=prediction_view(t);fp=data['fps'][t['mol_idx']];g,_,_=affine(t['xs'],v['yc'],t['c8'],al,data['mean'],data['sd'])
        cp=cal['rbf'];dp=cal['tanimoto']
        if t['split']=='validation':pc,pd=cached_val['rbf'][t['task_id']],cached_val['tanimoto'][t['task_id']]
        else:
            pc=simple_predict('rbf',v,fp,data['mean'],data['sd'],al,cp['lambda'],cp['bandwidth'])
            pd=simple_predict('tanimoto',v,fp,data['mean'],data['sd'],al,dp['lambda']);desc_fits+=2
        cmet=per_task_metrics(t,pc);dmet=per_task_metrics(t,pd)
        diagnostics.append({'task_id':t['task_id'],'group_id':t['group_id'],'split':t['split'],
          'mean_query_bridge_similarity':float(tanimoto(fp[t['q']],fp[t['c8']]).max(1).mean()),
          'context_residual_sd':float(np.std(v['yc']-g[t['c8']])),
          'D_minus_C_L':dmet['L']-cmet['L']})
    tsv(OUT/'DEVELOPMENT_CONDITION_DIAGNOSTICS.tsv',list(diagnostics[0]),diagnostics)
    strata=[]
    for f in ['mean_query_bridge_similarity','context_residual_sd']:
        threshold=float(np.median([r[f] for r in diagnostics if r['split']=='development']))
        for sp in ['development','validation']:
            for high in [False,True]:
                rr=[r for r in diagnostics if r['split']==sp and (r[f]>=threshold)==high];gg=co.defaultdict(list)
                for r in rr:gg[r['group_id']].append(r['D_minus_C_L'])
                strata.append({'feature':f,'development_median':threshold,'split':sp,'stratum':'high' if high else 'low',
                               'tasks':len(rr),'groups':len(gg),'group_equal_D_minus_C_L':float(np.mean([np.mean(v) for v in gg.values()])) if gg else None})
    dump(OUT/'DEVELOPMENT_CONDITION_STRATA.json',strata)
    rules={'status':'RULE_UNAVAILABLE','frozen_before_evaluation':True,'rules':[],
      'reason':'No specific direction-predicting chemical/experimental-condition rule established. The available bridge-similarity/residual-heterogeneity stratifications are generic local-interpolation diagnostics, not a distinct transferable chemical claim. No outcome-driven fragment or threshold search was performed.',
      'evidence':['DEVELOPMENT_CONDITION_DIAGNOSTICS.tsv','DEVELOPMENT_CONDITION_STRATA.json'],'rule_fits':0,'descriptive_development_medians':2}
    dump(OUT/'RULES.json',rules)
    selected={'common_affine_lambda':al,'calibration':cal,'best_no_structure':no_struct,'best_formal_comparator':comparator,
      'validation_per_seed':results,'validation_seed_mean_metrics':averaged,'neural_kernel_grids':grids,'isotonic_auxiliary_only':True}
    dump(OUT/'VALIDATION_SELECTION.json',selected)
    files.update({str(p.relative_to(ROOT)):sha(p) for p in [OUT/'VALIDATION_SELECTION.json',OUT/'RULES.json',pathlib.Path(__file__),
      ROOT/'scripts/sar_transport_methods.py',ROOT/'scripts/report_sar_transport.py',ROOT/'scripts/check_sar_transport_evaluation.py',ROOT/'scripts/sar_transport_strict.py',
      ROOT/'protocol/SAR_TRANSPORT_v1.json',ROOT/'protocol/SAR_TRANSPORT_STRICT_STRUCTURE_v1.json',OUT/'ELIGIBILITY_MANIFEST.json',OUT/'FEATURES.json',OUT/'C8.tsv',OUT/'QUERY_POINTS.tsv',
      OUT/'QUERY_EDGES.tsv.gz',DATA/'features/dataset.pkl',DATA/'features/preprocessors.pkl']})
    freeze={'frozen_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'files':files,'selected':{'best_no_structure':no_struct,'best_formal_comparator':comparator},
       'rules':'RULE_UNAVAILABLE','evaluation_C8_opened':False,'evaluation_queries_scored':False,'eligible_evaluation_groups':92,
       'extra_offline_kernel_fits':kernel_fits+desc_fits,'extra_offline_affine_fits':methods.AFFINE_MISSES,
       'descriptive_threshold_estimates':2,**duration(start,cpu)}
    dump(OUT/'EVALUATION_FREEZE.json',freeze)
    print(json.dumps({'stage':'frozen','best_no_structure':no_struct,'best_formal_comparator':comparator,'validation':averaged,'evaluation_opened':False}),flush=True)

def verify_eval_freeze():
    f=json.loads((OUT/'EVALUATION_FREEZE.json').read_text())
    for p,h in f['files'].items():assert sha(ROOT/p)==h,p
    return f

def project_evaluation_inputs(data):
    """Trusted projection: parses potency only for needed source keys / target C8 keys.

    The resulting predictor payload has no target-query label field. Shared assays
    can serve different roles in distinct tasks, but inputs remain task-local.
    """
    ts=[t for t in data['tasks'] if t['split']=='evaluation'];needed=set();source_keys=set();context_keys=set();query_target_keys=set()
    for t in ts:
        source_keys.update((t['source_assay'],c) for c in t['ids'])
        context_keys.update((t['target_assay'],t['ids'][j]) for j in t['c8'])
        query_target_keys.update((t['target_assay'],t['ids'][j]) for j in t['q'])
    needed=source_keys|context_keys
    allowed={}
    with gzip.open(VAULT/'UNIQUE_OBSERVATIONS.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=(r['assay_id'],r['canonical_compound_id'])
            if key in needed:allowed[key]=float(r['pchembl'])
    assert set(allowed)==needed
    dump(OUT/'EVALUATION_INPUT_ACCESS.json',{'stage':'frozen_prediction_inputs','source_observation_keys':len(source_keys),
      'context_target_keys':len(context_keys),'task_query_target_keys':len(query_target_keys),
      'query_target_keys_also_authorized_in_other_task_roles':len(query_target_keys&needed),
      'task_local_query_targets_passed_to_predictors':False,'global_numeric_unseen_claim':False,
      'explanation':'Shared-assay observations can be source inputs or C8 labels in another task. Predictors receive only their own allowed source values and eight target references; no adaptation state is shared across tasks. Qualification already accessed values for eligibility.'})
    views=[]
    for t in ts:
        v={k:x for k,x in t.items() if k not in ('yt','yc','xs')}
        v['xs']=np.array([allowed[(t['source_assay'],c)] for c in t['ids']]);v['yc']=np.array([allowed[(t['target_assay'],t['ids'][j])] for j in t['c8']])
        assert 'yt' not in v and len(v['yc'])==8
        views.append(v)
    return views

def predict(data):
    verify_eval_freeze();dest=RUNS/'evaluation_predictions';dest.mkdir(exist_ok=True)
    freeze_hash=sha(OUT/'EVALUATION_FREEZE.json')
    if (dest/'PREDICTIONS_FROZEN.json').exists():raise RuntimeError('Predictions already frozen; do not rerun adaptation')
    start=time.perf_counter();cpu=time.process_time();views=project_evaluation_inputs(data)
    with open(dest/'C8_ONLY_INPUTS.pkl','wb') as f:pickle.dump(views,f,protocol=5)
    selection=json.loads((OUT/'VALIDATION_SELECTION.json').read_text());cal=selection['calibration'];al=selection['common_affine_lambda']
    outputs={};adapts=[];fits=co.Counter()
    for name in ['source','affine','rbf','tanimoto','isotonic']:
        params=cal.get(name,{});pp=[]
        for v in views:
            pp.append(simple_predict(name,v,data['fps'][v['mol_idx']],data['mean'],data['sd'],al,params.get('lambda',1),params.get('bandwidth',1)))
            if name in ('rbf','tanimoto'):fits[name]+=1
            if name=='isotonic':fits[name]+=2
        path=dest/f'{name}.tsv';save_predictions(path,views,pp);outputs[name]={'deterministic':str(path.relative_to(ROOT))}
    for name in ['rf','fscap','actfound','no_condition','conditional']:
        outputs[name]={}
        for seed in SEEDS:
            d=RUNS/f'{name}_seed{seed}';pp=[]
            if name=='rf':
                model=joblib.load(d/'model.joblib')
                expected=np.random.RandomState(seed).randint(np.iinfo(np.int32).max,size=500)
                assert len(model.estimators_)==500 and np.array_equal([tree.random_state for tree in model.estimators_],expected)
                assert model.min_samples_leaf==5 and model.max_features==1/3
            else:
                model=build_model(name,data,seed);model.load_state_dict(torch.load(d/'best.pt',map_location='cpu',weights_only=False)['model'])
            for v in views:
                if name=='rf':p=model.predict(rf_features(v,data['fps'][v['mol_idx']],data['mean'],data['sd'],al));adapt=None
                elif name=='actfound':
                    cached=dest/'actfound_task_cache'/str(seed)/(v['task_id']+'.json');cached.parent.mkdir(parents=True,exist_ok=True)
                    started=cached.with_suffix('.started.json')
                    if cached.exists():
                        saved=json.loads(cached.read_text());assert saved['evaluation_freeze_sha256']==freeze_hash
                        p=np.array(saved['prediction']);adapt=saved['adaptation']
                    else:
                        if started.exists():raise RuntimeError('A task adaptation began without a completed cache; inspect before any repeat: '+str(started))
                        dump(started,{'evaluation_freeze_sha256':freeze_hash,'task_id':v['task_id'],'seed':seed})
                        p,adapt=neural_predict(model,name,v,data,al)
                        dump(cached,{'evaluation_freeze_sha256':freeze_hash,'prediction':p.tolist(),'adaptation':adapt})
                else:p,adapt=neural_predict(model,name,v,data,al,cal.get(name,{}).get('lambda',1))
                if adapt:adapts.append({'model':name,'seed':seed,'split':'evaluation','task_id':v['task_id'],'checkpoint':'fixed_best',**adapt})
                if name in ('no_condition','conditional'):fits[name]+=1
                pp.append(p)
            path=dest/f'{name}_{seed}.tsv';save_predictions(path,views,pp);outputs[name][str(seed)]=str(path.relative_to(ROOT))
            print(f'Predictions saved: {name} seed{seed}, {len(views)} tasks; query targets unscored',flush=True)
            del model
    fits['affine']=methods.AFFINE_MISSES
    tsv(dest/'TASK_ADAPTATIONS.tsv',list(adapts[0]),adapts)
    ledger={'frozen_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'prediction_files':outputs,
      'sha256':{p:sha(ROOT/p) for ss in outputs.values() for p in ss.values()},'input_sha256':sha(dest/'C8_ONLY_INPUTS.pkl'),
      'query_target_labels_passed_to_predictors':False,'target_values_per_task':8,'task_adaptations':len(adapts),
      'rf500_tree_seed_sequences_verified':True,
      'adaptation_steps_attempted':sum(r['steps_attempted'] for r in adapts),'adaptation_steps_applied':sum(r['steps_applied'] for r in adapts),
      'offline_fits':dict(fits),'total_offline_fits':sum(fits.values()),'new_global_neural_fits':0,'new_rf_fits':0,**duration(start,cpu)}
    dump(dest/'PREDICTIONS_FROZEN.json',ledger)
    print('All predictions hashed. Query scoring for this prediction batch has not started.',flush=True)

def load_eval_truth(data):
    ts=[t.copy() for t in data['tasks'] if t['split']=='evaluation'];needed={(t['target_assay'],c) for t in ts for c in t['ids']}
    needed.update((t['source_assay'],c) for t in ts for c in t['ids']);obs={};exact={}
    with gzip.open(VAULT/'UNIQUE_OBSERVATIONS.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            k=(r['assay_id'],r['canonical_compound_id'])
            if k in needed:obs[k]=float(r['pchembl']);exact[k]=Decimal(r['pchembl'])
    for t in ts:
        t['yt']=np.array([obs[(t['target_assay'],c)] for c in t['ids']]);t['xs']=np.array([obs[(t['source_assay'],c)] for c in t['ids']])
        t['yt_exact']=[exact[(t['target_assay'],c)] for c in t['ids']]
        t['xs_exact']=[exact[(t['source_assay'],c)] for c in t['ids']]
    return ts

def score(data):
    verify_eval_freeze();dest=RUNS/'evaluation_predictions';manifest=json.loads((dest/'PREDICTIONS_FROZEN.json').read_text())
    for p,h in manifest['sha256'].items():assert sha(ROOT/p)==h,p
    if (OUT/'PRIMARY_CONTRASTS.json').exists():raise RuntimeError('Evaluation already scored; do not overwrite')
    start=time.perf_counter();cpu=time.process_time();ts=load_eval_truth(data);by_model={};seed_metrics={};all_metrics={};edge_rows=[]
    for name,files in manifest['prediction_files'].items():
        by_model[name]={};seed_metrics[name]={}
        for seed,p in files.items():
            pp=load_predictions(ROOT/p,ts);mm=[]
            for t,pred in zip(ts,pp):
                mm.append(per_task_metrics(t,pred));pos={int(j):i for i,j in enumerate(t['q'])}
                for a,b,w in zip(t['edge_a'],t['edge_b'],t['weights']):
                    pd=pred[pos[int(a)]]-pred[pos[int(b)]];ds=t['xs'][a]-t['xs'][b];dy=t['yt'][a]-t['yt'][b]
                    dse=t['xs_exact'][a]-t['xs_exact'][b];dye=t['yt_exact'][a]-t['yt_exact'][b]
                    choosea=t['ids'][a]<t['ids'][b] if abs(pd)<=1e-8 else pd>0
                    L=max(t['yt'][a],t['yt'][b])-(t['yt'][a] if choosea else t['yt'][b])
                    edge_rows.append({'model':name,'seed':seed,'task_id':t['task_id'],'group_id':t['group_id'],'target_id':t['target_id'],
                     'source_unknown':t['source_unknown'],'compound_a':t['ids'][a],'compound_b':t['ids'][b],'source_gap':ds,'target_gap':dy,
                     'source_gap_exact':str(dse),'target_gap_exact':str(dye),
                     'predicted_target_gap':pd,'selected_compound':t['ids'][a] if choosea else t['ids'][b],'L':float(L),'edge_weight':float(w),
                     'observed_reversal_03':dse*dye<0 and abs(dse)>=Decimal('.3') and abs(dye)>=Decimal('.3'),
                     'observed_reversal_06':dse*dye<0 and abs(dse)>=Decimal('.6') and abs(dye)>=Decimal('.6')})
            by_model[name][seed]=mm;seed_metrics[name][seed]=aggregate(mm)
        all_metrics[name]={k:float(np.mean([v[k] for v in seed_metrics[name].values()])) for k in ['L','MAE','delta_MSE']}
    tsv(OUT/'ALL_EDGE_DECISIONS.tsv.gz',list(edge_rows[0]),edge_rows)
    task_rows=[{'model':name,'seed':seed,**r} for name,ss in by_model.items() for seed,rr in ss.items() for r in rr]
    tsv(OUT/'PER_TASK_METRICS.tsv',list(task_rows[0]),task_rows)
    dump(OUT/'ALL_METRICS.json',{'aggregation':'tasks equal per final group; groups equal; mean seed metrics, no prediction ensemble','models':all_metrics})
    dump(OUT/'PER_SEED_METRICS.json',seed_metrics)
    gids=sorted({t['group_id'] for t in ts});selection=json.loads((OUT/'VALIDATION_SELECTION.json').read_text())
    def groupvals(name,seed=None):
        ss=by_model[name]
        if seed is not None and seed in ss:ss={seed:ss[seed]}
        return np.array([np.mean([np.mean([r['L'] for r in rows if r['group_id']==g]) for rows in ss.values()]) for g in gids])
    rng=np.random.default_rng(20261031);draws=rng.integers(0,len(gids),(5000,len(gids)));contrasts=[];influence=[]
    for label,a,b in [('H1','tanimoto',selection['best_no_structure']),('model_increment','conditional',selection['best_formal_comparator']),('H2','conditional','no_condition')]:
        av,bv=groupvals(a),groupvals(b);delta=av-bv;boot=delta[draws].mean(1);low,high=np.percentile(boot,[100*.05/3/2,100*(1-.05/3/2)])
        ds={str(s):float((groupvals(a,str(s))-groupvals(b,str(s))).mean()) for s in SEEDS} if a=='conditional' else {}
        contrasts.append({'comparison':label,'first':a,'second':b,'groups':len(gids),'first_L':float(av.mean()),'second_L':float(bv.mean()),
          'difference_first_minus_second':float(delta.mean()),'family_95_interval':[float(low),float(high)],
          'individual_interval_level':1-.05/3,'relative_improvement':float((bv.mean()-av.mean())/bv.mean()) if bv.mean()!=0 else None,
          'per_seed_difference':ds,'both_seed_improvement':all(x<0 for x in ds.values()) if ds else None,
          'group_wins':int((delta< -1e-12).sum()),'group_losses':int((delta>1e-12).sum()),'group_ties':int((np.abs(delta)<=1e-12).sum())})
        for i,g in enumerate(gids):influence.append({'comparison':label,'omitted_group':g,'omitted_group_delta':float(delta[i]),'leave_one_group_mean_difference':float(np.delete(delta,i).mean())})
    dump(OUT/'PRIMARY_CONTRASTS.json',{'family_size':3,'bootstrap_draws':5000,'unit':'final group','scope':'fixed fitted model uncertainty only','comparisons':contrasts})
    tsv(OUT/'LEAVE_ONE_GROUP_INFLUENCE.tsv',list(influence[0]),influence)
    # Complete prespecified strata. Each stratum renormalizes available edge weights
    # within task, then averages task/group/seed; denominators explicitly reported.
    strata=[]
    funcs=[('source_gap',lambda r:'[0,.3)' if abs(Decimal(r['source_gap_exact']))<Decimal('.3') else '[.3,.6)' if abs(Decimal(r['source_gap_exact']))<Decimal('.6') else '[.6,1)' if abs(Decimal(r['source_gap_exact']))<1 else '[1,inf)'),
           ('target',lambda r:r['target_id']),('source_status',lambda r:'unknown' if r['source_unknown'] else 'known')]
    for typ,fn in funcs:
        buckets=co.defaultdict(list)
        for r in edge_rows:buckets[(r['model'],r['seed'],fn(r),r['group_id'],r['task_id'])].append((r['L'],r['edge_weight']))
        groupbucket=co.defaultdict(list)
        for (m,s,st,g,t),rr in buckets.items():groupbucket[(m,s,st,g)].append(sum(l*w for l,w in rr)/sum(w for l,w in rr))
        end=co.defaultdict(list)
        for (m,s,st,g),rr in groupbucket.items():end[(m,s,st)].append(float(np.mean(rr)))
        for (m,s,st),rr in end.items():strata.append({'stratification':typ,'model':m,'seed':s,'stratum':st,'groups':len(rr),'L':float(np.mean(rr))})
    tsv(OUT/'STRATIFIED_METRICS.tsv',list(strata[0]),strata)
    reversals=[]
    for t in ts:
        a,b=t['edge_a'],t['edge_b'];ds=t['xs'][a]-t['xs'][b];dy=t['yt'][a]-t['yt'][b]
        exact_ds=[t['xs_exact'][i]-t['xs_exact'][j] for i,j in zip(a,b)];exact_dy=[t['yt_exact'][i]-t['yt_exact'][j] for i,j in zip(a,b)]
        for margin in [.3,.6]:
            cutoff=Decimal(str(margin));eligible=np.array([abs(x)>=cutoff and abs(y)>=cutoff for x,y in zip(exact_ds,exact_dy)])
            flip=eligible&np.array([x*y<0 for x,y in zip(exact_ds,exact_dy)])
            reversals.append({'task_id':t['task_id'],'group_id':t['group_id'],'margin':margin,'query_edges':len(a),
              'margin_eligible_edges':int(eligible.sum()),'observed_reversals':int(flip.sum()),'unconditional_weighted_reversal_rate':float(t['weights']@flip),
              'conditional_weighted_reversal_rate':float(t['weights']@flip/(t['weights']@eligible)) if eligible.any() else None})
    tsv(OUT/'OBSERVED_REVERSALS.tsv',list(reversals[0]),reversals)
    h1,mi,h2=contrasts;model_continue=(mi['relative_improvement'] is not None and mi['relative_improvement']>=.10 and mi['family_95_interval'][1]<0 and h2['family_95_interval'][1]<0 and mi['both_seed_improvement'] and h2['both_seed_improvement'])
    h1supported=h1['family_95_interval'][1]<0
    trigger=(h1supported and h1['relative_improvement'] is not None and h1['relative_improvement']>=.10) or model_continue
    decision={'H1_supported_in_this_resource':h1supported,'H2_supported':h2['family_95_interval'][1]<0,
      'model_line_continue':model_continue,'knowledge_line_continue':False,'condition_rules':'RULE_UNAVAILABLE',
      'strict_structure_sensitivity_triggered':bool(trigger),'strict_structure_sensitivity_status':'ELIGIBILITY_CHECK_REQUIRED' if trigger else 'NOT_TRIGGERED',
      'decision':'MODEL_LINE_GATE_PASSED_REVIEW_ABSOLUTE_VALUE' if model_continue else 'CONTRACT_OR_STOP_NEW_CONTRIBUTION_NOT_ESTABLISHED',
      'JCIM_submission_readiness':False,'independent_database_confirmation':False,'causal_assay_mechanism':False}
    dump(OUT/'DECISION.json',decision)
    checks={'all_10_method_labels_scored':set(all_metrics)==set(NAMES),'prediction_hashes_verified':True,'same_query_lists':True,
      'query_target_labels_used_for_scoring_only_after_prediction_freeze':True,
      'task_local_target_query_labels_passed_to_predictors':False,
      'cross_task_authorized_source_or_C8_roles_reported':'EVALUATION_INPUT_ACCESS.json',
      'prediction_tie_rule':'absolute gap <=1e-8 canonical ID',
      'seed_predictions_ensembled':False,'statistical_units':len(gids),'old_sealed_outcomes_accessed':False,
      'primary_comparisons':3,'bootstrap_draws':5000,'new_fits_in_scoring':0,**duration(start,cpu)}
    dump(OUT/'RESULT_CHECKS.json',checks)
    print(json.dumps({'comparisons':contrasts,'decision':decision},indent=2),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['select','predict','score'],required=True);args=p.parse_args()
    check_freeze();torch.set_num_threads(4);data=load_data()
    {'select':selection,'predict':predict,'score':score}[args.stage](data)

if __name__=='__main__':main()
