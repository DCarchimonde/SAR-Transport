"""Render completed fixed evaluation results; no fitting or outcome selection."""
from __future__ import annotations
import csv,datetime as dt,gzip,json,math,collections as co
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sar_transport_eligibility import ROOT,OUT,dump,tsv,sha

RUNS=ROOT/'runs/sar_transport_v1'
LABELS={'source':'A  Source','affine':'B  Affine','rbf':'C  Potency RBF','tanimoto':'D  Tanimoto',
 'rf':'E  Random forest','fscap':'F  FS-CAP adapted','actfound':'G  ActFound adapted',
 'no_condition':'H  Deep kernel, no assay context','conditional':'I  Conditioned deep kernel','isotonic':'Auxiliary isotonic'}

def fit_ledger():
    strict=RUNS.name=='sar_transport_strict_v1'
    records=[];offline=0 if strict else 2;cpu=0.;global_neural=global_rf=episodes=kernelops=affineops=0
    for p in sorted(RUNS.rglob('COMPLETE.json')):
        r=json.loads(p.read_text());path=str(p.relative_to(ROOT));cpu+=r.get('cpu_seconds',0)
        typ='other';n=0
        if 'total_offline_fits' in r:n=r['total_offline_fits'];typ='calibration'
        elif 'offline_affine_fits' in r:n=r['offline_affine_fits'];typ='rf_features' if 'episodes' in r else 'rf'
        elif 'offline_validation_affine_fits' in r:n=r['offline_validation_affine_fits'];typ='neural'
        offline+=n
        if typ=='rf':global_rf+=r.get('global_fits',0)
        if typ=='neural':global_neural+=r.get('global_fits',0)
        kernelops+=r.get('episode_kernel_solves',0);affineops+=r.get('episode_affine_operations',0)
        episodes+=r.get('train_task_adaptations',0)+r.get('validation_task_adaptations',0)
        records.append({'path':path,'kind':typ,'offline_fits':n,'cpu_seconds':r.get('cpu_seconds',0),
          'global_fits':r.get('global_fits',0),'engineering_discarded':'engineering_attempts' in path})
    for name in ['ELIGIBILITY_SUMMARY.json','FEATURES.json']:
        r=json.loads((OUT/name).read_text());cpu+=r.get('cpu_seconds',0)
    freeze=json.loads((OUT/'EVALUATION_FREEZE.json').read_text())
    n=freeze['extra_offline_kernel_fits']+freeze['extra_offline_affine_fits'];offline+=n;cpu+=freeze.get('cpu_seconds',0)
    pred=json.loads((RUNS/'evaluation_predictions/PREDICTIONS_FROZEN.json').read_text())
    offline+=pred['total_offline_fits'];cpu+=pred.get('cpu_seconds',0);episodes+=pred['task_adaptations']
    check=json.loads((OUT/'RESULT_CHECKS.json').read_text());cpu+=check.get('cpu_seconds',0)
    records.extend([{'path':'EVALUATION_FREEZE.json','kind':'validation_selection','offline_fits':n,'cpu_seconds':freeze.get('cpu_seconds',0),'global_fits':0,'engineering_discarded':False},
                    {'path':'PREDICTIONS_FROZEN.json','kind':'evaluation_adaptation','offline_fits':pred['total_offline_fits'],'cpu_seconds':pred.get('cpu_seconds',0),'global_fits':0,'engineering_discarded':False}])
    steps=attempted=0
    for p in RUNS.rglob('TASK_ADAPTATIONS.jsonl.gz'):
        with gzip.open(p,'rt') as f:
            for line in f:
                r=json.loads(line);steps+=r['steps_applied'];attempted+=r['steps_attempted']
    steps+=pred['adaptation_steps_applied'];attempted+=pred['adaptation_steps_attempted']
    monitored=[]
    mp=OUT/'RESOURCE_MONITOR.jsonl'
    if mp.exists():
        with open(mp) as f:
            for line in f:
                try:monitored.append(json.loads(line))
                except json.JSONDecodeError:pass
    discarded=sum(r['global_fits'] for r in records if r['kind']=='neural' and r['engineering_discarded'])
    recovery_file=OUT/'RECOVERY_CPU_CHARGES.json';recovery=json.loads(recovery_file.read_text()) if recovery_file.exists() else {}
    known_interruption_cpu=recovery.get('known_interrupted_cpu_seconds',0.)
    cpu+=known_interruption_cpu
    conservative=recovery.get('conservative_unobserved_gap_charge_seconds',0.)
    rf_calls=sum(json.loads(p.read_text()).get('incremental_fit_calls',1) for p in RUNS.glob('rf_seed*/COMPLETE.json'))
    result={'global_neural_attempts_including_discarded':global_neural,'valid_global_neural_fits':global_neural-discarded,
      'engineering_restarts':discarded,'global_rf_fits':global_rf,
      'global_rf_attempts_including_interrupted':global_rf+recovery.get('failed_global_rf_attempts',0),
      'completed_rf_job_fit_api_calls':rf_calls,'failed_rf_fit_api_calls':recovery.get('failed_global_rf_attempts',0),
      'neural_checkpoint_continuations':recovery.get('neural_checkpoint_continuations',0),
      'offline_calibration_fits_including_synthetic_checks':offline,
      'synthetic_algebraic_calibration_checks':0 if strict else 2,'descriptive_rule_median_estimates':0 if strict else 2,'rule_predictors_fitted':0,
      'task_support_adaptations':episodes,'support_steps_attempted':attempted,'support_steps_applied':steps,
      'training_validation_episode_kernel_solves':kernelops,'training_validation_episode_affine_operations':affineops,
      'recorded_scientific_cpu_seconds':cpu,'recorded_scientific_cpu_hours':cpu/3600,
      'known_interrupted_cpu_seconds':known_interruption_cpu,'conservative_unobserved_gap_charge_seconds':conservative,
      'budget_charged_cpu_hours':(cpu+conservative)/3600,'conservative_gap_charge_is_measured_cpu':False,
      'uncheckpointed_adaptation_steps_lost_in_interruption':bool(recovery),
      'adaptation_step_count_scope':'Persisted completed/resumed and discarded-attempt logs; any unflushed work after the interrupted checkpoint has unknown step count and separately retained CPU cost.' if recovery else 'All persisted task logs',
      'gpu_hours':0,'paid_compute':False,'peak_observed_active_rss_gib':max((x['aggregate_active_rss_bytes'] for x in monitored),default=0)/1024**3,
      'cpu_accounting_scope':'Process timers for scientific preprocessing/fits/selection/prediction/scoring; includes discarded fit. Download/network wait, environment install and report rendering are separate from these timers.',
      'records':records}
    dump(OUT/'RUN_AND_BUDGET_LEDGER.json',result)
    return result

def figures(metrics,contrasts,seeds):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,ax=plt.subplots(figsize=(8,5.6));names=list(LABELS)
    values=[metrics[n]['L'] for n in names];yy=np.arange(len(names))
    ax.barh(yy,values,color=['#246a8c' if n=='conditional' else '#8fabb8' for n in names],height=.62)
    for y,n in zip(yy,names):
        vv=list(seeds[n].values())
        if len(vv)>1:ax.scatter([v['L'] for v in vv],[y-.12,y+.12],c='#263b42',s=18,zorder=3)
    ax.set_yticks(yy,[LABELS[n] for n in names]);ax.invert_yaxis();ax.set_xlabel('Group-equal selection loss L (pIC50 units; lower is better)')
    ax.set_title('Held-out analogue choice: all frozen methods');ax.grid(axis='x',alpha=.2);fig.tight_layout()
    fig.savefig(OUT/'FIGURE_1_ALL_METHODS.png');fig.savefig(OUT/'FIGURE_1_ALL_METHODS.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,3.5))
    for i,c in enumerate(contrasts):
        d=c['difference_first_minus_second'];lo,hi=c['family_95_interval']
        ax.errorbar(d,i,xerr=[[d-lo],[hi-d]],fmt='o',color='#246a8c',capsize=5)
    ax.axvline(0,color='#666',linestyle='--',linewidth=1)
    ax.set_yticks(range(3),[f"{c['comparison']}: {LABELS[c['first']]} minus {LABELS[c['second']]}" for c in contrasts])
    ax.invert_yaxis();ax.set_xlabel('Selection-loss difference (negative favors the first method)')
    ax.set_title('Three prespecified comparisons\nBonferroni family 95% intervals; fixed-model group bootstrap')
    fig.tight_layout();fig.savefig(OUT/'FIGURE_2_PRIMARY_CONTRASTS.png');fig.savefig(OUT/'FIGURE_2_PRIMARY_CONTRASTS.pdf');plt.close(fig)

def supplemental():
    with open(OUT/'PER_TASK_METRICS.tsv') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    tm=co.defaultdict(list)
    for r in rows:tm[(r['model'],r['task_id'])].append(float(r['L']))
    vals={k:float(np.mean(v)) for k,v in tm.items()};wins=[]
    for name in LABELS:
        tasks=sorted(t for n,t in vals if n==name);diff=np.array([vals[name,t]-vals['source',t] for t in tasks])
        wins.append({'model':name,'tasks':len(tasks),'wins_vs_source':int((diff< -1e-12).sum()),'losses_vs_source':int((diff>1e-12).sum()),'ties_vs_source':int((abs(diff)<=1e-12).sum())})
    tsv(OUT/'TASK_WIN_LOSS_VS_SOURCE.tsv',list(wins[0]),wins)
    # Descriptive case policy fixed before evaluation: highest task loss of I,
    # largest I-versus-frozen-comparator improvement/deterioration, then task ID.
    selected=json.loads((OUT/'VALIDATION_SELECTION.json').read_text());comp=selected['best_formal_comparator']
    tasks=sorted(t for n,t in vals if n=='conditional');cases=[]
    categories={'highest_conditional_loss':sorted(tasks,key=lambda t:(-vals['conditional',t],t))[:3],
      'largest_conditional_deterioration':sorted(tasks,key=lambda t:(-(vals['conditional',t]-vals[comp,t]),t))[:3],
      'largest_conditional_improvement':sorted(tasks,key=lambda t:(vals['conditional',t]-vals[comp,t],t))[:3]}
    for category,tt in categories.items():
        for task in tt:cases.append({'category':category,'task_id':task,'conditional_L':vals['conditional',task],
          'comparator':comp,'comparator_L':vals[comp,task],'difference':vals['conditional',task]-vals[comp,task],
          'interpretation':'Outcome-selected descriptive example; not independent rule evidence.'})
    tsv(OUT/'DESCRIPTIVE_CASES.tsv',list(cases[0]),cases)
    for c in cases:
        g=next(r for r in rows if r['task_id']==c['task_id']);c['group_id']=g['group_id'];c['target_id']=g['target_id']
    dump(OUT/'DESCRIPTIVE_CASES.json',cases)
    # Report empty source-unknown strata explicitly.
    unknown_tasks={r['task_id'] for r in rows if r['source_unknown']=='True'}
    dump(OUT/'SOURCE_UNKNOWN_LAYER.json',{'eligible_evaluation_tasks':len(unknown_tasks),
       'status':'NO_ELIGIBLE_TASKS' if not unknown_tasks else 'SEE_STRATIFIED_METRICS','no_independence_claim_from_missing_metadata':True})

def main():
    metrics=json.loads((OUT/'ALL_METRICS.json').read_text())['models']
    contrasts=json.loads((OUT/'PRIMARY_CONTRASTS.json').read_text())['comparisons']
    seeds=json.loads((OUT/'PER_SEED_METRICS.json').read_text())
    assert set(metrics)==set(LABELS)
    figures(metrics,contrasts,seeds);supplemental();budget=fit_ledger()
    rows=[{'method':LABELS[n],**metrics[n]} for n in LABELS]
    tsv(OUT/'ALL_METHODS_SUMMARY.tsv',list(rows[0]),rows)
    dump(OUT/'REPORT_ARTIFACTS.json',{'generated_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
      'figure_count':2,'new_fits':0,'new_bootstrap_draws':0,'budget_ledger':'RUN_AND_BUDGET_LEDGER.json',
      'figures':[{'path':p.name,'sha256':sha(p)} for p in sorted(OUT.glob('FIGURE_*.png'))]})
    print(json.dumps({'reported_models':len(metrics),'budget':{k:v for k,v in budget.items() if k!='records'}},indent=2))

if __name__=='__main__':main()
