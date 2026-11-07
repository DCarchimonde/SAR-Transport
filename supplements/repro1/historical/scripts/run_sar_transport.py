"""Bounded official SAR transport jobs; evaluation remains closed here."""
from __future__ import annotations
import argparse,collections as co,csv,gzip,hashlib,json,math,pathlib,pickle,resource,time,traceback,os,signal,threading
import joblib
import numpy as np
import torch
from sklearn.ensemble import RandomForestRegressor
from sar_transport_eligibility import ROOT,DATA,OUT,dump,sha,tsv
from sar_transport_methods import *
import sar_transport_methods as methods

RUNS=ROOT/'runs/sar_transport_v1';RUNS.mkdir(parents=True,exist_ok=True)
GRID=[.01,.1,1.,10.]

def check_freeze():
    f=OUT/'MODEL_IMPLEMENTATION_FREEZE_RECOVERY.json'
    if not f.exists():f=OUT/'MODEL_IMPLEMENTATION_FREEZE_JOB_LOCK.json'
    if not f.exists():f=OUT/'MODEL_IMPLEMENTATION_FREEZE_FSCAP_FIX.json'
    freeze=json.loads((f if f.exists() else OUT/'MODEL_IMPLEMENTATION_FREEZE.json').read_text())
    for p,h in freeze['files'].items():assert sha(ROOT/p)==h,p
    p=json.loads((OUT/'PROTOCOL_FREEZE.json').read_text());assert sha(ROOT/p['path'])==p['sha256']

def save_predictions(path,ts,preds):
    rows=[]
    for t,p in zip(ts,preds):
        assert len(p)==len(t['q']) and np.isfinite(p).all()
        rows.extend({'task_id':t['task_id'],'group_id':t['group_id'],'canonical_compound_id':t['ids'][int(q)],'prediction':float(x)} for q,x in zip(t['q'],p))
    tsv(path,['task_id','group_id','canonical_compound_id','prediction'],rows)

def duration(start,cpu):
    return {'wall_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}

def atomic_torch_save(value,path):
    temporary=path.with_suffix(path.suffix+'.pending');torch.save(value,temporary);os.replace(temporary,path)

def atomic_joblib_dump(value,path):
    temporary=path.with_suffix(path.suffix+'.pending');joblib.dump(value,temporary,compress=3);os.replace(temporary,path)

def start_budget_watch():
    from monitor_sar_transport_budget import snapshot
    finished=threading.Event()
    def watch():
        while not finished.wait(5):
            r=snapshot()
            if r['budget_charged_cpu_seconds']>48*3600-600 or r['aggregate_active_rss_bytes']>15.5*1024**3:
                dump(OUT/'RESOURCE_STOP.json',r)
                # This signal is sent from the fit's own process namespace.
                os.kill(os.getpid(),signal.SIGINT);return
    threading.Thread(target=watch,daemon=True).start();return finished

def calibrate(data):
    dest=RUNS/'calibration';dest.mkdir(exist_ok=True)
    if (dest/'COMPLETE.json').exists():return
    start=time.perf_counter();cpu=time.process_time();val=[t for t in data['tasks'] if t['split']=='validation']
    records=[];fits=co.Counter();saved={}
    def evaluate(kind,al,lam=1,bw=1):
        pp=[];mm=[]
        for t in val:
            p=simple_predict(kind,prediction_view(t),data['fps'][t['mol_idx']],data['mean'],data['sd'],al,lam,bw)
            pp.append(p);mm.append(per_task_metrics(t,p))
        score=aggregate(mm);records.append({'kind':kind,'affine_lambda':al,'lambda':lam,'bandwidth':bw,'group_metrics':score,'task_metrics':mm})
        return score,pp
    _,saved['source']=evaluate('source',1.)
    ag=[]
    for al in GRID:
        score,p=evaluate('affine',al);ag.append((score['L'],al,p));fits['affine']+=len(val)
    _,al,saved['affine']=min(ag,key=lambda x:(x[0],GRID.index(x[1])))
    # Reuse the selected affine C8 map; repetitions of its closed-form computation are cache-equivalent operations, not new fits.
    _,saved['isotonic']=evaluate('isotonic',al);fits['isotonic']+=2*len(val)
    selected={'affine_lambda':al}
    for kind in ['rbf','tanimoto']:
        configs=[]
        for lam in GRID:
            for bw in ([.5,1.,2.] if kind=='rbf' else [1.]):
                sc,p=evaluate(kind,al,lam,bw);configs.append((sc['L'],lam,bw,p));fits[kind]+=len(val)
        win=min(enumerate(configs),key=lambda x:(x[1][0],x[0]))[1]
        selected[kind]={'lambda':win[1],'bandwidth':win[2]};saved[kind]=win[3]
    metrics={}
    for name,p in saved.items():
        save_predictions(dest/f'validation_{name}.tsv',val,p);metrics[name]=aggregate([per_task_metrics(t,x) for t,x in zip(val,p)])
    dump(dest/'VALIDATION_GRID.json',records);dump(dest/'SELECTED.json',selected);dump(dest/'VALIDATION_METRICS.json',metrics)
    fits['affine']=methods.AFFINE_MISSES
    dump(dest/'COMPLETE.json',{'status':'COMPLETE','offline_fits':dict(fits),'total_offline_fits':sum(fits.values()),'evaluation_opened':False,**duration(start,cpu)})
    print(json.dumps({'stage':'calibration','selected':selected,'validation_metrics':metrics,'fits':dict(fits)}),flush=True)

def prepare_rf(data,al):
    dest=RUNS/'rf_features';dest.mkdir(exist_ok=True)
    if (dest/'COMPLETE.json').exists():return dest
    start=time.perf_counter();cpu=time.process_time();train=[t for t in data['tasks'] if t['split']=='development']
    rng=np.random.default_rng(20261031);episodes=[]
    for t in train:
        for i in range(8):
            c=rng.choice(len(t['ids']),8,replace=False);q=np.array(sorted(set(range(len(t['ids'])))-set(c)))
            episodes.append((t,c,q,i))
    n=sum(len(q) for t,c,q,i in episodes)
    first=episodes[0];width=rf_features(prediction_view(first[0],first[1],first[2]),data['fps'][first[0]['mol_idx']],data['mean'],data['sd'],al).shape[1]
    X=np.lib.format.open_memmap(dest/'X.npy',mode='w+',dtype='float32',shape=(n,width));y=np.zeros(n);w=np.zeros(n);at=0;index=[]
    for t,c,q,i in episodes:
        end=at+len(q);X[at:end]=rf_features(prediction_view(t,c,q),data['fps'][t['mol_idx']],data['mean'],data['sd'],al)
        y[at:end]=t['yt'][q];w[at:end]=1/(8*len(q))
        index.append({'task_id':t['task_id'],'episode':i,'context_ids':[t['ids'][j] for j in c],
                      'query_ids':[t['ids'][j] for j in q],'row_start':at,'row_end':end});at=end
    X.flush();np.save(dest/'y.npy',y);np.save(dest/'weights.npy',w)
    tsv(dest/'EPISODES.tsv',list(index[0]),index)
    dump(dest/'COMPLETE.json',{'status':'COMPLETE','rows':n,'columns':width,'tasks':len(train),'episodes':len(episodes),
      'offline_affine_fits':len(episodes),'shared_features_for_two_seeds':True,**duration(start,cpu)})
    return dest

def rf_job(data,seed):
    dest=RUNS/f'rf_seed{seed}';dest.mkdir(exist_ok=True)
    if (dest/'COMPLETE.json').exists():return
    al=json.loads((RUNS/'calibration/SELECTED.json').read_text())['affine_lambda'];featuredir=prepare_rf(data,al)
    start=time.perf_counter();cpu=time.process_time();before_affine=methods.AFFINE_MISSES;previous_cpu=previous_wall=0.;fit_calls=0
    checkpoint=dest/'forest_checkpoint.joblib'
    if checkpoint.exists():
        state=joblib.load(checkpoint);model=state['model'];previous_cpu=state['cpu_seconds'];previous_wall=state['wall_seconds'];fit_calls=state['incremental_fit_calls']
    else:
        if (dest/'STARTED.json').exists():raise RuntimeError('Prior RF attempt has no forest checkpoint; explicit retry record required')
        retry=seed==43 and (OUT/'RF_RETRY_AUTHORIZATION_20261103.json').exists() and RUNS.name=='sar_transport_v1'
        dump(dest/'STARTED.json',{'model':'rf','seed':seed,'global_fit_attempt':2 if retry else 1,'trees':500,'checkpoint_every_trees':50})
        model=RandomForestRegressor(n_estimators=50,min_samples_leaf=5,max_features=1/3,n_jobs=4,random_state=seed,warm_start=True)
    dump(dest/'CPU_CONTINUATION.json',{'prior_saved_cpu_seconds':previous_cpu})
    X=np.load(featuredir/'X.npy',mmap_mode='r');y=np.load(featuredir/'y.npy');w=np.load(featuredir/'weights.npy')
    done=len(getattr(model,'estimators_',[]))
    for trees in range(done+50,501,50):
        model.set_params(n_estimators=trees);model.fit(X,y,sample_weight=w);fit_calls+=1
        elapsed=duration(start,cpu);elapsed['cpu_seconds']+=previous_cpu;elapsed['wall_seconds']+=previous_wall
        atomic_joblib_dump({'model':model,'incremental_fit_calls':fit_calls,**elapsed},checkpoint)
        dump(dest/'FOREST_PROGRESS.json',{'trees_saved':trees,'target_trees':500,'incremental_fit_calls':fit_calls,**elapsed})
        print(f'RF seed{seed}: {trees}/500 trees saved; incremental fit calls={fit_calls}',flush=True)
    assert len(model.estimators_)==500
    atomic_joblib_dump(model,dest/'model.joblib')
    val=[t for t in data['tasks'] if t['split']=='validation'];pred=[];mm=[]
    for t in val:
        p=model.predict(rf_features(prediction_view(t),data['fps'][t['mol_idx']],data['mean'],data['sd'],al));pred.append(p);mm.append(per_task_metrics(t,p))
    save_predictions(dest/'validation_predictions.tsv',val,pred);dump(dest/'VALIDATION_METRICS.json',{'aggregate':aggregate(mm),'tasks':mm})
    elapsed=duration(start,cpu);elapsed['cpu_seconds']+=previous_cpu;elapsed['wall_seconds']+=previous_wall
    dump(dest/'COMPLETE.json',{'status':'COMPLETE','model':'rf','seed':seed,'global_fits':1,'trees':500,'incremental_fit_calls':fit_calls,'offline_affine_fits':methods.AFFINE_MISSES-before_affine,'evaluation_opened':False,**elapsed})
    print(json.dumps({'model':'rf','seed':seed,'validation':aggregate(mm),**duration(start,cpu)}),flush=True)

def validation(model,name,data,val,al,log,seed,epoch):
    preds=[];metrics=[];adapts=solves=0
    for t in val:
        p,adapt=neural_predict(model,name,prediction_view(t),data,al)
        preds.append(p);metrics.append(per_task_metrics(t,p))
        if adapt:
            log.write(json.dumps({'model':name,'seed':seed,'split':'validation','task_id':t['task_id'],'checkpoint':'epoch'+str(epoch),**adapt})+'\n');adapts+=1
        if name in ('no_condition','conditional'):solves+=1
    # Checkpoint selection is TASK equal MAE, not group equal L.
    return float(np.mean([m['MAE'] for m in metrics])),preds,metrics,adapts,solves

def neural_job(data,name,seed):
    dest=RUNS/f'{name}_seed{seed}';dest.mkdir(exist_ok=True)
    if (dest/'COMPLETE.json').exists():return
    al=json.loads((RUNS/'calibration/SELECTED.json').read_text())['affine_lambda']
    start=time.perf_counter();cpu=time.process_time();torch.set_num_threads(4);torch.manual_seed(seed);np.random.seed(seed)
    model=build_model(name,data,seed);lr=.00015 if name=='actfound' else .001
    opt=torch.optim.Adam(model.parameters(),lr=lr)
    scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(opt,T_max=150,eta_min=.0001) if name=='actfound' else None
    rng=np.random.default_rng(seed);train=[t for t in data['tasks'] if t['split']=='development'];val=[t for t in data['tasks'] if t['split']=='validation']
    best=float('inf');bestepoch=-1;firstepoch=0;stale=0;steps=0;trainadapt=valadapt=kernels=affineops=0;history=[]
    previous_cpu=previous_wall=0.
    if (dest/'last.pt').exists():
        state=torch.load(dest/'last.pt',weights_only=False,map_location='cpu');model.load_state_dict(state['model']);opt.load_state_dict(state['optimizer'])
        if scheduler:scheduler.load_state_dict(state['scheduler'])
        rng.bit_generator.state=state['numpy_rng'];torch.set_rng_state(state['torch_rng'])
        best=state['best'];bestepoch=state['best_epoch'];firstepoch=state['epoch']+1;stale=state['stale'];steps=state['optimizer_steps']
        trainadapt=state['train_adaptations'];valadapt=state['validation_adaptations'];kernels=state['episode_kernel_solves'];affineops=state['episode_affine_operations']
        history=state.get('history',json.loads((dest/'HISTORY.json').read_text()));assert len(history)==firstepoch
        previous_cpu=state['cpu_seconds'];previous_wall=state['wall_seconds']
    else:
        if (dest/'STARTED.json').exists():raise RuntimeError('Prior attempt has no resumable checkpoint; explicit engineering restart record required')
        dump(dest/'STARTED.json',{'model':name,'seed':seed,'global_fit_attempt':1,'parameter_count':sum(p.numel() for p in model.parameters() if p.requires_grad),'device':'cpu','threads':4})
    dump(dest/'CPU_CONTINUATION.json',{'prior_saved_cpu_seconds':previous_cpu,'resume_from_epoch_one_based':firstepoch+1})
    log=gzip.open(dest/'TASK_ADAPTATIONS.jsonl.gz','at')
    try:
        for epoch in range(firstepoch,150):
            epstart=time.perf_counter();model.train();total_loss=0.;order=rng.permutation(len(train))
            for b in range(0,len(order),8):
                opt.zero_grad(set_to_none=True);losses=[];fsviews=[];fstruth=[]
                for ix in order[b:b+8]:
                    t=train[int(ix)];cc=rng.choice(len(t['ids']),8,replace=False);rest=np.array(sorted(set(range(len(t['ids'])))-set(cc)))
                    qq=rng.choice(rest,min(16,len(rest)),replace=False);v=prediction_view(t,cc,qq)
                    if name=='fscap':
                        fsviews.append(v);fstruth.append(torch.tensor((t['yt'][qq]-data['mean'])/data['sd'],dtype=torch.float32));affineops+=1
                        continue
                    elif name=='actfound':
                        pred,loss,adapt=model.episode(v,data,epoch=epoch,query_targets=t['yt'][qq]);trainadapt+=1
                        log.write(json.dumps({'model':name,'seed':seed,'split':'development','task_id':t['task_id'],'checkpoint':'epoch'+str(epoch),**adapt})+'\n')
                    else:
                        pred=model.episode(v,data,al);truth=torch.tensor((t['yt'][qq]-data['mean'])/data['sd'],dtype=torch.float32)
                        loss=(pred-truth).square().mean();affineops+=1
                        if name in ('no_condition','conditional'):kernels+=1
                    if not torch.isfinite(loss):raise RuntimeError('Nonfinite training loss')
                    losses.append(loss)
                if name=='fscap':
                    losses=[(p-y).square().mean() for p,y in zip(model.batch_episodes(fsviews,data,al),fstruth)]
                loss=torch.stack(losses).mean()
                if name=='actfound' and epoch>=50:loss=loss+35*model.inner_lrs.square().sum()
                loss.backward()
                if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in model.parameters()):raise RuntimeError('Nonfinite gradients')
                opt.step();steps+=1;total_loss+=float(loss.detach())
            if scheduler:scheduler.step()
            mae,pp,mm,na,nk=validation(model,name,data,val,al,log,seed,epoch);valadapt+=na;kernels+=nk
            if name!='actfound':affineops+=len(val)
            if mae<best:
                best=mae;bestepoch=epoch;stale=0
                atomic_torch_save({'model':model.state_dict(),'epoch':epoch,'seed':seed,'name':name,'validation_MAE':mae},dest/'best.pt')
                save_predictions(dest/'validation_best.tsv',val,pp);dump(dest/'VALIDATION_METRICS.json',{'aggregate':aggregate(mm),'tasks':mm,'task_equal_MAE':mae})
            else:stale+=1
            row={'epoch':epoch,'task_equal_val_MAE':mae,'best_epoch':bestepoch,'mean_batch_training_loss':total_loss/math.ceil(len(train)/8),'epoch_wall_seconds':time.perf_counter()-epstart}
            history.append(row);dump(dest/'HISTORY.json',history);log.flush()
            elapsed=duration(start,cpu);elapsed['cpu_seconds']+=previous_cpu;elapsed['wall_seconds']+=previous_wall
            atomic_torch_save({'model':model.state_dict(),'optimizer':opt.state_dict(),'scheduler':scheduler.state_dict() if scheduler else None,'history':history,
              'epoch':epoch,'best':best,'best_epoch':bestepoch,'stale':stale,'numpy_rng':rng.bit_generator.state,'torch_rng':torch.get_rng_state(),
              'optimizer_steps':steps,'train_adaptations':trainadapt,'validation_adaptations':valadapt,'episode_kernel_solves':kernels,'episode_affine_operations':affineops,**elapsed},dest/'last.pt')
            print(f'{name} seed{seed} epoch{epoch+1}/150 valMAE={mae:.6f} best={bestepoch+1} seconds={row["epoch_wall_seconds"]:.2f}',flush=True)
            if (OUT/'RESOURCE_STOP.json').exists():raise RuntimeError('Main CPU budget reached')
            if stale>=20:break
        save_predictions(dest/'validation_last.tsv',val,pp)
        dump(dest/'COMPLETE.json',{'status':'COMPLETE','name':name,'seed':seed,'global_fits':1,'epochs':len(history),'best_epoch':bestepoch,
          'optimizer_steps':steps,'train_task_adaptations':trainadapt,'validation_task_adaptations':valadapt,
          'episode_kernel_solves':kernels,'episode_affine_operations':affineops,'offline_validation_affine_fits':methods.AFFINE_MISSES,'evaluation_opened':False,**elapsed})
    except BaseException as e:
        dump(dest/'INTERRUPTED.json',{'status':'INTERRUPTED','error':str(e),'traceback':traceback.format_exc(),**duration(start,cpu)})
        raise
    finally:log.close()

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['calibrate','rf','neural'],required=True)
    p.add_argument('--model',choices=['fscap','actfound','no_condition','conditional']);p.add_argument('--seed',type=int,choices=[17,43]);args=p.parse_args()
    check_freeze();torch.set_num_threads(4)
    key='calibration' if args.stage=='calibrate' else f'{args.model if args.stage=="neural" else "rf"}_seed{args.seed}'
    dest=RUNS/key;dest.mkdir(parents=True,exist_ok=True);lock=dest/'RUN_LOCK'
    while True:
        if (dest/'COMPLETE.json').exists():return
        if (OUT/'RESOURCE_STOP.json').exists():raise RuntimeError('Resource stop prevents opening another job')
        try:lock.mkdir();break
        except FileExistsError:
            if (dest/'INTERRUPTED.json').exists():raise RuntimeError('Existing attempt failed; no duplicate execution')
            time.sleep(5)
    finished=None
    try:
        data=load_data()
        finished=start_budget_watch()
        if args.stage=='calibrate':calibrate(data)
        elif args.stage=='rf':assert args.seed is not None;rf_job(data,args.seed)
        else:assert args.model and args.seed is not None;neural_job(data,args.model,args.seed)
    finally:
        if finished:finished.set()
        lock.rmdir()

if __name__=='__main__':main()
