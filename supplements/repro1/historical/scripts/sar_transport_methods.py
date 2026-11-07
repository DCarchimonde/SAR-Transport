"""Frozen SAR transport predictors. Prediction interfaces never require query labels."""
from __future__ import annotations
import importlib.util,json,pathlib,pickle,types
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from sklearn.isotonic import IsotonicRegression
from sar_transport_eligibility import ROOT,DATA,OUT

def load_data():
    with open(DATA/'features/dataset.pkl','rb') as f:return pickle.load(f)

AFFINE_CACHE={}
AFFINE_MISSES=0

def affine(xs, yc, c, lam, mean, sd, use_cache=True):
    global AFFINE_MISSES
    x=(np.asarray(xs)[c]-mean)/sd;y=(np.asarray(yc)-mean)/sd
    key=(tuple(x),tuple(y),lam,mean,sd)
    if use_cache and key in AFFINE_CACHE:a,b=AFFINE_CACHE[key]
    else:
        xx=x-x.mean();yy=y-y.mean();b=(xx@yy+lam)/(xx@xx+lam);a=y.mean()-b*x.mean()
        if use_cache:
            AFFINE_MISSES+=1;AFFINE_CACHE[key]=(float(a),float(b))
    return (a+b*(np.asarray(xs)-mean)/sd)*sd+mean, float(a),float(b)

def tanimoto(a,b):
    a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64)
    inter=a@b.T;return inter/np.maximum(a.sum(1)[:,None]+b.sum(1)[None,:]-inter,1e-12)

def simple_predict(kind,t,fp,mean,sd,aff_lam,lam=1.,bw=1.):
    c,q=t['c8'],t['q'];xs=t['xs'];yc=t['yc']
    if kind=='source':return xs[q].copy()
    g,_,_=affine(xs,yc,c,aff_lam,mean,sd)
    if kind=='affine':return g[q]
    if kind=='isotonic':
        fits=[]
        for inc in (True,False):
            reg=IsotonicRegression(increasing=inc,out_of_bounds='clip').fit(xs[c],yc)
            fits.append((float(np.mean((reg.predict(xs[c])-yc)**2)),reg))
        reg=min(enumerate(fits),key=lambda p:(p[1][0],p[0]))[1][1]
        return reg.predict(xs[q])
    if kind=='rbf':
        x=xs/sd;K=np.exp(-.5*((x[c,None]-x[None,c])/bw)**2);k=np.exp(-.5*((x[q,None]-x[None,c])/bw)**2)
    elif kind=='tanimoto':K=tanimoto(fp[c],fp[c]);k=tanimoto(fp[q],fp[c])
    else:raise ValueError(kind)
    return g[q]+k@np.linalg.solve(K+lam*np.eye(8),yc-g[c])

def prediction_view(t, c=None, q=None):
    c=t['c8'] if c is None else c;q=t['q'] if q is None else q
    # Query target values intentionally omitted from the predictor object.
    return {k:v for k,v in t.items() if k not in ('yt','yc','c8','q')} | {'c8':c,'q':q,'yc':t['yt'][c] if 'yt' in t and t['yt'] is not None else t['yc']}

def rf_features(v,fp,mean,sd,aff_lam):
    c,q=v['c8'],v['q'];g,_,_=affine(v['xs'],v['yc'],c,aff_lam,mean,sd)
    sim=tanimoto(fp[q],fp[c]);res=(v['yc']-g[c])/sd
    triples=[]
    for row in sim:
        order=sorted(range(8),key=lambda j:(-row[j],v['ids'][c[j]]))
        triples.append(np.array([(row[j],(v['xs'][c[j]]-mean)/sd,res[j]) for j in order]).reshape(-1))
    common=np.concatenate([v['common_target'],v['condition']])
    return np.concatenate([fp[q].astype(np.float32),((v['xs'][q]-mean)/sd)[:,None],np.tile(common,(len(q),1)),
                           np.array(triples),((g[q]-mean)/sd)[:,None]],axis=1).astype(np.float32)

def per_task_metrics(t,pred):
    q=t['q'];where={int(j):i for i,j in enumerate(q)}
    a=np.array([where[int(x)] for x in t['edge_a']]);b=np.array([where[int(x)] for x in t['edge_b']])
    yt=t['yt'][q];diff=pred[a]-pred[b]
    choose_a=diff>1e-8
    tie=np.abs(diff)<=1e-8
    choose_a=np.where(tie,np.array([t['ids'][j]<t['ids'][k] for j,k in zip(t['edge_a'],t['edge_b'])]),choose_a)
    loss=np.maximum(yt[a],yt[b])-np.where(choose_a,yt[a],yt[b])
    return {'task_id':t['task_id'],'group_id':t['group_id'],'target_id':t['target_id'],
            'source_unknown':t['source_unknown'],'L':float(loss@t['weights']),
            'MAE':float(np.mean(np.abs(pred-yt))),'delta_MSE':float(((diff-(yt[a]-yt[b]))**2)@t['weights'])}

def aggregate(metrics):
    groups={}
    for r in metrics:groups.setdefault(r['group_id'],[]).append(r)
    return {key:float(np.mean([np.mean([r[key] for r in rs]) for rs in groups.values()])) for key in ['L','MAE','delta_MSE']}

class DMPNN(nn.Module):
    def __init__(self,atomdim,bonddim):
        super().__init__();self.wi=nn.Linear(atomdim+bonddim,128,bias=False);self.wh=nn.Linear(128,128,bias=False);self.wo=nn.Linear(atomdim+128,128)
    def forward(self,graphs):
        aa=[];bb=[];ss=[];dd=[];rr=[];owner=[];no=bo=0
        for i,(a,b,s,d,r) in enumerate(graphs):
            aa.append(a);bb.append(b);ss.append(s+no);dd.append(d+no);rr.append(r+bo);owner.extend([i]*len(a));no+=len(a);bo+=len(b)
        device=self.wi.weight.device
        a=torch.as_tensor(np.concatenate(aa),device=device);b=torch.as_tensor(np.concatenate(bb),device=device)
        s=torch.as_tensor(np.concatenate(ss),device=device);d=torch.as_tensor(np.concatenate(dd),device=device);r=torch.as_tensor(np.concatenate(rr),device=device)
        own=torch.tensor(owner,device=device);h0=F.relu(self.wi(torch.cat([a[s],b],1)));h=h0
        for _ in range(2):
            ag=torch.zeros((len(a),128),device=device).index_add(0,d,h)
            h=F.relu(h0+self.wh(ag[s]-h[r]))
        ag=torch.zeros((len(a),128),device=device).index_add(0,d,h)
        atoms=F.relu(self.wo(torch.cat([a,ag],1)))
        pooled=torch.zeros((len(graphs),128),device=device).index_add(0,own,atoms)
        counts=torch.bincount(own,minlength=len(graphs)).to(atoms.dtype)
        return pooled/counts[:,None]

class DeepKernel(nn.Module):
    def __init__(self,data,conditioned):
        super().__init__();self.encoder=DMPNN(data['atom_dim'],data['bond_dim']);self.conditioned=conditioned
        self.target_emb=nn.Embedding(len(data['target_vocab'])+1,32)
        self.modulator=nn.Linear(data['condition_dim'],256)
        self.constant_condition=nn.Parameter(torch.zeros(data['condition_dim']))
        self.phi=nn.Linear(128+1+32+data['target_text_dim'],64)
        self.logalpha=nn.Parameter(torch.tensor(0.))
    def episode(self,v,data,aff_lam,kernel_lam=1.):
        inds=np.concatenate([v['c8'],v['q']]);mids=v['mol_idx'][inds]
        embed=self.encoder([data['graphs'][int(j)] for j in mids]);dev=embed.device
        cond=torch.as_tensor(v['condition'],device=dev) if self.conditioned else self.constant_condition
        gb=self.modulator(cond);gamma,beta=gb[:128],gb[128:]
        embed=embed*(1+.1*torch.tanh(gamma))+beta
        ys=torch.as_tensor((v['xs'][inds]-data['mean'])/data['sd'],dtype=torch.float32,device=dev)[:,None]
        target=self.target_emb(torch.tensor(v['target_index'],device=dev)).expand(len(inds),-1)
        name=torch.as_tensor(v['common_target'],device=dev).expand(len(inds),-1)
        phi=self.phi(torch.cat([embed,ys,target,name],1))/8.
        fps=data['fps'][mids];similarity=torch.as_tensor(tanimoto(fps,fps),dtype=torch.float32,device=dev)
        K=F.softplus(self.logalpha)*similarity+phi@phi.T
        g,_,_=affine(v['xs'],v['yc'],v['c8'],aff_lam,data['mean'],data['sd'],use_cache=not self.training)
        rc=torch.as_tensor((v['yc']-g[v['c8']])/data['sd'],dtype=torch.float32,device=dev)
        solution=torch.linalg.solve(K[:8,:8]+kernel_lam*torch.eye(8,device=dev),rc)
        return torch.as_tensor((g[v['q']]-data['mean'])/data['sd'],dtype=torch.float32,device=dev)+K[8:,:8]@solution

class MLP(nn.Module):
    def __init__(self,indim,outdim,width=512):
        super().__init__();layers=[]
        for _ in range(5):layers.extend([nn.Linear(indim,width),nn.BatchNorm1d(width),nn.ReLU()]);indim=width
        layers.append(nn.Linear(indim,outdim));self.layers=nn.Sequential(*layers)
    def forward(self,x):return self.layers(x)

class FSCAP(nn.Module):
    def __init__(self,data):
        super().__init__();extra=1+data['condition_dim']+data['target_text_dim']
        self.context=MLP(2048+extra+1,128);self.query=MLP(2048+extra,128)
        dims=[256,512,512,512,256,256,1];layers=[]
        for i,(a,b) in enumerate(zip(dims,dims[1:])):
            layers.append(nn.Linear(a,b))
            if i<len(dims)-2:
                layers.append(nn.ReLU())
                if i<3:layers.append(nn.Dropout(.1))
        self.predictor=nn.Sequential(*layers)
    def batch_episodes(self,views,data,aff_lam):
        """Normalize across mixed tasks, as in the official mini-batch implementation."""
        dev=next(self.parameters()).device;contexts=[];queries=[];baselines=[];sizes=[]
        for v in views:
            c,q=v['c8'],v['q'];fp=data['fps'][v['mol_idx']]
            yc=(v['yc']-data['mean'])/data['sd'];xs=(v['xs']-data['mean'])/data['sd']
            common=np.concatenate([v['common_target'],v['condition']])
            contexts.append(np.concatenate([fp[c]*yc[:,None],xs[c,None],np.tile(common,(8,1)),yc[:,None]],1))
            queries.append(np.concatenate([fp[q],xs[q,None],np.tile(common,(len(q),1))],1));sizes.append(len(q))
            g,_,_=affine(v['xs'],v['yc'],c,aff_lam,data['mean'],data['sd'],use_cache=False)
            baselines.append((g[q]-data['mean'])/data['sd'])
        ce=self.context(torch.as_tensor(np.concatenate(contexts),dtype=torch.float32,device=dev)).reshape(len(views),8,128).mean(1)
        qe=self.query(torch.as_tensor(np.concatenate(queries),dtype=torch.float32,device=dev))
        ce=torch.repeat_interleave(ce,torch.tensor(sizes,device=dev),dim=0)
        pred=self.predictor(torch.cat([ce,qe],1))[:,0]+torch.as_tensor(np.concatenate(baselines),dtype=torch.float32,device=dev)
        return list(torch.split(pred,sizes))
    def episode(self,v,data,aff_lam):
        dev=next(self.parameters()).device;c=v['c8'];q=v['q'];fp=data['fps'][v['mol_idx']]
        yc=(v['yc']-data['mean'])/data['sd'];xs=(v['xs']-data['mean'])/data['sd']
        common=np.concatenate([v['common_target'],v['condition']])
        cf=np.concatenate([fp[c]*yc[:,None],xs[c,None],np.tile(common,(8,1)),yc[:,None]],1)
        qf=np.concatenate([fp[q],xs[q,None],np.tile(common,(len(q),1))],1)
        context=self.context(torch.as_tensor(cf,dtype=torch.float32,device=dev)).mean(0,keepdim=True)
        query=self.query(torch.as_tensor(qf,dtype=torch.float32,device=dev))
        res=self.predictor(torch.cat([context.expand(len(q),-1),query],1))[:,0]
        g,_,_=affine(v['xs'],v['yc'],c,aff_lam,data['mean'],data['sd'],use_cache=not self.training)
        return res+torch.as_tensor((g[q]-data['mean'])/data['sd'],dtype=torch.float32,device=dev)

def official_actfound_module():
    path=ROOT/'vendor/sar_transport/ActFound/models/meta_neural_network_architectures.py'
    spec=importlib.util.spec_from_file_location('official_actfound_arch',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

class ActFound(nn.Module):
    """Official backbone and native objective; clean separation of adaptation and scoring.

    See vendor source for the original GPL-3.0 implementation. GPU-specific wrapper
    is replaced with device-neutral task orchestration; no pretrained weights.
    """
    def __init__(self,data,seed):
        super().__init__();width=2048+1+data['condition_dim']+data['target_text_dim']
        args=types.SimpleNamespace(hid_dim=2048,num_stages=2,per_step_bn_statistics=True,learnable_bn_gamma=True,
          learnable_bn_beta=True,enable_inner_loop_optimizable_bn_params=False,num_updates=5,norm_layer='batch_norm',
          input_celline=False,rng=np.random.RandomState(seed))
        self.net=official_actfound_module().FCNReLUNormNetworkQSAR((2,width),args,meta=True)
        self.inner_lrs=nn.Parameter(torch.full((5,),.001));self.temp=.2

    def features(self,v,data):
        inds=np.concatenate([v['c8'],v['q']]);fps=data['fps'][v['mol_idx'][inds]]
        common=np.concatenate([v['common_target'],v['condition']])
        x=np.concatenate([fps,((v['xs'][inds]-data['mean'])/data['sd'])[:,None],np.tile(common,(len(inds),1))],1)
        dev=next(self.parameters()).device
        return torch.as_tensor(x,dtype=torch.float32,device=dev),torch.as_tensor(fps,dtype=torch.float32,device=dev)

    @staticmethod
    def robust(a,b,idx):
        diff=(a-b).abs();err=torch.where(diff<=2,diff.square(),diff-2+4)
        return torch.gather(err,0,idx).mean()

    def weights_indices(self,e,fp):
        u=F.normalize(e,p=2,dim=1,eps=1e-8);cos=u@u.T/self.temp
        inter=fp@fp.T;sim=inter/(fp.sum(1)[:,None]+fp.sum(1)[None,:]-inter).clamp_min(1e-8)
        si=torch.topk(sim[:8,:8]-torch.eye(8,device=e.device),k=4,dim=0).indices
        qi=torch.topk(sim[:8,8:],k=4,dim=0).indices
        sw=torch.softmax(torch.gather(cos[:8,:8]-1e9*torch.eye(8,device=e.device),0,si),0)
        qw=torch.softmax(torch.gather(cos[:8,8:],0,qi),0)
        return sim,si,qi,sw,qw

    def target_objective(self,values,yc,yq,sim,qi,qw):
        cv,qv=values[:8],values[8:];preds=(torch.gather(yc[:,None]-(cv[:,None]-qv[None,:]),0,qi)*qw).sum(0)
        n=len(yq);ii=torch.topk(sim[8:,8:]-torch.eye(n,device=sim.device),k=n//2,dim=0).indices
        loss2=self.robust(qv[:,None]-qv[None,:],yq[:,None]-yq[None,:],ii)
        loss1=self.robust(cv[:,None]-qv[None,:],yc[:,None]-yq[None,:],qi)
        return (.25*loss2+.5*loss1+(preds-yq).square().mean()).clamp_min(1e-12).sqrt()

    def episode(self,v,data,epoch=0,query_targets=None):
        # query_targets is training-only. Validation/evaluation calls omit it.
        training=query_targets is not None;x,fp=self.features(v,data);dev=x.device
        yc=torch.as_tensor(v['yc'],dtype=torch.float32,device=dev)
        scale=torch.tensor(1.,device=dev) if training else yc.std(unbiased=True).clamp_min(.2)
        w=self.net.layer_dict['linear'].weights
        if not training:w=w.detach().clone().requires_grad_()
        losses=[];applied=0;support_losses=[];initial=None
        # Preserve official per-step batch statistics during validation independently per task.
        saved={n:p.detach().clone() for n,p in self.net.named_parameters() if not p.requires_grad} if not training else {}
        for step in range(5):
            with torch.set_grad_enabled(training):
                e,_=self.net(x,params={'layer_dict.linear.weights':w},num_step=step,training=True,backup_running_statistics=False)
            if not training:e=e.detach()
            sim,si,qi,sw,qw=self.weights_indices(e,fp)
            values=F.linear(e,w)[:,0];cv=values[:8]
            sp=(torch.gather(yc[:,None]-(cv[:,None]-cv[None,:])*scale,0,si)*sw).sum(0)
            sl=((sp-yc).square().mean()/scale.square()).clamp_min(1e-12).sqrt()
            value=float(sl.detach());support_losses.append(value)
            if initial is None:initial=value
            if value<10*initial:
                grad=torch.autograd.grad(sl,w,create_graph=training and epoch>10,retain_graph=training)[0]
                if training and epoch<=10:grad=grad.detach()
                w=w-self.inner_lrs[step]*grad
                applied+=1
                if not training:w=w.detach().requires_grad_()
            values=F.linear(e,w)[:,0]
            if training and (epoch<5 or step==4):
                yq=torch.as_tensor(query_targets,dtype=torch.float32,device=dev)
                tl=self.target_objective(values,yc,yq,sim,qi,qw)
                if epoch<5:
                    weight=max(.2-epoch*.2/5,.03/5) if step<4 else min(.2+epoch*4*.2/5,1-4*.03/5)
                    tl=tl*weight
                losses.append(tl)
        pred=(torch.gather(yc[:,None]-(values[:8,None]-values[None,8:])*scale,0,qi)*qw).sum(0)
        if not training:
            for n,p in self.net.named_parameters():
                if n in saved:p.data.copy_(saved[n])
        return pred,torch.stack(losses).sum() if losses else None,{'steps_attempted':5,'steps_applied':applied,'support_loss_start':support_losses[0],'support_loss_end':support_losses[-1]}

def build_model(name,data,seed):
    if name=='fscap':return FSCAP(data)
    if name=='actfound':return ActFound(data,seed)
    if name in ('no_condition','conditional'):return DeepKernel(data,name=='conditional')
    raise ValueError(name)

def neural_predict(model,name,v,data,aff_lam,kernel_lam=1.):
    model.eval()
    if name=='actfound':
        with torch.enable_grad():pred,_,adapt=model.episode(v,data)
        return pred.detach().cpu().numpy().astype('float64'),adapt
    with torch.no_grad():
        pred=model.episode(v,data,aff_lam,kernel_lam) if name in ('no_condition','conditional') else model.episode(v,data,aff_lam)
    return pred.cpu().numpy().astype('float64')*data['sd']+data['mean'],None
