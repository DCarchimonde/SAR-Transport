"""Optional NumPy/RDKit adapters; no neural dependencies. P0 formulas preserved."""
import hashlib
import numpy as np
from scipy.special import ndtr

MEAN=6.956174367283088
SD=1.2545393278285155
TAU=.3115195915190363
LAM=.1
POLICIES=('maxmin','facility','point_vopt','contrast_vopt','pairwise_kg','source_spread','random17','random43','random71')

def fingerprints(smiles):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    return np.array([DataStructs.BulkTanimotoSimilarity(f,fps) for f in fps])

def maxmin(ids,K,pool,n):
    c=[min(pool,key=lambda j:ids[j])]
    while len(c)<n:
        c.append(min((j for j in pool if j not in c),key=lambda j:(-float(np.min(1-K[j,c])),ids[j])))
    return c

def affine(xs,yc,c,ledger):
    xx=(xs[c]-MEAN)/SD; yy=(yc-MEAN)/SD
    xc=xx-xx.mean(); ycc=yy-yy.mean()
    b=(xc@ycc+1)/(xc@xc+1); a=yy.mean()-b*xx.mean()
    ledger['affine_fits']+=1
    return (a+b*(xs-MEAN)/SD)*SD+MEAN

def fixed_predict(view,ledger):
    ids=[x.compound_id for x in view.compounds]; ix={x:j for j,x in enumerate(ids)}
    xs=np.array([x.source_potency for x in view.compounds]); c=[ix[x] for x,y in view.reference_labels]
    yc=np.array([y for x,y in view.reference_labels]); q=[ix[x] for x in view.query_ids]
    K=fingerprints([x.canonical_smiles for x in view.compounds]); g=affine(xs,yc,c,ledger)
    result={'source':xs[q],'affine':g[q]}
    for model in ['rbf','tanimoto']:
        M=np.exp(-.5*((xs[:,None]/SD-xs[None,:]/SD)/2.)**2) if model=='rbf' else K
        result[model]=g[q]+M[np.ix_(q,c)]@np.linalg.solve(M[np.ix_(c,c)]+LAM*np.eye(8),yc-g[c])
        ledger['kernel_solves']+=1
    return {m:dict(zip(view.query_ids,map(float,p))) for m,p in result.items()}

def gain(m,q):
    m=np.abs(np.asarray(m));q=np.asarray(q);s=np.sqrt(np.maximum(q,0));z=np.minimum(m/np.maximum(s,1e-100),1e100)
    g=s*np.exp(-.5*z*z)/np.sqrt(2*np.pi)-m*ndtr(-z)
    return np.where(q>1e-20,np.maximum(g,0),0.)

def acquire(view,oracle_factory,ledger):
    """Returns predictions at 4..8. Provider is task-local and pool-only.

    Cache entries keyed by ordered support IDs are reused only after each
    independent oracle has paid for that identical support. No target truth is
    passed to the scoring-free acquisition function.
    """
    ids=[x.compound_id for x in view.compounds];index={x:j for j,x in enumerate(ids)}
    pool=np.array([index[x] for x in view.reference_pool_ids],int)
    query=np.array([index[x] for x in view.query_ids],int)
    ea=np.array([index[e.compound_a] for e in view.query_edges]);eb=np.array([index[e.compound_b] for e in view.query_edges]);w=np.array([e.weight for e in view.query_edges])
    xs=np.array([x.source_potency for x in view.compounds]);K=fingerprints([x.canonical_smiles for x in view.compounds])
    initial=tuple(maxmin(ids,K,pool,4));cache={};output=[];traces=[]
    def fit(c,oracle):
        c=tuple(c);yc=np.array([oracle.reveal(ids[j]) for j in c])
        if c in cache:
            ledger['state_cache_hits']+=1;return cache[c]
        g=affine(xs,yc,list(c),ledger)
        inverse=np.linalg.solve(K[np.ix_(c,c)]+LAM*np.eye(len(c)),np.eye(len(c)))
        ledger['kernel_solves']+=1
        pred=g+K[:,c]@inverse@(yc-g[list(c)])
        V=TAU**2*(K-K[:,c]@inverse@K[list(c),:]);V=.5*(V+V.T)
        if np.min(np.diag(V)) < -1e-8:raise ValueError('material negative working variance')
        cache[c]=(pred,V);return pred,V
    for policy in POLICIES:
        oracle=oracle_factory(policy);c=list(initial);trace=[]
        for j in c:oracle.reveal(ids[j])
        for count in range(4,9):
            pred,V=fit(c,oracle)
            output.append({'policy':policy,'budget':count,'cost_used':oracle.spent,'reference_ids':[ids[j] for j in c],
                           'predictions':dict(zip(view.query_ids,map(float,pred[query])))})
            if count==8:break
            remaining=np.array([j for j in pool if j not in c],int)
            if policy=='maxmin':score=np.min(1-K[np.ix_(remaining,c)],axis=1)
            elif policy=='facility':
                before=np.max(K[np.ix_(query,c)],axis=1)
                score=np.maximum(K[np.ix_(query,remaining)]-before[:,None],0).sum(0)
            elif policy=='source_spread':
                spread=max(float(np.ptp(xs[pool])),1e-12)
                score=np.min(.5*(1-K[np.ix_(remaining,c)])+.5*np.abs(xs[remaining,None]-xs[None,c])/spread,axis=1)
            elif policy.startswith('random'):
                seed=int(policy[6:]);h=hashlib.sha256((view.legacy_task_id+'|'+str(seed)).encode()).hexdigest()
                rng=np.random.default_rng(int(h[:16],16))
                order=rng.permutation([j for j in pool if j not in initial]);j=int(order[count-4]);score=None
            else:
                den=np.maximum(np.diag(V)[remaining],0)+TAU**2*LAM
                if policy=='point_vopt':score=(V[np.ix_(query,remaining)]**2/den).mean(0)
                else:
                    cross=V[np.ix_(ea,remaining)]-V[np.ix_(eb,remaining)];q=cross**2/den
                    edge_var=np.maximum(np.diag(V)[ea]+np.diag(V)[eb]-2*V[ea,eb],0)
                    if np.any(q>edge_var[:,None]+1e-8):raise ValueError('information variance exceeds contrast variance')
                    score=w@q if policy=='contrast_vopt' else w@gain((pred[ea]-pred[eb])[:,None],q)
            if score is not None:
                k=min(range(len(remaining)),key=lambda k:(-float(score[k]),ids[int(remaining[k])]))
                j=int(remaining[k])
            trace.append({'step':count+1,'selected_id':ids[j],'policy_score':None if score is None else float(score[list(remaining).index(j)])})
            oracle.reveal(ids[j]);c.append(j);ledger['acquisition_steps']+=1
        assert oracle.spent==8
        traces.append({'policy':policy,'initial4':[ids[j] for j in initial],'steps':trace,'query_signature':view.query_signature})
    return output,traces
