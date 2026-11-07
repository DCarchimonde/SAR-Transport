"""Build development-only preprocessing and task inputs; no evaluation target values."""
from __future__ import annotations
import csv,gzip,hashlib,json,pathlib,pickle,time,resource,re
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import OneHotEncoder
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sar_transport_eligibility import ROOT,DATA,OUT,VAULT,dump,sha

def readtsv(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))

def oh(x, choices):return [float(x==v) for v in choices]+[float(x not in choices)]

def atomfeat(a):
    return np.array(oh(a.GetAtomicNum(),list(range(1,101)))+oh(a.GetTotalDegree(),list(range(6)))+
        oh(a.GetFormalCharge(),[-2,-1,0,1,2])+oh(int(a.GetChiralTag()),list(range(4)))+
        oh(a.GetTotalNumHs(),list(range(5)))+oh(int(a.GetHybridization()),[2,3,4,6,7])+
        [float(a.GetIsAromatic()),a.GetMass()*.01],dtype=np.float32)

def graph(m):
    atoms=np.stack([atomfeat(a) for a in m.GetAtoms()]);src=[];dst=[];rev=[];bf=[]
    for b in m.GetBonds():
        f=oh(str(b.GetBondType()),['SINGLE','DOUBLE','TRIPLE','AROMATIC'])+[float(b.GetIsConjugated()),float(b.IsInRing())]+oh(int(b.GetStereo()),list(range(6)))
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();k=len(src)
        src.extend([i,j]);dst.extend([j,i]);rev.extend([k+1,k]);bf.extend([f,f])
    return (atoms,np.array(bf,dtype=np.float32).reshape(-1,14),np.array(src,dtype=np.int64),np.array(dst,dtype=np.int64),np.array(rev,dtype=np.int64))

def clean_text(s):
    # Do not use assay/document identifiers or URLs as predictive tokens.
    s=re.sub(r'CHEMBL\d+',' ',s or '',flags=re.I)
    return re.sub(r'https?://\S+|10\.\d{4,9}/\S+',' ',s)

def main():
    start=time.perf_counter();cpu=time.process_time()
    dest=DATA/'features';dest.mkdir(exist_ok=True)
    if (dest/'dataset.pkl').exists():raise RuntimeError('Features already exist; no silent refit')
    assert json.loads((OUT/'ELIGIBILITY_SUMMARY.json').read_text())['gate']['passed']
    tasks=[r for r in readtsv(OUT/'TASKS.tsv') if r['status']=='eligible']
    cmap={r['raw_compound_id']:r for r in readtsv(OUT/'COMPOUND_IDENTITY.tsv')}
    cs={};qs={};edges={};members={}
    for r in readtsv(OUT/'C8.tsv'):cs.setdefault(r['task_id'],[]).append(r['canonical_compound_id'])
    for r in readtsv(OUT/'QUERY_POINTS.tsv'):qs.setdefault(r['task_id'],[]).append(r['canonical_compound_id'])
    for r in readtsv(OUT/'QUERY_EDGES.tsv.gz'):edges.setdefault(r['task_id'],[]).append((r['compound_a'],r['compound_b'],float(r['edge_weight'])))
    for r in readtsv(OUT/'ALL_ELIGIBLE_EDGES.tsv.gz'):members.setdefault(r['task_id'],set()).update([r['compound_a'],r['compound_b']])
    ids=sorted({c for t in tasks for c in members[t['task_id']]});ci={c:i for i,c in enumerate(ids)}
    fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
    fps=[];graphs=[]
    for c in ids:
        mol=Chem.MolFromSmiles(cmap[c]['canonical_isomeric_smiles']);fp=fpgen.GetFingerprint(mol)
        ar=np.zeros(2048,dtype=np.uint8);DataStructs.ConvertToNumpyArray(fp,ar);fps.append(ar);graphs.append(graph(mol))
    fp=np.stack(fps)
    # Row filter occurs before parsing potency. Evaluation target values are not loaded.
    obs={};scale_values=[]
    with gzip.open(VAULT/'UNIQUE_OBSERVATIONS.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['split']=='evaluation':continue
            obs[(r['assay_id'],r['canonical_compound_id'])]=float(r['pchembl'])
            if r['split']=='development':scale_values.append(float(r['pchembl']))
    mean=float(np.mean(scale_values));sd=float(np.std(scale_values));assert sd>0
    metadata={r['assay_chembl_id']:r for r in json.loads((DATA/'metadata/ASSAYS.json').read_text())}
    targetmeta={r['target_chembl_id']:r for r in json.loads((DATA/'metadata/TARGETS.json').read_text())}
    aids=sorted({a for r in tasks for a in [r['source_assay'],r['target_assay']]});ai={a:i for i,a in enumerate(aids)}
    tids=sorted({json.loads(r['target_ids'])[0] for r in tasks});ti={t:i for i,t in enumerate(tids)}
    da=sorted({a for r in tasks if r['split']=='development' for a in [r['source_assay'],r['target_assay']]})
    dt=sorted({json.loads(r['target_ids'])[0] for r in tasks if r['split']=='development'})
    texts={a:clean_text(metadata.get(a,{}).get('description')) for a in aids}
    names={t:clean_text(targetmeta.get(t,{}).get('pref_name','')) for t in tids}
    vectorizer=TfidfVectorizer(lowercase=True,ngram_range=(1,2),max_features=20000,min_df=1,sublinear_tf=True)
    train_text=[texts[a] for a in da]+[names[t] for t in dt]
    matrix=vectorizer.fit_transform(train_text)
    dim=min(128,matrix.shape[0]-1,matrix.shape[1]-1)
    svd=TruncatedSVD(n_components=dim,random_state=20261031)
    svd.fit(matrix)
    desc=svd.transform(vectorizer.transform([texts[a] for a in aids])).astype('float32')
    target_text=svd.transform(vectorizer.transform([names[t] for t in tids])).astype('float32')
    fields=['assay_type','assay_organism','assay_tax_id','assay_strain','assay_tissue','assay_cell_type','assay_subcellular_fraction','bao_format']
    cat=[[str(metadata.get(a,{}).get(f) or 'MISSING') for f in fields] for a in aids]
    onehot=OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float32)
    onehot.fit([cat[ai[a]] for a in da]);structured=onehot.transform(cat)
    masks=np.array([[float(v=='MISSING') for v in row]+[float(not texts[a]),1.,1.,1.] for a,row in zip(aids,cat)],dtype=np.float32)
    # Last three masks: substrate, incubation time and temperature not supplied as standardized fields.
    assayfeatures=np.concatenate([desc,structured,masks],1)
    dev_target_vocab={t:i+1 for i,t in enumerate(dt)}
    dataset=[]
    for r in tasks:
        tid=r['task_id'];comp=sorted(members[tid]);local={c:i for i,c in enumerate(comp)};target=json.loads(r['target_ids'])[0]
        xs=yt=None
        if r['split']!='evaluation':
            xs=np.array([obs[(r['source_assay'],c)] for c in comp],dtype=np.float64)
            yt=np.array([obs[(r['target_assay'],c)] for c in comp],dtype=np.float64)
        e=edges[tid]
        dataset.append({'task_id':tid,'group_id':r['group_id'],'split':r['split'],'target_id':target,
          'target_name':targetmeta.get(target,{}).get('pref_name',''),'source_assay':r['source_assay'],'target_assay':r['target_assay'],
          'source_unknown':r['source_unknown']=='True','ids':comp,'mol_idx':np.array([ci[c] for c in comp]),
          'c8':np.array([local[c] for c in cs[tid]],dtype=np.int64),'q':np.array([local[c] for c in qs[tid]],dtype=np.int64),
          'edge_a':np.array([local[a] for a,b,w in e]),'edge_b':np.array([local[b] for a,b,w in e]),'weights':np.array([w for a,b,w in e]),
          'xs':xs,'yt':yt,'condition':np.concatenate([assayfeatures[ai[r['source_assay']]],assayfeatures[ai[r['target_assay']]]]),
          'common_target':target_text[ti[target]],'target_index':dev_target_vocab.get(target,0)})
    bundle={'tasks':dataset,'compound_ids':ids,'fps':fp,'graphs':graphs,'mean':mean,'sd':sd,'target_vocab':dev_target_vocab,
            'desc_dim':dim,'condition_dim':len(dataset[0]['condition']),'target_text_dim':target_text.shape[1],
            'atom_dim':graphs[0][0].shape[1],'bond_dim':14}
    with open(dest/'dataset.pkl','wb') as f:pickle.dump(bundle,f,protocol=5)
    with open(dest/'preprocessors.pkl','wb') as f:pickle.dump({'tfidf':vectorizer,'svd':svd,'onehot':onehot,'fields':fields},f,protocol=5)
    dump(OUT/'FEATURES.json',{'development_only_preprocessing':True,'preprocessing_fits':{'target_location_scale':1,'TFIDF':1,'SVD':1,'categorical_vocabulary':1},
      'evaluation_target_values_loaded':False,'evaluation_source_values_loaded':False,'development_assays':len(da),'development_targets':len(dt),
      'unique_compounds':len(ids),'task_counts':dict(__import__('collections').Counter(t['split'] for t in dataset)),
      'text_components':dim,'condition_width':bundle['condition_dim'],'common_target_name_width':bundle['target_text_dim'],
      'atom_features':bundle['atom_dim'],'bond_features':14,'assay_parameters_not_inferred':True,
      'target_name_available':sum(bool(names[t]) for t in tids),'missing_structured_fields':['substrate','incubation_time','temperature'],
      'wall_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'dataset_sha256':sha(dest/'dataset.pkl'),'preprocessors_sha256':sha(dest/'preprocessors.pkl')})
    print('Features ready; evaluation potency withheld. '+str(len(dataset))+' tasks',flush=True)

if __name__=='__main__':main()
