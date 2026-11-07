"""Ki metadata/structure precheck. Never reads potency columns into a table.
No model fitting, effect calculation, or numerical-conflict eligibility.
"""
from pathlib import Path
import urllib.request,urllib.parse,gzip,json,collections,hashlib,datetime,time
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
from rdkit import Chem,RDLogger,rdBase

ROOT=Path('/workspace/scratch/d12ba5af8303');OUT=ROOT/'ki_identity_check_20261104/provenance';OUT.mkdir(exist_ok=True)
BASE='https://raw.githubusercontent.com/Jnelen/ChEMBL_MatchedPairsAnalysis/0a77cdd28abf40a1562aadce727d032cd668815c/'
API='https://www.ebi.ac.uk/chembl/api/data/'
start=datetime.datetime.now(datetime.timezone.utc).isoformat()
ki=pd.read_csv(OUT.parent/'Ki_MinimalCuration_Unpruned_IDENTITY_ONLY.tsv.gz',sep='\t',dtype=str)
assays=set(ki.pair1_assay_ID)|set(ki.pair2_assay_ID)
mols=set(ki.pair1_compoundID_1)|set(ki.pair1_compoundID_2)|set(ki.pair2_compoundID_1)|set(ki.pair2_compoundID_2)
assert len(assays)==382 and len(mols)==2636

def fetch_json(url):
 with urllib.request.urlopen(url,timeout=45) as r:return json.load(r)
status=fetch_json(API+'status.json');(OUT/'CHEMBL_STATUS.json').write_text(json.dumps(status,indent=2))
def get_batch(kind,idfield,ids,fields):
 q=urllib.parse.urlencode({idfield+'__in':','.join(ids),'limit':1000,'only':','.join(fields)})
 obj=fetch_json(API+kind+'.json?'+q)
 return obj[kind+'s']

def get_meta(kind,idfield,ids,fields):
 ids=sorted(ids);parts=[ids[i:i+40] for i in range(0,len(ids),40)]
 allrows=[]
 with ThreadPoolExecutor(max_workers=4) as ex:
  futures=[ex.submit(get_batch,kind,idfield,part,fields) for part in parts]
  for f in as_completed(futures):allrows.extend(f.result())
 assert all(x[idfield] in set(ids) for x in allrows)
 return sorted(allrows,key=lambda x:x[idfield])

meta=get_meta('assay','assay_chembl_id',assays,['assay_chembl_id','document_chembl_id','target_chembl_id','assay_type','description','confidence_score','bao_format','assay_organism'])
(OUT/'ASSAYS.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False))
print('Assay metadata obtained',len(meta),flush=True)
docids={x['document_chembl_id'] for x in meta if x.get('document_chembl_id')}
docs=get_meta('document','document_chembl_id',docids,['document_chembl_id','doi','pubmed_id','journal','year','title'])
(OUT/'DOCUMENTS.json').write_text(json.dumps(docs,indent=2,ensure_ascii=False))
print('Document metadata obtained',len(docs),flush=True)
meta={x['assay_chembl_id']:x for x in meta};docs={x['document_chembl_id']:x for x in docs}

def read_structure_file(item):
 # The stream is never stored. Only ID and structure columns materialize.
 with urllib.request.urlopen(BASE+item['path'],timeout=45) as r:
  with gzip.GzipFile(fileobj=r) as g:
   t=pd.read_csv(g,usecols=['assay_chembl_id','compound_chembl_id','canonical_smiles'],dtype=str)
 return t[t.assay_chembl_id.isin(assays)&t.compound_chembl_id.isin(mols)].drop_duplicates()
files=json.loads((OUT.parent/'DATASET_LIST.json').read_text());pieces=[]
with ThreadPoolExecutor(max_workers=6) as ex:
 for f in as_completed([ex.submit(read_structure_file,x) for x in files]):pieces.append(f.result())
st=pd.concat(pieces,ignore_index=True).drop_duplicates()
st.to_csv(OUT/'SOURCE_ID_STRUCTURE_ONLY.tsv.gz',sep='\t',index=False)
print('Source identity/structure projection obtained',len(st),flush=True)

RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
def normalize_structure(smiles):
 m=Chem.MolFromSmiles(smiles)
 if m is None:return None,'invalid_smiles',[]
 frags=Chem.GetMolFrags(m,asMols=True,sanitizeFrags=True);removed=[]
 if len(frags)>1:
  organic=[f for f in frags if any(a.GetAtomicNum()==6 for a in f.GetAtoms())]
  if len(organic)!=1:return None,'ambiguous_multiple_components',[]
  m=organic[0]
  removed=[Chem.MolToSmiles(f,isomericSmiles=True) for f in frags if not any(a.GetAtomicNum()==6 for a in f.GetAtoms())]
 for atom in m.GetAtoms():atom.SetAtomMapNum(0)
 if not m.GetNumHeavyAtoms():return None,'empty_structure',removed
 return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True),'ok',sorted(removed)
structures=collections.defaultdict(set)
for c,s in st[['compound_chembl_id','canonical_smiles']].itertuples(index=False,name=None):structures[c].add(s)
mappings={};smiids=collections.defaultdict(list)
for c in sorted(mols):
 normalized=[];errors=[];removed=[]
 for s in sorted(structures[c]):
  ns,state,salt=normalize_structure(s)
  if state!='ok':errors.append(state)
  else:normalized.append(ns)
  removed.extend(salt)
 if not structures[c]:errors.append('missing_structure')
 if len(set(normalized))>1:errors.append('same_ID_ambiguous_normalized_structures')
 state=';'.join(sorted(set(errors))) if errors else 'eligible'
 smi=normalized[0] if state=='eligible' else ''
 mappings[c]={'raw_compound_id':c,'canonical_compound_id':'','canonical_isomeric_smiles':smi,'status':state,'removed_carbon_free_fragments':json.dumps(sorted(set(removed)))}
 if smi:smiids[smi].append(c)
for smi,ids in smiids.items():
 for c in ids:mappings[c]['canonical_compound_id']=min(ids)
pd.DataFrame(mappings.values()).to_csv(OUT/'COMPOUND_IDENTITY.tsv',sep='\t',index=False)

def normdoi(x):
 if not x:return ''
 s=urllib.parse.unquote(str(x)).strip().lower()
 for prefix in ('https://doi.org/','http://doi.org/','http://dx.doi.org/','https://dx.doi.org/','doi:'):
  if s.startswith(prefix):s=s[len(prefix):].strip()
 return s
provenance={}
for a in sorted(assays):
 mm=meta.get(a,{});d=mm.get('document_chembl_id') or '';doi=normdoi(docs.get(d,{}).get('doi'))
 provenance[a]={'assay_id':a,'document_chembl_id':d,'normalized_doi':doi,'current_target_id':mm.get('target_chembl_id') or '', 'metadata_returned':bool(mm),'source_unknown':not bool(d or doi)}
pd.DataFrame(provenance.values()).to_csv(OUT/'ASSAY_SOURCE_MAPPING.tsv',sep='\t',index=False)

# Exact original compound, canonical structure, assay and BOTH DOI/document exclusion.
icroot=ROOT/'sar_review_20261104/bundle/reports/sar_transport_v1'
icmap=pd.read_csv(icroot/'COMPOUND_IDENTITY.tsv',sep='\t',dtype=str).fillna('')
icass=pd.read_csv(icroot/'ASSAY_SOURCE_MAPPING.tsv',sep='\t',dtype=str).fillna('')
icraw=set(icmap.raw_compound_id)
icsmi=set(icmap.loc[icmap.status.eq('eligible'),'canonical_isomeric_smiles'])-{''}
ica=set(icass.assay_id);icdoi=set(icass.normalized_doi)-{''};icdoc=set(icass.document_chembl_id)-{''}

class UF:
 def __init__(self,nodes):self.p={x:x for x in nodes}
 def find(self,x):
  while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
  return x
 def union(self,a,b):
  a,b=self.find(a),self.find(b)
  if a!=b:self.p[max(a,b)]=min(a,b)
u=UF(assays);compound_assays=collections.defaultdict(set);raw_assays=collections.defaultdict(set)
for row in ki.itertuples(index=False):
 d=row._asdict()
 for side in (1,2):
  a=d[f'pair{side}_assay_ID']
  for k in (1,2):
   c=d[f'pair{side}_compoundID_{k}'];raw_assays[c].add(a)
   if mappings[c]['status']=='eligible':compound_assays[mappings[c]['canonical_compound_id']].add(a)
for aset in compound_assays.values():
 ss=sorted(aset)
 for a in ss[1:]:u.union(ss[0],a)
for field in ['document_chembl_id','normalized_doi']:
 groups=collections.defaultdict(list)
 for a,m in provenance.items():
  if m[field]:groups[m[field]].append(a)
 for aa in groups.values():
  for a in aa[1:]:u.union(aa[0],a)

overlap=collections.defaultdict(set)
for a in assays:
 g=u.find(a);m=provenance[a]
 if a in ica:overlap[g].add('IC50_assay_ID')
 if m['document_chembl_id'] and m['document_chembl_id'] in icdoc:overlap[g].add('IC50_document_ID')
 if m['normalized_doi'] and m['normalized_doi'] in icdoi:overlap[g].add('IC50_DOI')
for c,aa in raw_assays.items():
 reason=[]
 if c in icraw:reason.append('IC50_raw_compound_ID')
 if mappings[c]['canonical_isomeric_smiles'] in icsmi:reason.append('IC50_canonical_structure')
 for a in aa:overlap[u.find(a)].update(reason)
overlap={g:reasons for g,reasons in overlap.items() if reasons}

cache_keys=set(st[['assay_chembl_id','compound_chembl_id']].itertuples(index=False,name=None))
edges=collections.defaultdict(set);task_targets=collections.defaultdict(set);excluded_rows=collections.Counter()
for row in ki.itertuples(index=False):
 d=row._asdict();a,b=d['pair1_assay_ID'],d['pair2_assay_ID'];target=d['Target_Name']
 cs=[d['pair1_compoundID_1'],d['pair1_compoundID_2'],d['pair2_compoundID_1'],d['pair2_compoundID_2']]
 reason=''
 if a==b:reason='same_assay'
 elif cs[:2]!=cs[2:]:reason='pair_identity_mismatch'
 elif any(mappings[c]['status']!='eligible' for c in cs):reason='structure_ineligible'
 elif any((aid,c) not in cache_keys for aid,cc in [(a,cs[:2]),(b,cs[2:])] for c in cc):reason='missing_identity_in_source_cache'
 elif any(provenance[aid]['current_target_id'] and provenance[aid]['current_target_id']!=target for aid in [a,b]):reason='current_target_conflict'
 if reason:excluded_rows[reason]+=1;continue
 c1,c2=[mappings[c]['canonical_compound_id'] for c in cs[:2]]
 if c1==c2:excluded_rows['collapsed_pair_identity']+=1;continue
 pair=tuple(sorted((a,b)));edges[pair].add(tuple(sorted((c1,c2))));task_targets[pair].add(target)

tasks=[];qrows=[];prows=[]
for (a,b),ee in sorted(edges.items()):
 ids=set(x for e in ee for x in e);ordered=sorted(ids,key=lambda c:(hashlib.sha256(('sentinel_v0|'+c).encode()).hexdigest(),c))
 q=set(ordered[:len(ordered)//2]);pool=set(ordered[len(ordered)//2:]);qe=sorted(e for e in ee if set(e)<=q)
 gid=u.find(a);assert gid==u.find(b)
 known=all(not provenance[x]['source_unknown'] for x in [a,b])
 eligible=len(ids)>=20 and len(pool)>=8 and bool(qe) and len(task_targets[(a,b)])==1
 task_id=a+'__'+b
 tasks.append({'task_id':task_id,'source_assay':a,'target_assay':b,'target_ids':json.dumps(sorted(task_targets[(a,b)])),'group_id':gid,'total_participants':len(ids),'sentinel_count':len(q),'pool_count':len(pool),'sentinel_edges':len(qe),'identity_sentinel_eligible':eligible,'touched_IC50_group':gid in overlap,'overlap_reasons':json.dumps(sorted(overlap.get(gid,set()))),'task_provenance_known':known})
 if eligible and gid not in overlap and known:
  qrows.extend({'task_id':task_id,'compound_a':x,'compound_b':y} for x,y in qe)
  prows.extend({'task_id':task_id,'canonical_compound_id':c,'role':'sentinel' if c in q else 'reference_pool'} for c in ordered)
t=pd.DataFrame(tasks);t.to_csv(OUT/'TASKS_IDENTITY_SENTINEL.tsv',sep='\t',index=False)
pd.DataFrame(qrows).to_csv(OUT/'SENTINEL_EDGES.tsv.gz',sep='\t',index=False)
pd.DataFrame(prows).to_csv(OUT/'SENTINEL_POOL_ASSIGNMENT.tsv',sep='\t',index=False)
allgroups=collections.defaultdict(list)
for a in sorted(assays):allgroups[u.find(a)].append(a)
grows=[{'group_id':g,'assay_ids':json.dumps(aa),'assays':len(aa),'removed_for_IC50_overlap':g in overlap,'overlap_reasons':json.dumps(sorted(overlap.get(g,set()))),'unknown_provenance_assays':sum(provenance[a]['source_unknown'] for a in aa)} for g,aa in sorted(allgroups.items())]
pd.DataFrame(grows).to_csv(OUT/'GROUPS.tsv',sep='\t',index=False)
keep=t[t.identity_sentinel_eligible&~t.touched_IC50_group&t.task_provenance_known]
summary={
 'scope':'Identity, canonical structures, official assay/document metadata only; no Ki potency column materialized/analyzed; numerical qualification pending.',
 'started_utc':start,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':'0a77cdd28abf40a1562aadce727d032cd668815c','chembl_status':status,'rdkit_version':rdBase.rdkitVersion,
 'official_assays_requested':len(assays),'official_assays_returned':len(meta),'official_documents_requested':len(docids),'official_documents_returned':len(docs),'metadata_missing_assays':sorted(assays-set(meta)),
 'structure_cache_files':len(files),'structure_projection_records':len(st),'structure_status':dict(collections.Counter(x['status'] for x in mappings.values())),
 'structure_alias_groups':sum(len(v)>1 for v in smiids.values()),'unique_eligible_canonical_structures':len(smiids),
 'final_groups_before_IC50_overlap_exclusion':len(allgroups),'groups_removed_for_any_IC50_overlap':len(overlap),'remaining_groups':len(allgroups)-len(overlap),
 'overlap_group_reason_counts':dict(collections.Counter(reason for rr in overlap.values() for reason in rr)),
 'identity_row_exclusions':dict(excluded_rows),'tasks_before_overlap_exclusion':len(t),
 'sentinel_rule':"Sort participant canonical IDs by SHA256('sentinel_v0|'+ID), tie ID; first floor(n/2) queries, remaining references; retain original edges within queries.",
 'tasks_with_sentinel_eligibility_before_overlap_exclusion':int(t.identity_sentinel_eligible.sum()),
 'provisional_confirmation_tasks':len(keep),'provisional_confirmation_groups':keep.group_id.nunique(),'provisional_confirmation_targets':len(set(x for s in keep.target_ids for x in json.loads(s))),
 'query_edges':len(qrows),'reference_pool_size_min':int(keep.pool_count.min()) if len(keep) else None,
 'source_unknown_assays':sum(x['source_unknown'] for x in provenance.values()),
 'ki_potency_accessed':False,'conflict_finite_exact_value_qualification_completed':False,
 'not_independent_database':True,'not_confirmation_performance':True,
 'caution':'The counts may decrease after frozen numeric conflict/finite qualification. Same ChEMBL source; published Ki aggregate outcomes and pruning history were already known. Metadata release is current ChEMBL37, not original ChEMBL32.'}
(OUT/'KI_PROVENANCE_SENTINEL_PRECHECK.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)
