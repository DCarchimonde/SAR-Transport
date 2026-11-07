"""Deterministic SAR transport eligibility, identity groups and fixed C8.

This program reads potency ONLY to check finite exact values and conflicting keys.
It never calculates effects, model predictions, rankings or performance metrics.
Numeric observations are written separately under data/sar_transport/restricted.
"""
from __future__ import annotations
import collections as co
import csv
import datetime as dt
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import io
import json
import math
import pathlib
import pickle
import platform
import resource
import sys
import time
import urllib.parse

import numpy as np
from rdkit import Chem, DataStructs, rdBase, RDLogger
from rdkit.Chem import rdFingerprintGenerator

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/sar_transport'
OUT = ROOT / 'reports/sar_transport_v1'
VAULT = DATA / 'restricted'
SOURCE = DATA / 'source/IC50_MinimalCuration'


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def dump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def tsv(p, cols, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if p.suffix == '.gz' else open
    with opener(p, 'wt', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v, ensure_ascii=False, sort_keys=True) if isinstance(v, (list,dict)) else v for k,v in r.items()})


class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        allowed = {('collections','defaultdict'):co.defaultdict,
                   ('builtins','set'):set, ('decimal','Decimal'):Decimal}
        if (module,name) not in allowed:
            raise pickle.UnpicklingError(f'Forbidden cache class {module}.{name}')
        return allowed[module,name]


class UnionFind:
    def __init__(self, nodes):
        self.p = {x:x for x in nodes}
    def find(self, x):
        r=x
        while self.p[r] != r:
            r=self.p[r]
        while self.p[x] != x:
            y=self.p[x]; self.p[x]=r; x=y
        return r
    def union(self,a,b):
        a,b=self.find(a),self.find(b)
        if a!=b:self.p[max(a,b)]=min(a,b)
    def groups(self):
        groups=co.defaultdict(list)
        for x in sorted(self.p):groups[self.find(x)].append(x)
        return dict(groups)


def exact_decimal(value):
    try:
        d=Decimal(str(value).strip())
        return d if d.is_finite() else None
    except (InvalidOperation,ValueError):return None


def value_text(value):
    if value == 0:return '0'
    s=format(value,'f')
    return s.rstrip('0').rstrip('.') if '.' in s else s


def normalize_structure(smiles):
    m=Chem.MolFromSmiles(smiles)
    if m is None:return None,'invalid_smiles',[]
    frags=Chem.GetMolFrags(m,asMols=True,sanitizeFrags=True)
    removed=[]
    if len(frags)>1:
        organic=[f for f in frags if any(a.GetAtomicNum()==6 for a in f.GetAtoms())]
        if len(organic)!=1:return None,'ambiguous_multiple_components',[]
        m=organic[0]
        removed=[Chem.MolToSmiles(f,isomericSmiles=True) for f in frags if not any(a.GetAtomicNum()==6 for a in f.GetAtoms())]
    for atom in m.GetAtoms():atom.SetAtomMapNum(0)
    if not m.GetNumHeavyAtoms():return None,'empty_structure',removed
    return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True),'ok',sorted(removed)


def normalize_doi(x):
    if not x:return ''
    s=urllib.parse.unquote(str(x)).strip().lower()
    for prefix in ('https://doi.org/','http://doi.org/','http://dx.doi.org/','https://dx.doi.org/','doi:'):
        if s.startswith(prefix):s=s[len(prefix):].strip()
    return s


def choose_c8(ids, fps):
    ids=sorted(ids)
    chosen=[ids[0]];remaining=set(ids[1:])
    nearest={c:DataStructs.TanimotoSimilarity(fps[c],fps[chosen[0]]) for c in remaining}
    while len(chosen)<8:
        # Max minimum distance = minimum maximum similarity. Lexical ID breaks ties.
        nxt=min(remaining,key=lambda c:(nearest[c],c))
        chosen.append(nxt);remaining.remove(nxt)
        for c in remaining:
            nearest[c]=max(nearest[c],DataStructs.TanimotoSimilarity(fps[c],fps[nxt]))
    return chosen


def main():
    start=time.perf_counter(); cpu=time.process_time(); started=dt.datetime.now(dt.timezone.utc).isoformat()
    if (OUT/'ELIGIBILITY_SUMMARY.json').exists():
        raise RuntimeError('Existing run preserved. Review existing results instead of silent overwrite.')
    OUT.mkdir(parents=True,exist_ok=True);VAULT.mkdir(parents=True,exist_ok=True)
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    pf=ROOT/'protocol/SAR_TRANSPORT_v1_PREFLIGHT.json'
    assert sha(pf)==json.loads((OUT/'PREFLIGHT_FREEZE.json').read_text())['sha256']
    manifest=json.loads((DATA/'source/ACQUISITION.json').read_text())
    for r in manifest:
        assert sha(DATA/'source'/r['path'])==r['sha256'],r['path']
    meta_summary=json.loads((DATA/'metadata/SUMMARY.json').read_text())
    meta={x['assay_chembl_id']:x for x in json.loads((DATA/'metadata/ASSAYS.json').read_text())}
    docs={x['document_chembl_id']:x for x in json.loads((DATA/'metadata/DOCUMENTS.json').read_text())}
    raw_values=co.defaultdict(set); raw_targets=co.defaultdict(set); raw_counts=co.Counter(); bad_raw=set()
    rows=[]; assays=set(); compounds=set(); targets=set()
    with gzip.open(SOURCE/'IC50_MinimalCuration_Final.csv.gz','rt') as f:
        for n,r in enumerate(csv.DictReader(f),2):
            record={'source_row':n,'target':r['Target_Name'], 'reason':''}
            targets.add(r['Target_Name'])
            for arm in (1,2):
                assay=r[f'pair{arm}_assay_ID'];assays.add(assay);record[f'assay{arm}']=assay
                for side in (1,2):
                    c=r[f'pair{arm}_compoundID_{side}'];compounds.add(c);record[f'c{arm}{side}']=c
                    key=(assay,c);value=exact_decimal(r[f'pair{arm}_pchembl_value_{side}'])
                    raw_counts[key]+=1;raw_targets[key].add(r['Target_Name'])
                    if value is None:bad_raw.add(key)
                    else:raw_values[key].add(value)
            if record['c11']!=record['c21'] or record['c12']!=record['c22']:
                record['reason']='source_pair_identity_mismatch'
            if record['assay1']==record['assay2']:record['reason']='same_assay'
            rows.append(record)
    print('Final identity and exact-value consistency loaded; numeric values not displayed',flush=True)
    with gzip.open(SOURCE/'IC50_MinimalCuration_assaysDict.pkl.gz','rb') as f:
        cache=RestrictedUnpickler(f).load()
    structures=co.defaultdict(set);cache_values=co.defaultdict(set);cache_targets=co.defaultdict(set)
    for aid,entries in cache.items():
        if aid not in assays:continue
        for entry in entries:
            if len(entry)!=5:raise RuntimeError('Unexpected source cache tuple schema')
            c,s,v,a,t=entry
            if c not in compounds:continue
            structures[c].add(s)
            if (a,c) in raw_counts:
                dv=exact_decimal(v)
                if dv is None:bad_raw.add((a,c))
                else:cache_values[(a,c)].add(dv)
                cache_targets[(a,c)].add(t)
    del cache
    mappings={};smiles_ids=co.defaultdict(list);normalization_cache={}
    for c in sorted(compounds):
        normalized=[];errors=[];removed=[]
        for s in sorted(structures[c]):
            if s not in normalization_cache:normalization_cache[s]=normalize_structure(s)
            ns,status,salt=normalization_cache[s]
            if status!='ok':errors.append(status)
            else:normalized.append(ns)
            removed.extend(salt)
        if not structures[c]:errors.append('missing_structure')
        if len(set(normalized))>1:errors.append('same_ID_ambiguous_normalized_structures')
        status=';'.join(sorted(set(errors))) if errors else 'eligible'
        smi=normalized[0] if status=='eligible' else ''
        mappings[c]={'raw_compound_id':c,'canonical_compound_id':'','canonical_isomeric_smiles':smi,
                     'status':status,'source_smiles':sorted(structures[c]),'removed_carbon_free_fragments':sorted(set(removed))}
        if smi:smiles_ids[smi].append(c)
    for smi,ids in smiles_ids.items():
        for c in ids:mappings[c]['canonical_compound_id']=min(ids)
    tsv(OUT/'COMPOUND_IDENTITY.tsv',list(next(iter(mappings.values()))),mappings.values())
    print('Structures normalized with fixed salt and stereo rule',flush=True)
    canonical_values=co.defaultdict(set);canonical_targets=co.defaultdict(set);keyraw=co.defaultdict(set);key_reasons=co.defaultdict(set)
    final_values=co.defaultdict(set);cache_mismatches=[]
    for key,ct in raw_counts.items():
        a,c=key;m=mappings[c]
        if m['status']!='eligible':continue
        ck=(a,m['canonical_compound_id']);keyraw[ck].add(c)
        vals=raw_values[key];cv=cache_values[key]
        canonical_values[ck].update(vals|cv);final_values[ck].update(vals)
        canonical_targets[ck].update(raw_targets[key]|cache_targets[key])
        if key in bad_raw:key_reasons[ck].add('nonfinite_or_non_numeric')
        if not cv or not vals<=cv:
            key_reasons[ck].add('Final_value_not_in_source_cache')
            cache_mismatches.append({'assay_id':a,'raw_compound_id':c,'final_unique_values':len(vals),'cache_unique_values':len(cv)})
        if len(raw_targets[key])!=1 or not raw_targets[key]<=cache_targets[key]:
            key_reasons[ck].add('target_mapping_conflict')
    assay_sources={};assay_targets=co.defaultdict(set)
    for (a,c),ts in canonical_targets.items():assay_targets[a].update(ts)
    for a in sorted(assays):
        meta.setdefault(a, {'assay_chembl_id':a,'metadata_status':'not_returned_by_current_API'})
        mm=meta[a];doc=mm.get('document_chembl_id') or '';doi=normalize_doi(docs.get(doc,{}).get('doi'))
        source_key='doi:'+doi if doi else ('doc:'+doc if doc else '')
        assay_sources[a]={'assay_id':a,'source_key':source_key,'normalized_doi':doi,'document_chembl_id':doc,
                          'source_unknown':not bool(source_key),'current_target_id':mm.get('target_chembl_id') or '',
                          'original_target_ids':sorted(assay_targets[a]),'assay_type':mm.get('assay_type') or '',
                          'description_present':bool(mm.get('description'))}
    valid={};keyledger=[];conflicts=[]
    for ck in sorted(canonical_values):
        a,c=ck;vals=canonical_values[ck];ts=canonical_targets[ck]
        if len(vals)!=1:key_reasons[ck].add('different_reported_values')
        if len(ts)!=1:key_reasons[ck].add('target_mapping_conflict')
        if meta[a].get('target_chembl_id') and meta[a]['target_chembl_id'] not in ts:key_reasons[ck].add('current_metadata_target_conflict')
        reasons=sorted(key_reasons[ck]);status='excluded' if reasons else 'eligible'
        row={'assay_id':a,'canonical_compound_id':c,'raw_compound_ids':sorted(keyraw[ck]),'target_ids':sorted(ts),
             'distinct_values_in_Final':len(final_values[ck]),'distinct_values_with_same_key_cache':len(vals),
             'status':status,'reasons':reasons,'expanded_occurrences':sum(raw_counts[(a,x)] for x in keyraw[ck])}
        keyledger.append(row)
        if len(vals)>1:conflicts.append({**row,'reported_values':[value_text(v) for v in sorted(vals)]})
        if not reasons:valid[ck]={'value':next(iter(vals)),'target':next(iter(ts))}
    tsv(OUT/'OBSERVATION_ELIGIBILITY.tsv',list(keyledger[0]),keyledger)
    tsv(VAULT/'CONFLICT_VALUES.tsv.gz',list(conflicts[0]) if conflicts else ['assay_id','canonical_compound_id','reported_values'],conflicts)
    dump(OUT/'CACHE_MAPPING_CHECK.json',{'Final_keys_missing_or_inconsistent_in_cache':cache_mismatches})
    tsv(OUT/'ASSAY_SOURCE_MAPPING.tsv',list(next(iter(assay_sources.values()))),assay_sources.values())
    # Group every eligible observation before filtering tasks by size or query availability.
    active_assays={a for a,c in valid};uf=UnionFind(active_assays);compound_assays=co.defaultdict(list)
    for a,c in valid:compound_assays[c].append(a)
    for als in compound_assays.values():
        for a in als[1:]:uf.union(als[0],a)
    identity_components=uf.groups();source_assays=co.defaultdict(list)
    for a in sorted(active_assays):
        sk=assay_sources[a]['source_key']
        if sk:source_assays[sk].append(a)
    merge_witness=[]
    for sk,als in sorted(source_assays.items()):
        for a in als[1:]:
            before=(uf.find(als[0]),uf.find(a))
            if before[0]!=before[1]:merge_witness.append({'source_key':sk,'assay_a':als[0],'assay_b':a,'before_group_a':before[0],'before_group_b':before[1]})
            uf.union(als[0],a)
    groups=uf.groups();gids=sorted(groups)
    shuffled=np.random.default_rng(20261031).permutation(gids).tolist()
    nd=math.floor(.70*len(gids));nv=math.floor(.15*len(gids))
    split={g:('development' if i<nd else 'validation' if i<nd+nv else 'evaluation') for i,g in enumerate(shuffled)}
    group_order={g:i for i,g in enumerate(shuffled)}
    tsv(OUT/'SOURCE_MERGE_WITNESSES.tsv',['source_key','assay_a','assay_b','before_group_a','before_group_b'],merge_witness)
    # Canonicalize original edges and keep complete source row provenance.
    task_edges=co.defaultdict(dict);task_targets=co.defaultdict(set);excluded_row_counts=co.Counter();rowledger=[]
    task_rawcounts=co.Counter();task_original_compounds=co.defaultdict(set)
    for r in rows:
        s,t=sorted([r['assay1'],r['assay2']]);tid=s+'__'+t;task_rawcounts[tid]+=1
        task_original_compounds[tid].update([r['c11'],r['c12']]);task_targets[tid].add(r['target'])
        reason=r['reason'];a=b=''
        if not reason:
            if any(mappings[r[k]]['status']!='eligible' for k in ('c11','c12','c21','c22')):
                reason='structure_ineligible'
            else:
                a,b=sorted([mappings[r['c11']]['canonical_compound_id'],mappings[r['c12']]['canonical_compound_id']])
                if a==b:reason='canonical_self_pair'
                elif any((assay,c) not in valid for assay in (s,t) for c in (a,b)):reason='observation_key_excluded'
                elif {valid[(assay,c)]['target'] for assay in (s,t) for c in (a,b)}!={r['target']}:reason='target_mapping_conflict'
        if not reason:
            task_edges[tid].setdefault((a,b),[]).append(r['source_row'])
        else:excluded_row_counts[reason]+=1
        rowledger.append({'source_row':r['source_row'],'task_id':tid,'compound_a':a,'compound_b':b,'status':'excluded' if reason else 'eligible','reason':reason})
    tsv(OUT/'SOURCE_ROW_ELIGIBILITY.tsv.gz',list(rowledger[0]),rowledger)
    del rowledger,rows
    mols={min(ids):Chem.MolFromSmiles(smi) for smi,ids in smiles_ids.items()}
    fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
    fps={c:fpgen.GetFingerprint(m) for c,m in mols.items()}
    tasks=[];contexts=[];queries=[];points=[];all_edges=[]
    for tid in sorted(task_rawcounts):
        s,t=tid.split('__');edges=task_edges.get(tid,{})
        members=sorted({c for e in edges for c in e});reasons=[]
        group=uf.find(s) if s in active_assays else ''
        if s in active_assays and t in active_assays and uf.find(s)!=uf.find(t) and edges:
            raise RuntimeError('Cross-group eligible task')
        if len(task_targets[tid])!=1:reasons.append('task_target_conflict')
        if len(members)<20:reasons.append('fewer_than_20_common_participating_compounds')
        c8=choose_c8(members,fps) if not reasons else []
        qe=sorted(e for e in edges if not set(e)&set(c8)) if c8 else []
        if c8 and not qe:reasons.append('no_query_edge_after_fixed_C8')
        ok=not reasons
        rec={'task_id':tid,'source_assay':s,'target_assay':t,'target_ids':sorted(task_targets[tid]),'group_id':group,
             'split':split.get(group,''),'source_unknown':any(assay_sources[x]['source_unknown'] for x in (s,t)),
             'source_expanded_rows':task_rawcounts[tid],'original_compound_count':len(task_original_compounds[tid]),
             'eligible_compounds':len(members),'eligible_analogue_edges':len(edges),'context_count':len(c8),
             'query_edges':len(qe),'query_points':len({c for e in qe for c in e}),
             'status':'eligible' if ok else 'excluded','reasons':reasons}
        tasks.append(rec)
        qset=set(qe)
        rec['eligible_expanded_rows']=sum(len(v) for v in edges.values())
        for e,source_rows in sorted(edges.items()):
            all_edges.append({'task_id':tid,'compound_a':e[0],'compound_b':e[1],'source_rows':source_rows,'is_query':ok and e in qset})
        if not ok:continue
        contexts.extend({'task_id':tid,'group_id':group,'split':split[group],'context_order':i+1,'canonical_compound_id':c} for i,c in enumerate(c8))
        deg=co.Counter(c for e in qe for c in e);total=sum(1/deg[a]+1/deg[b] for a,b in qe)
        for a,b in qe:
            queries.append({'task_id':tid,'group_id':group,'split':split[group],'compound_a':a,'compound_b':b,
                            'degree_a':deg[a],'degree_b':deg[b],'edge_weight':(1/deg[a]+1/deg[b])/total})
        points.extend({'task_id':tid,'group_id':group,'split':split[group],'canonical_compound_id':c,
                       'allowed_inputs':'structure; source observed pChEMBL; source/target metadata; fixed C8 target values only',
                       'exposure_status':'qualification-program value access; query outcomes not scored or inspected by analyst'} for c in sorted(deg))
    for name,rs,cols in [('TASKS.tsv',tasks,None),('C8.tsv',contexts,['task_id','group_id','split','context_order','canonical_compound_id']),
                         ('QUERY_EDGES.tsv.gz',queries,['task_id','group_id','split','compound_a','compound_b','degree_a','degree_b','edge_weight']),
                         ('QUERY_POINTS.tsv',points,['task_id','group_id','split','canonical_compound_id','allowed_inputs','exposure_status']),
                         ('ALL_ELIGIBLE_EDGES.tsv.gz',all_edges,['task_id','compound_a','compound_b','source_rows','is_query'])]:
        tsv(OUT/name,cols or list(rs[0]),rs)
    eligible_tasks=[r for r in tasks if r['status']=='eligible'];groups_by_split={}
    group_rows=[]
    for g,als in sorted(groups.items()):
        gt=[r for r in tasks if r['group_id']==g];ge=[r for r in gt if r['status']=='eligible']
        gc={c for a,c in valid if a in set(als)}
        group_rows.append({'group_id':g,'split':split[g],'permutation_order':group_order[g],'assay_ids':als,
                           'assays':len(als),'canonical_compounds':len(gc),'tasks':len(gt),'eligible_tasks':len(ge),
                           'query_edges':sum(r['query_edges'] for r in ge),'query_points_across_tasks':sum(r['query_points'] for r in ge),
                           'source_keys':sorted({assay_sources[a]['source_key'] for a in als if assay_sources[a]['source_key']}),
                           'source_unknown_assays':sum(assay_sources[a]['source_unknown'] for a in als)})
    tsv(OUT/'GROUPS.tsv',list(group_rows[0]),group_rows)
    observation_rows=[]
    for (a,c),o in sorted(valid.items()):
        g=uf.find(a)
        observation_rows.append({'assay_id':a,'canonical_compound_id':c,'target_id':o['target'],'group_id':g,'split':split[g],'pchembl':value_text(o['value'])})
    tsv(VAULT/'UNIQUE_OBSERVATIONS.tsv.gz',list(observation_rows[0]),observation_rows)
    tsv(OUT/'UNIQUE_OBSERVATION_IDENTITIES.tsv',list(observation_rows[0])[:-1],observation_rows)
    for sp in ('development','validation','evaluation'):
        et=[r for r in eligible_tasks if r['split']==sp];egs={r['group_id'] for r in et}
        groups_by_split[sp]={'all_groups':sum(split[g]==sp for g in groups),'eligible_final_groups':len(egs),
                            'eligible_tasks':len(et),'unique_query_compounds':len({r['canonical_compound_id'] for r in points if r['split']==sp}),
                            'query_point_task_instances':sum(r['query_points'] for r in et),'query_edges':sum(r['query_edges'] for r in et),
                            'targets':len({t for r in et for t in r['target_ids']}),
                            'unique_observations':sum(r['split']==sp for r in observation_rows),
                            'eligible_expanded_rows':sum(r['eligible_expanded_rows'] for r in et)}
    # Mechanical identity/leakage checks, no resampling or outcome-based checks.
    assert len({(r['task_id'],r['compound_a'],r['compound_b']) for r in queries})==len(queries)
    cb=co.defaultdict(set);qb=co.defaultdict(set);wb=co.defaultdict(float)
    for r in contexts:cb[r['task_id']].add(r['canonical_compound_id'])
    for r in queries:qb[r['task_id']].update([r['compound_a'],r['compound_b']]);wb[r['task_id']]+=r['edge_weight']
    assert all(len(cb[t])==8 and not(cb[t]&qb[t]) for t in cb)
    assert all(abs(w-1)<1e-10 for w in wb.values())
    assert all(len({split[uf.find(a)] for a in als})==1 for als in compound_assays.values())
    assert all(len({split[uf.find(a)] for a in als})==1 for als in source_assays.values())
    gate=groups_by_split['development']['eligible_final_groups']>=100 and groups_by_split['evaluation']['eligible_final_groups']>=30
    summary={'status':'ELIGIBILITY_GATE_PASS' if gate else 'HARD_STOP_INSUFFICIENT_FINAL_GROUPS',
             'source_rows':sum(task_rawcounts.values()),'original_assays':len(assays),'original_compounds':len(compounds),'original_targets':len(targets),
             'raw_assay_compound_keys':len(raw_counts),'canonical_keys_before_exclusion':len(canonical_values),
             'eligible_unique_observations':len(valid),'distinct_value_conflict_keys':len(conflicts),
             'excluded_canonical_keys':sum(r['status']=='excluded' for r in keyledger),
             'canonical_key_exclusion_reasons':dict(co.Counter(reason for r in keyledger for reason in r['reasons'])),
             'structure_status_counts':dict(co.Counter(r['status'] for r in mappings.values())),
             'canonical_compound_identities':len(smiles_ids),'aliased_identity_groups':sum(len(ids)>1 for ids in smiles_ids.values()),
             'carbon_free_counterion_removed_ids':sum(bool(r['removed_carbon_free_fragments']) for r in mappings.values()),
             'excluded_row_counts':dict(excluded_row_counts),'eligible_analogue_edges':len(all_edges),
             'assay_compound_components_before_documents':len(identity_components),'final_groups_after_documents':len(groups),
             'document_merge_operations':len(merge_witness),'largest_group_assays':max(len(x) for x in groups.values()),
             'tasks_total':len(tasks),'tasks_eligible':len(eligible_tasks),'task_exclusion_reasons':dict(co.Counter(x for r in tasks for x in r['reasons'])),
             'split_coverage':groups_by_split,'gate':{'development_at_least':100,'evaluation_at_least':30,'passed':gate},
             'api_version':meta_summary['api_status'],'missing_source_assays':sum(r['source_unknown'] for r in assay_sources.values()),
             'description_present_assays':sum(r['description_present'] for r in assay_sources.values()),
             'new_neural_fits':0,'new_RF_fits':0,'low_dimensional_fits':0,'model_inferences':0,'resampling_draws':0,'effect_analyses':0,
             'old_sealed_outcomes_accessed':False,'wall_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu,
             'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'started_utc':started,
             'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
             'software':{'python':sys.version,'rdkit':rdBase.rdkitVersion,'numpy':np.__version__,'platform':platform.platform()}}
    dump(OUT/'ELIGIBILITY_SUMMARY.json',summary)
    dump(OUT/'ELIGIBILITY_CHECKS.json',{'duplicate_query_edges':0,'context_query_overlap':0,'cross_split_assays':0,'cross_split_canonical_compounds':0,
          'cross_split_known_source_keys':0,'edge_weights_sum_to_one':True,'cache_value_mapping_mismatches':len(cache_mismatches),
          'source_sha256_verified':True,'preflight_sha256_verified':True,'potency_effects_computed':False})
    dump(OUT/'ELIGIBILITY_MANIFEST.json',{'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pf, pathlib.Path(__file__),DATA/'metadata/ASSAYS.json',DATA/'metadata/DOCUMENTS.json']},
          'outputs':{str(p.relative_to(ROOT)):sha(p) for p in sorted(OUT.glob('*')) if p.is_file() and p.name not in ['ELIGIBILITY_MANIFEST.json','ELIGIBILITY_RUN.log']},
          'restricted_outputs':{str(p.relative_to(ROOT)):sha(p) for p in sorted(VAULT.glob('*'))},
          'numeric_access':'All Final values and relevant same-key source-cache values read programmatically for exactness/conflicts; not used to select split/C8 or inspect scientific effects.'})
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
