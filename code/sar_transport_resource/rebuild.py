"""Reconstruct a supplied task from portable tables, not raw ChEMBL curation.

No reference is reselected, no eligibility rule rerun and no scientific score
computed. Original manifests and tabular exports supply the frozen design.
"""
import csv,gzip
from pathlib import Path
from .core import Compound,Edge,TaskView
from .catalog import list_tasks

def rebuild_exported_task(data_root,endpoint,task_id):
    root=Path(data_root)
    rows=[r for r in list_tasks(root) if r['endpoint']==endpoint and r['task_id']==task_id]
    if len(rows)!=1:raise KeyError((endpoint,task_id))
    meta=rows[0]
    def table(name):
        with gzip.open(root/'tables'/name,'rt',newline='') as f:
            return [r for r in csv.DictReader(f,delimiter='\t') if r['endpoint']==endpoint and r['task_id']==task_id]
    compounds=table('OBSERVABLE_COMPOUNDS.tsv.gz');refs=table('REFERENCE_LABELS.tsv.gz');edges=table('QUERY_EDGES.tsv.gz')
    view=TaskView(task_id=task_id,track_id=meta['track_id'],legacy_task_id=meta['legacy_task_id'],
        group_id=meta['group_id'],source_assay=meta['source_assay'],target_assay=meta['target_assay'],
        endpoint=endpoint,potency_unit=meta['potency_unit'],target_chembl_id=meta['target_id'],
        compounds=tuple(Compound(r['compound_id'],r['canonical_smiles'],float(r['source_potency'])) for r in compounds),
        reference_pool_ids=(tuple(r['compound_id'] for r in refs) if meta['track_id']=='fixed_c8_v1' else
                            tuple(r['compound_id'] for r in compounds if r['role']=='reference_pool')),
        reference_labels=tuple((r['compound_id'],float(r['target_potency'])) for r in refs),
        query_ids=tuple(r['compound_id'] for r in compounds if r['role']=='query'),
        query_edges=tuple(Edge(r['compound_a'],r['compound_b'],float(r['weight'])) for r in edges))
    if view.query_signature!=meta['query_signature']:raise ValueError('Rebuilt query does not match frozen manifest')
    return view
