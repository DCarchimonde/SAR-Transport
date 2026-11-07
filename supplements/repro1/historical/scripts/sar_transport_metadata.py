"""Fetch metadata only for the fixed Final-table ID universe. No activity API calls."""
from __future__ import annotations
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import pathlib
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/sar_transport/metadata'
OUT.mkdir(parents=True, exist_ok=True)
API = 'https://www.ebi.ac.uk/chembl/api/data/'


def fetch(kind, ids=None, fields=None):
    ids = sorted(ids or [])
    params = {'limit': 1000}
    if ids:
        params[f'{kind}_chembl_id__in'] = ','.join(ids)
    if fields:
        params['only'] = ','.join(fields)
    url = API + kind + '.json?' + urllib.parse.urlencode(params)
    key = hashlib.sha256(url.encode()).hexdigest()[:20]
    dest = OUT / f'{kind}_{key}.json'
    log = dest.with_suffix('.request.json')
    if dest.exists() and log.exists():
        return json.loads(dest.read_text())
    for attempt in range(3):
        try:
            start = time.time()
            with urllib.request.urlopen(url, timeout=45) as r:
                data = r.read()
            obj = json.loads(data)
            if obj.get('page_meta', {}).get('next'):
                raise RuntimeError('Unexpected pagination: stop rather than miss metadata')
            dest.write_bytes(data)
            log.write_text(json.dumps({
                'url': url, 'allowed_ids': ids, 'fields': fields,
                'retrieved_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data),
                'wall_seconds': time.time() - start, 'attempt': attempt + 1,
                'activity_records_requested': 0,
            }, indent=2) + '\n')
            return obj
        except Exception:
            if attempt == 2:
                raise
            time.sleep(1 + attempt)


def batches(kind, ids, fields):
    ids = sorted(set(ids))
    chunks = [ids[i:i+50] for i in range(0, len(ids), 50)]
    records = []
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for n, result in enumerate(pool.map(lambda x: fetch(kind, x, fields), chunks), 1):
            records.extend(result[kind + 's'])
            print(f'{kind} metadata {n}/{len(chunks)} batches; no activity request', flush=True)
    returned = {x[f'{kind}_chembl_id'] for x in records}
    if not returned <= set(ids):
        raise RuntimeError('API returned unrequested IDs')
    return records, sorted(set(ids) - returned)


def main():
    start = time.perf_counter()
    universe = json.loads((ROOT / 'data/sar_transport/IDENTITY_UNIVERSE.json').read_text())
    status = fetch('status')
    assays, missing_a = batches('assay', universe['assays'], [
        'assay_chembl_id', 'description', 'document_chembl_id', 'target_chembl_id',
        'assay_type', 'assay_organism', 'assay_tax_id', 'assay_strain',
        'assay_tissue', 'assay_cell_type', 'assay_subcellular_fraction',
        'bao_format', 'confidence_score', 'variant_sequence', 'assay_parameters',
        'src_id', 'src_assay_id',
    ])
    docs = {r['document_chembl_id'] for r in assays if r.get('document_chembl_id')}
    documents, missing_d = batches('document', docs, [
        'document_chembl_id', 'doi', 'pubmed_id', 'journal', 'year',
    ])
    (OUT / 'ASSAYS.json').write_text(json.dumps(sorted(assays, key=lambda x:x['assay_chembl_id']), indent=2) + '\n')
    (OUT / 'DOCUMENTS.json').write_text(json.dumps(sorted(documents, key=lambda x:x['document_chembl_id']), indent=2) + '\n')
    summary = {'api_status': status, 'snapshot_warning': 'Current official API metadata, not ChEMBL32 snapshot',
               'assays_requested': len(universe['assays']), 'assays_returned': len(assays),
               'documents_requested': len(docs), 'documents_returned': len(documents),
               'missing_assays': missing_a, 'missing_documents': missing_d,
               'activity_records_requested': 0, 'wall_seconds': time.perf_counter()-start,
               'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
    (OUT / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
