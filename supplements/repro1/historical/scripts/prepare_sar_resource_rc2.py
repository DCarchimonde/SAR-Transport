"""RC2 packaging preparation. Copies immutable evidence; no fitting or scoring."""
from pathlib import Path
import csv, hashlib, json, shutil, time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'resources/sar_transport_v1'
OUT = ROOT / 'reports/sar_resource_rc2'
PUB = ROOT / 'PUBLIC_RELEASE_CANDIDATE_RC2'
START = time.process_time()

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def main():
    if PUB.exists(): raise SystemExit('RC2 directory already exists; do not overwrite a prepared candidate')
    OUT.mkdir(parents=True, exist_ok=True)
    old_paths = []
    for base in (OLD, ROOT/'reports/sar_resource_v1', ROOT/'manuscripts/sar_resource_v1'):
        old_paths += [p for p in base.rglob('*') if p.is_file() and not any(x in p.parts for x in ('__pycache__', 'build')) and '.egg-info' not in str(p)]
    put(OUT/'ORIGINAL_ASSET_HASHES.json', {str(p.relative_to(ROOT)): digest(p) for p in sorted(old_paths)})
    put(OUT/'RESUMPTION.json', {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'head': 'df4dd25a2f71dd656ba6d760defe61579df6477a',
        'subsequent_changes_from_planning_head': [],
        'scientific_result': 'd5d3227b8a0269e2a896056ba41abf7d502d768f',
        'implementation_freeze': '6a2b2b1c6f3d580eb0c26db34d3608cda038b2b2',
        'input_package_sha256': digest(ROOT.parent/'deliverables/sar_resource_v1/SAR_Transport_Resource_v1_20261105.zip'),
        'revision_kit_sha256': digest(ROOT.parent/'upload/JCIM_SAR_Resource_Revision_Kit_20261105.zip'),
        'scope': 'RC2 interface, coverage, prose and private packaging only; scientific assets immutable'})
    code = PUB/'code'; code.mkdir(parents=True)
    for name in ('sar_transport_resource', 'tests', 'examples', 'private_fixtures', 'schemas'):
        shutil.copytree(OLD/name, code/name, ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(OLD/'pyproject.toml', code/'pyproject.toml')
    shutil.copytree(OLD/'data', PUB/'resource_data')
    shutil.copytree(OLD/'tables', PUB/'resource_data/tables')
    # Historical manifests and original bytes are preserved. A new manifest is the default.
    rows = list(csv.DictReader((OLD/'tables/CURRENT_TASK_MANIFEST.tsv').open(), delimiter='\t'))
    for r in rows:
        r['exposure_status'] = r['exposure_status_current']
        r['initial_reference_contract'] = ('contracts/'+r['endpoint']+'/'+r['legacy_task_id']+'.json'
                                            if r['track_id']=='acquisition_sentinel_v1' else '')
    with (PUB/'resource_data/CURRENT_TASK_MANIFEST.tsv').open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t');w.writeheader();w.writerows(rows)
    shutil.copy2(OLD/'DATA_EXPOSURE_STATE.json', PUB/'resource_data/DATA_EXPOSURE_STATE.json')
    shutil.copytree(ROOT/'manuscripts/sar_resource_v1', ROOT/'manuscripts/sar_resource_rc2',
                    ignore=shutil.ignore_patterns('*.aux','*.log','*.out','*.bbl','*.blg','*.toc','*.pdf','*.png'))
    # Figures are immutable scientific illustrations, copied byte-for-byte.
    shutil.copytree(ROOT/'manuscripts/sar_resource_v1/figures', ROOT/'manuscripts/sar_resource_rc2/figures', dirs_exist_ok=True)
    shutil.copytree(ROOT.parent/'rc2_inputs/internal_reuse_evidence', OUT/'planning_internal_reuse_evidence')
    shutil.copy2(ROOT.parent/'rc2_inputs/ENGLISH_REPLACEMENTS.txt', OUT/'ENGLISH_REPLACEMENTS.txt')
    shutil.copy2(ROOT.parent/'rc2_inputs/VERIFIED_COVERAGE.json', OUT/'PLANNING_VERIFIED_COVERAGE.json')
    put(OUT/'PREPARATION_BUDGET.json', {'cpu_seconds':time.process_time()-START,
        'scientific_fits':0,'scientific_matrix_scores':0,'resampling':0,'old_model_inference':0})
    print(json.dumps({'status':'RC2_COPY_PREPARED','original_files_snapshotted':len(old_paths)}))

if __name__=='__main__':main()
