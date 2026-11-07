"""Copy frozen assets to their historical paths. Does not execute scientific code."""
from pathlib import Path
import argparse,csv,hashlib,shutil,zipfile

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--resource-root',required=True,type=Path)
    p.add_argument('--destination',required=True,type=Path)
    a=p.parse_args();here=Path(__file__).resolve().parent
    if a.destination.exists():raise SystemExit('Use a new destination; existing data will not be overwritten.')
    with (here/'RESTORE_MAP.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    for r in rows:
        base=a.resource_root if r['source_distribution']=='v1.0.0' else here
        src=base/r['source_path']
        if not src.is_file() or digest(src)!=r['sha256']:raise SystemExit('Missing or changed source: '+str(src))
    a.destination.mkdir(parents=True)
    for r in rows:
        base=a.resource_root if r['source_distribution']=='v1.0.0' else here
        dst=a.destination/r['original_path'];dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(base/r['source_path'],dst)
    with zipfile.ZipFile(here/'TRAINING_INPUTS.zip') as z:
        for name in z.namelist():
            if Path(name).is_absolute() or '..' in Path(name).parts:raise SystemExit('Unsafe archive path')
        z.extractall(a.destination)
    # Task outputs remain the same bytes, with the original expected directory layout.
    shutil.copytree(a.resource_root/'resource_data',a.destination/'resources/sar_transport_v1/data',dirs_exist_ok=True)
    print('Frozen layout staged. No training, inference, scoring, or raw-data retrieval was executed.')
    print('Read README.md before using any historical script; original stop guards and paths are retained.')

if __name__=='__main__':main()
