"""Current-exposure entry point; raw historical task/manifests stay unchanged."""
from dataclasses import dataclass
import csv, hashlib, json
from pathlib import Path
from .core import TaskView
from .acquisition import AcquisitionContract

@dataclass(frozen=True)
class TaskEntry:
    view: TaskView
    exposure_status: str
    split: str
    initial_reference_contract: AcquisitionContract | None

def list_tasks(data_root):
    with (Path(data_root)/'CURRENT_TASK_MANIFEST.tsv').open(newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    if any(not r.get('exposure_status_current','').startswith('exposed:') for r in rows):
        raise ValueError('RC2 default manifest must state the current exposed status')
    return rows

def load_task_entry(data_root, endpoint, task_id):
    root=Path(data_root)
    rows=[r for r in list_tasks(root) if r['endpoint']==endpoint and r['task_id']==task_id]
    if len(rows)!=1: raise KeyError((endpoint,task_id))
    row=rows[0];p=root/row['task_file'];raw=p.read_bytes()
    view=TaskView.from_dict(json.loads(raw))
    contract=None
    if row.get('initial_reference_contract'):
        contract=AcquisitionContract.from_dict(json.loads((root/row['initial_reference_contract']).read_text()))
        if hashlib.sha256(raw).hexdigest()!=contract.task_definition_sha256:
            raise ValueError('Original task bytes do not match the acquisition contract')
        contract.validate(view)
    if row['query_signature']!=view.query_signature: raise ValueError('Manifest query mismatch')
    return TaskEntry(view,row['exposure_status_current'],row['split'],contract)
