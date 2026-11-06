"""Formal common-C4 acquisition, separate from the general-purpose oracle.

Contracts are exported from saved designs, never selected by this module.
Signatures detect accidental mixing; this is not an adversarial security boundary.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Mapping
from .core import ACQUISITION_TRACK, BudgetOracle, TaskView

PROTOCOL_ID = 'sar_transport_acquisition_common_c4_total8_v1'

def object_signature(obj: Mapping) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()

@dataclass(frozen=True)
class AcquisitionContract:
    protocol_id: str
    task_id: str
    endpoint: str
    initial_reference_ids: tuple[str, ...]
    total_budget: int
    additional_unique_references: int
    query_signature: str
    task_definition_sha256: str
    source_file: str
    source_sha256: str
    protocol_signature: str

    def __post_init__(self):
        object.__setattr__(self, 'initial_reference_ids', tuple(self.initial_reference_ids))
        if self.protocol_id != PROTOCOL_ID or self.total_budget != 8 or self.additional_unique_references != 4:
            raise ValueError('The formal track is the frozen common-C4 plus four protocol')
        if len(self.initial_reference_ids) != 4 or len(set(self.initial_reference_ids)) != 4:
            raise ValueError('Exactly four distinct saved initial references are required')
        if self.protocol_signature != object_signature(self.to_dict(include_signature=False)):
            raise ValueError('Contract signature mismatch')

    def to_dict(self, include_signature=True):
        result = dict(self.__dict__)
        result['initial_reference_ids'] = list(self.initial_reference_ids)
        if not include_signature: result.pop('protocol_signature')
        return result

    @classmethod
    def from_dict(cls, obj):
        if set(obj) != set(cls.__dataclass_fields__):
            raise ValueError('Unknown or missing acquisition-contract fields')
        return cls(**obj)

    def validate(self, view):
        if (view.track_id != ACQUISITION_TRACK or view.reference_labels or
            view.task_id != self.task_id or view.endpoint != self.endpoint or
            view.query_signature != self.query_signature):
            raise ValueError('Contract does not identify this empty-label acquisition task')
        if not set(self.initial_reference_ids) <= set(view.reference_pool_ids) or len(view.reference_pool_ids) < 8:
            raise ValueError('Initial references or total budget are outside the fixed pool')

def audit_acquisition_trace(view: TaskView, contract: AcquisitionContract, events: list[dict]) -> dict:
    """Classify a trajectory without labels, fitting, prediction or scoring.

    A generic oracle session has no formal header and cannot pass this check,
    even if it happens to spend eight. Rejected requests are allowed if they
    expose no label and incur no fee. Repeats must refer to an already paid ID.
    """
    contract.validate(view)
    errors, seen, paid = [], set(), []
    if not events or events[0].get('event') != 'start': errors.append('missing_start')
    elif any(events[0].get(k) != v for k,v in {
        'task_id':view.task_id, 'budget':8, 'query_signature':view.query_signature}.items()):
        errors.append('start_mismatch')
    header = events[1] if len(events)>1 else {}
    if (header.get('event') != 'formal_protocol' or
        header.get('protocol_signature') != contract.protocol_signature or
        header.get('initial_reference_ids') != list(contract.initial_reference_ids)):
        errors.append('missing_or_mismatched_formal_header')
    phase='initial'
    for event in events[2:]:
        if event.get('task_id') != view.task_id: errors.append('mixed_task')
        kind = event.get('event')
        if kind == 'initial_complete':
            if tuple(paid) != contract.initial_reference_ids or event.get('spent') != 4:
                errors.append('initial_four_not_paid_before_callback')
            phase='additional'
        elif kind == 'revealed':
            cid=event.get('compound_id')
            if cid in seen or cid not in view.reference_pool_ids or event.get('cost') != 1:
                errors.append('invalid_new_measurement')
            seen.add(cid);paid.append(cid)
            if phase=='initial' and tuple(paid) != contract.initial_reference_ids[:len(paid)]:
                errors.append('initial_order_mismatch')
            if event.get('spent') != len(paid) or event.get('remaining') != 8-len(paid):
                errors.append('billing_mismatch')
        elif kind == 'cached':
            if event.get('compound_id') not in seen or event.get('cost') != 0:
                errors.append('invalid_cached_measurement')
        elif kind == 'rejected':
            if event.get('cost',0) != 0 or 'target_potency' in event:
                errors.append('rejection_exposed_or_charged')
        else: errors.append('unexpected_trace_event')
    if phase != 'additional': errors.append('missing_initial_complete')
    if len(paid) != 8 or tuple(paid[:4]) != contract.initial_reference_ids or len(set(paid[4:])) != 4:
        errors.append('not_exactly_common_four_plus_four')
    return {'protocol_compliant':not errors, 'failure_reasons':sorted(set(errors)),
            'protocol_id':contract.protocol_id, 'protocol_signature':contract.protocol_signature,
            'query_signature':view.query_signature, 'initial_reference_ids':list(contract.initial_reference_ids),
            'cost_used':len(paid), 'total_budget':8, 'additional_unique_references':max(0,len(paid)-4),
            'measurement_ids':paid, 'trace_signature':object_signature({'events':events})}

class FormalAcquisitionRunner:
    """Pay saved common C4 before handing control to an external callback.

    callback(current_view, oracle) chooses only the remaining four measurements
    and returns finite query predictions. It never receives query target truth.
    Under-spending is returned as non-compliant, never a paper-track result.
    Each instance may run once; a fresh trace path is required.
    """
    def __init__(self, view:TaskView, contract:AcquisitionContract,
                 truth_provider:Callable[[str],float], trace_path:str|Path):
        contract.validate(view)
        self.view, self.contract = view, contract
        self.trace_path=Path(trace_path)
        self.oracle=BudgetOracle(view, truth_provider, 8, self.trace_path)
        self.oracle._append({'event':'formal_protocol', 'protocol_id':contract.protocol_id,
            'protocol_signature':contract.protocol_signature,
            'initial_reference_ids':list(contract.initial_reference_ids), 'query_signature':view.query_signature})
        for cid in contract.initial_reference_ids: self.oracle.reveal(cid)
        self.oracle._append({'event':'initial_complete', 'query_signature':view.query_signature})
        self._ran=False

    def run(self, callback:Callable[[TaskView,BudgetOracle],Mapping[str,float]]) -> dict:
        if self._ran: raise RuntimeError('A formal session can run only once')
        self._ran=True
        predictions = dict(callback(self.oracle.current_view(),self.oracle))
        if set(predictions) != set(self.view.query_ids):
            raise ValueError('Predictions must match the unchanged query IDs exactly')
        predictions={cid:float(value) for cid,value in predictions.items()}
        if not all(math.isfinite(x) for x in predictions.values()):
            raise ValueError('Predictions must be finite')
        current=self.oracle.current_view()
        if current.query_signature != self.view.query_signature:
            raise ValueError('Query graph changed during acquisition')
        events=[json.loads(line) for line in self.trace_path.read_text().splitlines()]
        report=audit_acquisition_trace(self.view,self.contract,events)
        return {'task_id':self.view.task_id,'endpoint':self.view.endpoint,
                **report,'predictions':predictions,
                'prediction_signature':object_signature(predictions),
                'scope':'information-access contract, not an adversarial security boundary'}
