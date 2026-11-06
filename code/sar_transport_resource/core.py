"""Minimal assay-transfer task API. Standard library only; no scientific fitting."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Iterable, Mapping

FIXED_TRACK = "fixed_c8_v1"
ACQUISITION_TRACK = "acquisition_sentinel_v1"


def _finite(value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Values must be finite")
    return value


@dataclass(frozen=True, slots=True)
class Compound:
    compound_id: str
    canonical_smiles: str
    source_potency: float

    def __post_init__(self):
        object.__setattr__(self, "source_potency", _finite(self.source_potency))
        if not self.compound_id or not self.canonical_smiles:
            raise ValueError("Compound identity and structure are required")


@dataclass(frozen=True, slots=True)
class Edge:
    compound_a: str
    compound_b: str
    weight: float

    def __post_init__(self):
        object.__setattr__(self, "weight", _finite(self.weight))
        if self.compound_a == self.compound_b or self.weight <= 0:
            raise ValueError("Edges require two distinct compounds and positive weight")


def degree_weighted_edges(pairs: Iterable[tuple[str, str]]) -> tuple[Edge, ...]:
    """Original query-graph rule: normalize 1/degree(a) + 1/degree(b)."""
    pairs = tuple(sorted(tuple(sorted(pair)) for pair in pairs))
    if not pairs or len(set(pairs)) != len(pairs) or any(a == b for a, b in pairs):
        raise ValueError("Provide a nonempty graph without duplicate/self edges")
    degree = Counter(x for pair in pairs for x in pair)
    weights = [1 / degree[a] + 1 / degree[b] for a, b in pairs]
    total = math.fsum(weights)
    return tuple(Edge(a, b, w / total) for (a, b), w in zip(pairs, weights))


@dataclass(frozen=True, slots=True)
class TaskView:
    """Predictor-visible observations only. Never contains query target truth.

    fixed_c8_v1 has exactly eight reference labels. acquisition_sentinel_v1
    retains a fixed query graph and a disjoint eligible measurement pool.
    """
    task_id: str
    track_id: str
    legacy_task_id: str
    group_id: str
    source_assay: str
    target_assay: str
    endpoint: str
    potency_unit: str
    target_chembl_id: str
    compounds: tuple[Compound, ...]
    reference_pool_ids: tuple[str, ...]
    reference_labels: tuple[tuple[str, float], ...]
    query_ids: tuple[str, ...]
    query_edges: tuple[Edge, ...]

    def __post_init__(self):
        for name in ("compounds", "reference_pool_ids", "query_ids", "query_edges"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(self, "reference_labels", tuple((str(k), _finite(v)) for k, v in self.reference_labels))
        if self.track_id not in (FIXED_TRACK, ACQUISITION_TRACK):
            raise ValueError("Unknown benchmark track")
        if self.task_id != self.track_id + "::" + self.legacy_task_id:
            raise ValueError("Task ID must include its distinct track identity")
        if not all((self.group_id, self.source_assay, self.target_assay)):
            raise ValueError("Assay and group identities are required")
        if self.endpoint not in ("IC50", "Ki") or self.potency_unit != {"IC50": "pIC50", "Ki": "pKi"}[self.endpoint]:
            raise ValueError("Endpoint and logarithmic potency unit must explicitly agree")
        if not self.target_chembl_id:
            raise ValueError("The shared target identity is required")
        ids = [c.compound_id for c in self.compounds]
        pool, query = set(self.reference_pool_ids), set(self.query_ids)
        refs = [k for k, _ in self.reference_labels]
        if len(ids) != len(set(ids)) or len(pool) != len(self.reference_pool_ids) or len(query) != len(self.query_ids):
            raise ValueError("Compound, query and pool IDs must be unique")
        if pool & query or pool | query != set(ids):
            raise ValueError("Reference pool and fixed query IDs must partition compounds")
        if len(refs) != len(set(refs)) or not set(refs) <= pool:
            raise ValueError("Reference labels must be unique eligible pool observations")
        if self.track_id == FIXED_TRACK and (len(refs) != 8 or set(refs) != pool):
            raise ValueError("The fixed-C8 track requires exactly the original eight labeled references")
        expected = degree_weighted_edges((e.compound_a, e.compound_b) for e in self.query_edges)
        if set(x for e in expected for x in (e.compound_a, e.compound_b)) != query:
            raise ValueError("Query IDs must equal the vertices of the fixed query graph")
        supplied = {(e.compound_a, e.compound_b): e.weight for e in self.query_edges}
        if any((e.compound_a, e.compound_b) not in supplied or not math.isclose(supplied[(e.compound_a, e.compound_b)], e.weight, rel_tol=1e-10, abs_tol=1e-12) for e in expected):
            raise ValueError("Query edge weights must follow the fixed degree-weighting rule")
        object.__setattr__(self, "query_edges", expected)

    def to_dict(self) -> dict:
        return {
            "task_id": self.task_id, "track_id": self.track_id,
            "legacy_task_id": self.legacy_task_id, "group_id": self.group_id,
            "source_assay": self.source_assay, "target_assay": self.target_assay,
            "endpoint": self.endpoint, "potency_unit": self.potency_unit,
            "target_chembl_id": self.target_chembl_id,
            "compounds": [{"compound_id": c.compound_id, "canonical_smiles": c.canonical_smiles, "source_potency": c.source_potency} for c in self.compounds],
            "reference_pool_ids": list(self.reference_pool_ids),
            "reference_labels": [{"compound_id": k, "target_potency": v} for k, v in self.reference_labels],
            "query_ids": list(self.query_ids),
            "query_edges": [{"compound_a": e.compound_a, "compound_b": e.compound_b, "weight": e.weight} for e in self.query_edges],
        }

    @classmethod
    def from_dict(cls, obj: Mapping) -> "TaskView":
        allowed = {"task_id", "track_id", "legacy_task_id", "group_id", "source_assay", "target_assay", "endpoint", "potency_unit", "target_chembl_id", "compounds", "reference_pool_ids", "reference_labels", "query_ids", "query_edges"}
        if set(obj) != allowed:
            raise ValueError("TaskView fields differ from the allowed observable schema")
        return cls(**{k: v for k, v in obj.items() if k not in ("compounds", "reference_labels", "query_edges")},
                   compounds=tuple(Compound(**c) for c in obj["compounds"]),
                   reference_labels=tuple((x["compound_id"], x["target_potency"]) for x in obj["reference_labels"]),
                   query_edges=tuple(Edge(**e) for e in obj["query_edges"]))

    @property
    def query_signature(self) -> str:
        obj = {"query_ids": list(self.query_ids), "edges": [(e.compound_a, e.compound_b, e.weight) for e in self.query_edges]}
        return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True, slots=True)
class TaskScore:
    task_id: str
    group_id: str
    query_signature: str
    L: float
    MAE: float
    delta_MSE: float


def score_predictions(view: TaskView, predictions: Mapping[str, float], query_truth: Mapping[str, float]) -> TaskScore:
    """Score a fixed query graph; the evaluator, not TaskView, supplies truth."""
    expected = set(view.query_ids)
    if set(predictions) != expected or set(query_truth) != expected:
        raise ValueError("Predictions and truth must match exactly the fixed query IDs")
    p, y = ({k: _finite(v) for k, v in obj.items()} for obj in (predictions, query_truth))
    losses, errors = [], []
    for edge in view.query_edges:
        a, b, w = edge.compound_a, edge.compound_b, edge.weight
        gap = p[a] - p[b]
        chosen = min(a, b) if abs(gap) <= 1e-8 else (a if gap > 0 else b)
        losses.append(w * (max(y[a], y[b]) - y[chosen]))
        errors.append(w * (gap - (y[a] - y[b])) ** 2)
    return TaskScore(view.task_id, view.group_id, view.query_signature,
                     math.fsum(losses), math.fsum(abs(p[k] - y[k]) for k in view.query_ids) / len(view.query_ids),
                     math.fsum(errors))


def aggregate_scores(scores: Iterable[TaskScore]) -> dict:
    """Tasks equal within group; final groups equal. Never pool edge counts."""
    scores = tuple(scores)
    if not scores or len({s.task_id for s in scores}) != len(scores):
        raise ValueError("Supply one score per task in a nonempty method run")
    groups = defaultdict(list)
    for score in scores:
        groups[score.group_id].append(score)
    result = {key: math.fsum(math.fsum(getattr(s, key) for s in group) / len(group) for group in groups.values()) / len(groups)
              for key in ("L", "MAE", "delta_MSE")}
    return {"groups": len(groups), "tasks": len(scores), **result}


class BudgetOracle:
    """Logical measurement API, not a security boundary or blinded server.

    The caller retains truth externally and supplies a provider. Each unique
    allowed pool ID costs one; cached repeats cost zero. Query access is rejected.
    New sessions require a fresh trace path. Resume/recovery is not implemented.
    """
    def __init__(self, view: TaskView, truth_provider: Callable[[str], float], budget: int, trace_path: str | Path):
        if view.track_id != ACQUISITION_TRACK or view.reference_labels:
            raise ValueError("Start acquisition with an empty-label acquisition TaskView")
        if isinstance(budget, bool) or not isinstance(budget, int) or budget < 0 or budget > len(view.reference_pool_ids):
            raise ValueError("Budget must be an integer within the reference-pool size")
        self._base = view
        self._provider = truth_provider
        self._budget = budget
        self._cache: dict[str, float] = {}
        self._trace_path = Path(trace_path)
        self._trace_path.parent.mkdir(parents=True, exist_ok=True)
        with self._trace_path.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps({"event": "start", "task_id": view.task_id, "budget": budget, "query_signature": view.query_signature}) + "\n")

    @property
    def spent(self) -> int:
        return len(self._cache)

    @property
    def remaining(self) -> int:
        return self._budget - self.spent

    def _append(self, event: dict):
        event.update(task_id=self._base.task_id, spent=self.spent, remaining=self.remaining)
        with self._trace_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event, sort_keys=True) + "\n")
            handle.flush()

    def reveal(self, compound_id: str) -> float:
        if compound_id in self._base.query_ids:
            self._append({"event": "rejected", "compound_id": compound_id, "reason": "query_is_not_measurement_pool"})
            raise ValueError("Query/sentinel targets cannot be requested")
        if compound_id not in self._base.reference_pool_ids:
            self._append({"event": "rejected", "compound_id": compound_id, "reason": "unknown_pool_id"})
            raise KeyError(compound_id)
        if compound_id in self._cache:
            self._append({"event": "cached", "compound_id": compound_id, "cost": 0})
            return self._cache[compound_id]
        if self.remaining == 0:
            self._append({"event": "rejected", "compound_id": compound_id, "reason": "budget_exhausted"})
            raise RuntimeError("Measurement budget exhausted")
        value = _finite(self._provider(compound_id))
        self._cache[compound_id] = value
        self._append({"event": "revealed", "compound_id": compound_id, "target_potency": value, "cost": 1})
        return value

    def current_view(self) -> TaskView:
        return replace(self._base, reference_labels=tuple(self._cache.items()))
