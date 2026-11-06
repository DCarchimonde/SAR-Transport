from dataclasses import FrozenInstanceError, replace
import json
from pathlib import Path
import tempfile
import unittest

from sar_transport_resource import (ACQUISITION_TRACK, BudgetOracle, Compound, TaskScore,
                                    TaskView, aggregate_scores, degree_weighted_edges,
                                    score_predictions)


def fixture():
    return TaskView(
        task_id=ACQUISITION_TRACK + "::synthetic", track_id=ACQUISITION_TRACK,
        legacy_task_id="synthetic", group_id="g1", source_assay="s", target_assay="t",
        endpoint="IC50", potency_unit="pIC50", target_chembl_id="CHEMBL_SYNTHETIC_TARGET",
        compounds=tuple(Compound(k, "C", 1.0) for k in ("a", "b", "c", "p1", "p2")),
        reference_pool_ids=("p1", "p2"), reference_labels=(), query_ids=("a", "b", "c"),
        query_edges=degree_weighted_edges([("a", "b"), ("b", "c")]),
    )


class CoreTests(unittest.TestCase):
    def test_endpoint_unit_identity_and_roundtrip(self):
        view = fixture()
        self.assertEqual(TaskView.from_dict(view.to_dict()), view)
        with self.assertRaises(ValueError): replace(view, potency_unit="pKi")
        with self.assertRaises(ValueError): replace(view, target_chembl_id="")
        ki = replace(view, endpoint="Ki", potency_unit="pKi")
        self.assertEqual(ki.potency_unit, "pKi")

    def test_budget_cache_and_forbidden_queries(self):
        calls = []
        def provider(key):
            calls.append(key)
            return {"p1": 4.0, "p2": 5.0}[key]
        with tempfile.TemporaryDirectory() as d:
            trace = Path(d) / "trace.jsonl"
            oracle = BudgetOracle(fixture(), provider, 1, trace)
            self.assertEqual(oracle.reveal("p1"), 4)
            self.assertEqual(oracle.reveal("p1"), 4)
            self.assertEqual((oracle.spent, oracle.remaining), (1, 0))
            with self.assertRaises(ValueError): oracle.reveal("a")
            with self.assertRaises(KeyError): oracle.reveal("unknown")
            with self.assertRaises(RuntimeError): oracle.reveal("p2")
            self.assertEqual(calls, ["p1"])
            events = [json.loads(line) for line in trace.read_text().splitlines()]
            self.assertEqual(sum(e.get("cost", 0) for e in events), 1)
            self.assertEqual([e["event"] for e in events], ["start", "revealed", "cached", "rejected", "rejected", "rejected"])
            updated = oracle.current_view()
            self.assertEqual(updated.query_signature, fixture().query_signature)
            self.assertEqual(updated.reference_labels, (("p1", 4.0),))
            self.assertEqual(updated.query_edges, fixture().query_edges)
            self.assertNotIn("query_truth", updated.to_dict())
            with self.assertRaises(FileExistsError): BudgetOracle(fixture(), provider, 1, trace)

    def test_immutable_view_and_disjoint_pool(self):
        view = fixture()
        with self.assertRaises(FrozenInstanceError): view.query_ids = ("a",)
        with self.assertRaises(ValueError): replace(view, reference_pool_ids=("a", "p1", "p2"))
        bad = view.to_dict();bad["query_truth"] = {"a": 9}
        with self.assertRaises(ValueError): TaskView.from_dict(bad)

    def test_tie_tolerance_is_absolute_and_canonical(self):
        view = fixture()
        truth = {"a": 1, "b": 2, "c": 0}
        score = score_predictions(view, {"a": 0, "b": 0, "c": 0}, truth)
        # Equal edge weights: a-b tie chooses a (loss 1); b-c tie chooses b (loss 0).
        self.assertEqual(score.L, 0.5)
        negative_gap = score_predictions(view, {"a": -1, "b": 0, "c": -2}, truth)
        self.assertEqual(negative_gap.L, 0)
        near_tie = score_predictions(view, {"a": 0, "b": 1e-9, "c": -2}, truth)
        self.assertEqual(near_tie.L, 0.5)

    def test_fixed_queries_and_weights_cannot_be_silently_changed(self):
        view = fixture()
        with self.assertRaises(ValueError): score_predictions(view, {"a": 1, "b": 2}, {"a": 1, "b": 2, "c": 0})
        with self.assertRaises(ValueError): score_predictions(view, {"a": 1, "b": 2, "c": 0, "p1": 4}, {"a": 1, "b": 2, "c": 0})
        with self.assertRaises(ValueError): replace(view, query_edges=(replace(view.query_edges[0], weight=.9), view.query_edges[1]))
        with self.assertRaises(ValueError): score_predictions(view, {"a": float("nan"), "b": 2, "c": 0}, {"a": 1, "b": 2, "c": 0})

    def test_group_equal_aggregation_not_task_microaverage(self):
        scores = [TaskScore("t1", "g1", "q1", 0, 0, 0), TaskScore("t2", "g1", "q2", 2, 2, 2), TaskScore("t3", "g2", "q3", 5, 5, 5)]
        result = aggregate_scores(scores)
        self.assertEqual(result["L"], 3)  # ((0+2)/2 + 5)/2
        with self.assertRaises(ValueError): aggregate_scores(scores + [scores[0]])


if __name__ == "__main__":
    unittest.main()
