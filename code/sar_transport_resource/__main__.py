import argparse
from dataclasses import asdict
import json
from pathlib import Path
from .core import ACQUISITION_TRACK, BudgetOracle, Compound, TaskView, degree_weighted_edges, score_predictions


def main():
    parser = argparse.ArgumentParser(description="Score saved predictions; no model fitting")
    sub = parser.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo")
    demo.add_argument("--example", default="examples/fixed_c8_task.json")
    demo.add_argument("--private-truth", default="private_fixtures/task_truth.json")
    demo.add_argument("--predictions", default="examples/saved_predictions.json")
    acquisition = sub.add_parser("acquire-demo", help="Synthetic budget-oracle demo; no scientific result")
    acquisition.add_argument("--trace", default="acquisition_demo_trace.jsonl")
    args = parser.parse_args()
    if args.command == "acquire-demo":
        view = TaskView(ACQUISITION_TRACK + "::synthetic_demo", ACQUISITION_TRACK,
                        "synthetic_demo", "synthetic_group", "synthetic_source", "synthetic_target",
                        "IC50", "pIC50", "CHEMBL_SYNTHETIC_TARGET",
                        tuple(Compound(cid, "C", 1.0) for cid in ("q1", "q2", "p1", "p2")),
                        ("p1", "p2"), (), ("q1", "q2"), degree_weighted_edges([("q1", "q2")]))
        private_truth = {"p1": 4.0, "p2": 5.0}
        oracle = BudgetOracle(view, private_truth.__getitem__, 1, args.trace)
        first = oracle.reveal("p1")
        cached = oracle.reveal("p1")
        rejected = {}
        for cid in ("q1", "unknown", "p2"):
            try: oracle.reveal(cid)
            except (ValueError, KeyError, RuntimeError) as exc: rejected[cid] = type(exc).__name__
        print(json.dumps({"synthetic_only": True, "first": first, "cached": cached,
                          "spent": oracle.spent, "remaining": oracle.remaining,
                          "rejected": rejected, "trace": str(args.trace),
                          "query_signature_unchanged": oracle.current_view().query_signature == view.query_signature}, indent=2))
        return
    view = TaskView.from_dict(json.loads(Path(args.example).read_text(encoding="utf-8")))
    truth = json.loads(Path(args.private_truth).read_text(encoding="utf-8"))
    if truth["task_id"] != view.task_id:
        raise ValueError("Truth belongs to a different task")
    predictions = json.loads(Path(args.predictions).read_text(encoding="utf-8"))
    if predictions["task_id"] != view.task_id:
        raise ValueError("Predictions belong to a different task")
    results = {name: asdict(score_predictions(view, values, truth["query_truth"]))
               for name, values in predictions["methods"].items()}
    print(json.dumps({"task_id": view.task_id, "results": results}, indent=2))


if __name__ == "__main__":
    main()
