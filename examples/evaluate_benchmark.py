"""Evaluate a frozen independent dataset; no automatic tuning or invented labels.

Freeze before scores exist:
python examples/evaluate_benchmark.py freeze dataset.json parameters.json protocol.json

Later, evaluate completed-job scores:
python examples/evaluate_benchmark.py evaluate dataset.json protocol.json scores.json report.json
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biocentric.benchmark import evaluate, freeze_protocol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    freeze = commands.add_parser("freeze")
    freeze.add_argument("dataset", type=Path)
    freeze.add_argument("parameters", type=Path)
    freeze.add_argument("protocol", type=Path)
    run = commands.add_parser("evaluate")
    run.add_argument("dataset", type=Path)
    run.add_argument("protocol", type=Path)
    run.add_argument("scores", type=Path)
    run.add_argument("report", type=Path)
    args = parser.parse_args()
    records = json.loads(args.dataset.read_text())
    if args.command == "freeze":
        freeze_protocol(records, json.loads(args.parameters.read_text()), args.protocol)
        print(f"Protocol frozen: {args.protocol}")
    else:
        report = evaluate(records, json.loads(args.scores.read_text()), json.loads(args.protocol.read_text()))
        with args.report.open("x") as handle:
            json.dump(report, handle, indent=2)
        print(f"Evaluation {report['status']}: {args.report}")


if __name__ == "__main__":
    main()
