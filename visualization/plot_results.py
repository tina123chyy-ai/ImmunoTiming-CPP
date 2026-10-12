#!/usr/bin/env python3
"""Plot trajectories exported by the C++ engine. matplotlib is optional.

Usage: python3 visualization/plot_results.py results/baseline --output results/baseline.png
"""
import argparse
import csv
from pathlib import Path
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("folder", type=Path, help="Folder of C++ CSV output")
parser.add_argument("--output", type=Path, default=None)
args = parser.parse_args()
files = [x for x in sorted(args.folder.glob("*.csv")) if x.name not in {"summary.csv", "single_metrics.csv", "top_k.csv", "sensitivity.csv"}]
if not files:
    raise SystemExit("No trajectory CSV files found")
for file in files:
    with file.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    plt.plot([float(row["time"]) for row in rows],
             [float(row["tumor"]) for row in rows], label=file.stem.replace("_", " "))
plt.xlabel("Model time (arbitrary units)")
plt.ylabel("Normalized tumor burden")
plt.title("Conceptual model: hypothetical inputs, NOT patient data")
plt.legend()
plt.grid(True, alpha=0.25)
plt.tight_layout()
if args.output:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=160)
    print(f"Saved {args.output}")
else:
    plt.show()
