"""Summarize mg_compare histories and plot the roundoff plateau.

Usage: python3 analyze_stagnation.py OUTPUT_DIRECTORY
Requires matplotlib for the plot. Inputs are {macro,split}_diagnostic_history.csv
and {macro,split}_corrected.csv, produced by the report's reproduction commands.
Writes selected CSVs and a PNG in OUTPUT_DIRECTORY.
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


directory = Path(sys.argv[1])
selected, summary, corrected = [], [], []
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
for ax, kind in zip(axes, ("macro", "split")):
    history = read(directory / f"{kind}_diagnostic_history.csv")
    for level in sorted({int(row["refinement"]) for row in history}):
        rows = [r for r in history if int(r["refinement"]) == level]

        def crossing(column, tolerance):
            return next((r["iteration"] for r in rows
                         if float(r[column]) <= tolerance), "")

        summary.append(dict(
            refinement=level, smoother=kind,
            true_1e6=crossing("true_residual", 1e-6),
            preconditioned_1e6=crossing("preconditioned_residual", 1e-6),
            preconditioned_1e8=crossing("preconditioned_residual", 1e-8),
            peak_euclidean=max(float(r["true_residual"]) for r in rows),
        ))
    rows = [r for r in history if r["refinement"] == "6"]
    selected.extend(rows)
    for column, label in (("recursive_residual", "Recursive"),
                          ("true_residual", "True (double matvec)"),
                          ("extended_residual", "Extended accumulation")):
        ax.semilogy([int(r["iteration"]) for r in rows],
                    [float(r[column]) for r in rows], label=label)
    ax.axhline(1e-8, color="black", linestyle=":", label="Original target")
    ax.set(title=kind if kind == "macro" else "Split patches + HCT Jacobi",
           xlabel="PCG iteration", ylabel="Relative Euclidean residual",
           ylim=(1e-16, 1e4))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    corrected.extend(r for r in read(directory / f"{kind}_corrected.csv")
                     if r["refinement"] == "6")
fig.suptitle("Johnson–Mercier, 37,248 unknowns: the recursive residual loses accuracy")
fig.savefig(directory / "mg_stagnation.png", dpi=180)
deeper = directory / "split_corrected_deeper.csv"
if deeper.exists():
    corrected.extend(r for r in read(deeper) if r["refinement"] == "7")
write(directory / "mg_stagnation_history.csv", selected)
write(directory / "mg_stagnation_norms.csv", summary)
write(directory / "mg_stagnation_corrections.csv", corrected)
