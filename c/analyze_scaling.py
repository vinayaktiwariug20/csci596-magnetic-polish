#!/usr/bin/env python3

import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
INPUT = ROOT / "results" / "mpi_scaling_raw.txt"
CSV_OUTPUT = ROOT / "results" / "mpi_scaling.csv"
PLOT_OUTPUT = ROOT / "results" / "mpi_scaling.png"

pattern = re.compile(
    r"MPI_TIMING processes=(\d+) calls=(\d+) "
    r"dipole_seconds=([0-9.]+)"
)

times = defaultdict(list)
calls = {}

for line in INPUT.read_text().splitlines():
    match = pattern.search(line)
    if not match:
        continue

    processes = int(match.group(1))
    calls[processes] = int(match.group(2))
    times[processes].append(float(match.group(3)))

if not times:
    raise SystemExit(f"No MPI timing results found in {INPUT}")

process_counts = sorted(times)
medians = {
    processes: statistics.median(times[processes])
    for processes in process_counts
}

baseline = medians[1]
speedups = {
    processes: baseline / medians[processes]
    for processes in process_counts
}
efficiencies = {
    processes: speedups[processes] / processes
    for processes in process_counts
}

with CSV_OUTPUT.open("w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(
        [
            "mpi_processes",
            "calls",
            "run_1_seconds",
            "run_2_seconds",
            "run_3_seconds",
            "median_seconds",
            "speedup",
            "parallel_efficiency",
        ]
    )

    for processes in process_counts:
        runs = times[processes]
        writer.writerow(
            [
                processes,
                calls[processes],
                *[f"{value:.6f}" for value in runs],
                f"{medians[processes]:.6f}",
                f"{speedups[processes]:.6f}",
                f"{efficiencies[processes]:.6f}",
            ]
        )

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

axes[0].plot(
    process_counts,
    [medians[p] for p in process_counts],
    marker="o",
    linewidth=2,
)
axes[0].set_xlabel("MPI processes")
axes[0].set_ylabel("Magnetic-force time (seconds)")
axes[0].set_title("Force-kernel runtime")
axes[0].set_xticks(process_counts)
axes[0].grid(alpha=0.3)

axes[1].plot(
    process_counts,
    [speedups[p] for p in process_counts],
    marker="o",
    linewidth=2,
    label="Measured",
)
axes[1].plot(
    process_counts,
    process_counts,
    linestyle="--",
    color="gray",
    label="Ideal",
)
axes[1].set_xlabel("MPI processes")
axes[1].set_ylabel("Speedup")
axes[1].set_title("Strong scaling")
axes[1].set_xticks(process_counts)
axes[1].grid(alpha=0.3)
axes[1].legend()

fig.suptitle(
    "Magnetic flake-force MPI scaling\n"
    "N=2000, 3% loading, 50 force evaluations"
)
fig.tight_layout()
fig.savefig(PLOT_OUTPUT, dpi=200)

print(f"Wrote {CSV_OUTPUT}")
print(f"Wrote {PLOT_OUTPUT}")

for processes in process_counts:
    print(
        f"{processes} process(es): "
        f"median={medians[processes]:.3f} s, "
        f"speedup={speedups[processes]:.2f}x, "
        f"efficiency={100 * efficiencies[processes]:.1f}%"
    )
