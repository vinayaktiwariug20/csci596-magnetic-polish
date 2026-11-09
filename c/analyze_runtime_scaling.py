#!/usr/bin/env python3

import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

timing_pattern = re.compile(
    r"MPI_TIMING processes=(\d+) calls=(\d+) "
    r"dipole_seconds=([0-9.]+) total_seconds=([0-9.]+)"
)


def read_runs(path, run_pattern):
    records = defaultdict(list)
    current_key = None

    for line in path.read_text(errors="ignore").splitlines():
        run_match = run_pattern.search(line)
        if run_match:
            current_key = tuple(int(value) for value in run_match.groups())
            continue

        timing_match = timing_pattern.search(line)
        if timing_match and current_key is not None:
            records[current_key].append(
                (
                    float(timing_match.group(3)),
                    float(timing_match.group(4)),
                )
            )

    if not records:
        raise SystemExit(f"No timing records found in {path}")

    for key, runs in records.items():
        if len(runs) != 3:
            raise SystemExit(
                f"Expected 3 runs for {key} in {path}, found {len(runs)}"
            )

    return records


hybrid_runs = read_runs(
    RESULTS / "hybrid_scaling_total_raw.txt",
    re.compile(r"RUN processes=(\d+) threads=(\d+) repeat=\d+"),
)

multinode_runs = read_runs(
    RESULTS / "multinode_20k_scaling_raw.txt",
    re.compile(
        r"RUN nodes=(\d+) processes=(\d+) threads=(\d+) repeat=\d+"
    ),
)


with (RESULTS / "hybrid_runtime.csv").open("w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(
        [
            "mpi_processes",
            "openmp_threads",
            "total_cores",
            "median_force_seconds",
            "median_total_seconds",
            "force_fraction",
        ]
    )

    for (processes, threads), runs in sorted(hybrid_runs.items()):
        force = statistics.median(run[0] for run in runs)
        total = statistics.median(run[1] for run in runs)
        writer.writerow(
            [
                processes,
                threads,
                processes * threads,
                f"{force:.6f}",
                f"{total:.6f}",
                f"{force / total:.6f}",
            ]
        )


multi_summary = {}
for (nodes, processes, threads), runs in sorted(multinode_runs.items()):
    force = statistics.median(run[0] for run in runs)
    total = statistics.median(run[1] for run in runs)
    multi_summary[nodes] = {
        "processes": processes,
        "threads": threads,
        "force": force,
        "total": total,
    }

baseline_force = multi_summary[1]["force"]
baseline_total = multi_summary[1]["total"]

with (RESULTS / "multinode_20k_scaling.csv").open(
    "w", newline=""
) as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(
        [
            "nodes",
            "mpi_processes",
            "openmp_threads_per_process",
            "total_cores",
            "median_force_seconds",
            "force_speedup",
            "force_efficiency",
            "median_total_seconds",
            "total_speedup",
            "total_efficiency",
        ]
    )

    for nodes, values in sorted(multi_summary.items()):
        force_speedup = baseline_force / values["force"]
        total_speedup = baseline_total / values["total"]

        writer.writerow(
            [
                nodes,
                values["processes"],
                values["threads"],
                nodes * values["threads"],
                f"{values['force']:.6f}",
                f"{force_speedup:.6f}",
                f"{force_speedup / nodes:.6f}",
                f"{values['total']:.6f}",
                f"{total_speedup:.6f}",
                f"{total_speedup / nodes:.6f}",
            ]
        )


nodes = sorted(multi_summary)
force_times = [multi_summary[n]["force"] for n in nodes]
total_times = [multi_summary[n]["total"] for n in nodes]
force_speedups = [baseline_force / value for value in force_times]
total_speedups = [baseline_total / value for value in total_times]

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

axes[0].plot(nodes, force_times, marker="o", label="Force kernel")
axes[0].plot(nodes, total_times, marker="o", label="Total simulation")
axes[0].set_xlabel("Nodes")
axes[0].set_ylabel("Median time (seconds)")
axes[0].set_title("20,000-flake runtime")
axes[0].set_xticks(nodes)
axes[0].grid(alpha=0.3)
axes[0].legend()

axes[1].plot(nodes, force_speedups, marker="o", label="Force kernel")
axes[1].plot(nodes, total_speedups, marker="o", label="Total simulation")
axes[1].plot(nodes, nodes, linestyle="--", color="gray", label="Ideal")
axes[1].set_xlabel("Nodes")
axes[1].set_ylabel("Speedup")
axes[1].set_title("Multi-node strong scaling")
axes[1].set_xticks(nodes)
axes[1].grid(alpha=0.3)
axes[1].legend()

fig.suptitle(
    "Magnetic flake simulation scaling\n"
    "20,000 flakes, 8 OpenMP threads per MPI process"
)
fig.tight_layout()
fig.savefig(RESULTS / "multinode_20k_scaling.png", dpi=200)

print("Wrote results/hybrid_runtime.csv")
print("Wrote results/multinode_20k_scaling.csv")
print("Wrote results/multinode_20k_scaling.png")

for nodes, values in sorted(multi_summary.items()):
    print(
        f"{nodes} node(s): force={values['force']:.3f}s, "
        f"total={values['total']:.3f}s, "
        f"force speedup={baseline_force / values['force']:.2f}x, "
        f"total speedup={baseline_total / values['total']:.2f}x"
    )
