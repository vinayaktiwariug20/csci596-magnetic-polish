# MPI magnetic-force implementation

## Purpose

This implementation parallelizes the magnetic forces and torques between
nearby flakes. It lives in `c/soft.c` behind the `USE_MPI` compile flag:
without `-DUSE_MPI` the program builds and runs exactly as the serial/OpenMP
version (checked byte for byte), with it the same file builds the MPI version.

The scientific question is whether magnetic interactions between flakes can
cause chaining and reduce their alignment in magnetic nail polish.

## What is parallelized

The expensive `dipole_rates()` function calculates the magnetic force and
torque on each flake.

The MPI version:

1. Divides the target flakes into continuous blocks.
2. Assigns one block to each MPI process.
3. Uses OpenMP to divide each process's local block among threads.
4. Retains the existing cell list, so only nearby flakes inside the cutoff
   are examined.
5. Uses two `MPI_Allgatherv` calls to combine the partial rotational rates
   (`wdip`) and translational rates (`vdip`).

After communication, every process has the complete rates required for the
next simulation step.

## Current MPI design

This is a replicated-state MPI implementation:

- Every process stores all flake positions and orientations.
- Every process builds the same spatial cell list.
- Only the magnetic-force targets are divided among processes.
- Contact handling and time integration are repeated on every process.
- Magnetic rates are synchronized after each force update.

This design preserves the original simulation behavior and provides a clear
MPI baseline. Its limitation is that memory and non-force work are not
distributed. A future spatial-domain version could assign cells to processes
and exchange boundary flakes.

## Build

From the `c` directory:

```bash
mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm
```

The final program compiles with `-Wall -Wextra` without warnings.

## Validation

Four correctness checks passed:

1. `c/validate_mpi.sl`: original OpenMP versus two-process MPI.
2. `c/validate_multinode.sl`: one MPI process on each of two Discovery nodes
   (`a01-02` and `a01-03`).
3. `c/validate_large_mpi.sl`: one versus eight MPI processes with 2,000
   flakes.
4. `c/scale_hybrid.sl`: all MPI/OpenMP layouts produced identical scientific
   output.

The validation messages were:

```text
VALIDATION PASSED: OpenMP and MPI scientific results match.
MULTI-NODE VALIDATION PASSED: results match across two nodes.
LARGE VALIDATION PASSED: 1- and 8-process results match.
HYBRID VALIDATION PASSED: all configurations match.
```

## MPI strong scaling

Run with:

```bash
sbatch c/scale_mpi.sl
```

Test configuration: 2,000 flakes, 3% volume fraction, one-second static hold,
cutoff of four flake radii, 50 force evaluations, and three repeats.

| MPI processes | Median force-kernel time (s) | Force speedup | Force efficiency |
|---:|---:|---:|---:|
| 1 | 7.896 | 1.00 | 100.0% |
| 2 | 3.968 | 1.99 | 99.5% |
| 4 | 2.090 | 3.78 | 94.4% |
| 8 | 1.089 | 7.25 | 90.7% |

The timer includes the force calculation and both `MPI_Allgatherv`
operations. It excludes contact handling and the rest of the simulation.

## Hybrid MPI and OpenMP

Run with:

```bash
sbatch c/scale_hybrid.sl
```

All layouts use eight CPU cores and are repeated three times.

| MPI processes | Threads/process | Force time (s) | Total time (s) |
|---:|---:|---:|---:|
| 1 | 8 | 1.029 | 2.981 |
| 2 | 4 | 1.490 | 4.456 |
| 4 | 2 | 2.041 | 7.352 |
| 8 | 1 | 1.140 | 11.466 |

With eight total cores fixed, increasing the number of MPI processes leaves
fewer OpenMP threads for each copy of the non-force work. The median
non-force time therefore rises from 1.95 to 2.97, 5.31, and 10.33 seconds as
threads per process fall from 8 to 4, 2, and 1. Replicating the work does not
itself add the rank times together because ranks run simultaneously; the
slowdown comes from giving each copy fewer OpenMP threads.

The mixed-layout force timings were higher than in an earlier allocation,
although the force code and CPU model were unchanged. The difference is
therefore treated as shared-node, placement, NUMA, or system-load variability.

## 20,000-flake multi-node scaling

Run with:

```bash
sbatch c/scale_multinode_20k.sl
```

This experiment uses one MPI process and eight OpenMP threads per node. Each
configuration is repeated three times.

| Nodes | Total cores | Force time (s) | Force speedup | Total time (s) | Total speedup |
|---:|---:|---:|---:|---:|---:|
| 1 | 8 | 10.144 | 1.00 | 41.731 | 1.00 |
| 2 | 16 | 7.266 | 1.40 | 36.946 | 1.13 |
| 4 | 32 | 5.023 | 2.02 | 34.689 | 1.20 |

The placement-recorded run used AMD EPYC 7513 nodes and bound every rank to
eight consecutive cores within one socket. The output CSVs from 1, 2, and 4
nodes were byte-identical.

The force calculation scales across nodes, while non-force time remains
roughly 30 to 32 seconds because every node still gives its copy eight
OpenMP threads. On one node, the force kernel is about 24% of total runtime,
so Amdahl's law limits total speedup to about 1.32 times even if the force
calculation were made infinitely fast. The measured four-node total speedup
is 1.20 times: force time falls by about 50%, while total runtime falls by
about 17%.

## Future work

The current MPI design distributes only the magnetic-force calculation; every rank still stores all flakes and repeats contact handling and time integration. A spatial domain decomposition could assign cells and flakes to ranks, exchange only boundary data, and parallelize collisions as well as magnetic forces. This would reduce replicated work and improve full-simulation scaling beyond the current Amdahl limit.

## Results

- `results/mpi_scaling_raw.txt`: raw 2,000-flake MPI timings
- `results/mpi_scaling.csv`: force speedup and efficiency calculations
- `results/mpi_scaling.png`: 2,000-flake force scaling plot
- `results/hybrid_scaling_total_raw.txt`: hybrid force and total timings
- `results/hybrid_runtime.csv`: hybrid median force and total runtimes
- `results/multinode_20k_scaling_raw.txt`: raw multi-node timings
- `results/multinode_20k_scaling.csv`: multi-node scaling calculations
- `results/multinode_20k_scaling.png`: force and total scaling plot
- `c/analyze_scaling.py`: regenerates the original MPI plot
- `c/analyze_runtime_scaling.py`: regenerates total-runtime results
