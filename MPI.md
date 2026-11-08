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

For this problem size, one MPI process with eight OpenMP threads is fastest
for both the force kernel and the complete simulation. Increasing the number
of MPI processes makes total runtime worse when the total core count is
fixed, because collision handling and the rest of the simulation are repeated
on every rank.

## 20,000-flake multi-node scaling

Run with:

```bash
sbatch c/scale_multinode_20k.sl
```

This experiment uses one MPI process and eight OpenMP threads per node. Each
configuration is repeated three times.

| Nodes | Total cores | Force time (s) | Force speedup | Total time (s) | Total speedup |
|---:|---:|---:|---:|---:|---:|
| 1 | 8 | 10.116 | 1.00 | 42.706 | 1.00 |
| 2 | 16 | 5.565 | 1.82 | 37.705 | 1.13 |
| 4 | 32 | 4.149 | 2.44 | 36.092 | 1.18 |

The force calculation scales substantially across nodes. Total runtime
improves more modestly because collision handling is replicated. Four nodes
reduce force time by about 59% and total runtime by about 15.5%.

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
