# MPI magnetic-force implementation

## Purpose

This implementation parallelizes the calculation of magnetic forces and
torques between nearby flakes. The original simulation and physical model
remain in `c/soft.c`. The MPI version is in `c/soft_mpi.c`.

The scientific question is whether magnetic interactions between flakes can
cause chaining and reduce their alignment in magnetic nail polish.

## What is parallelized

The expensive `dipole_rates()` function calculates the magnetic force and
torque on every flake.

The MPI version divides the target flakes into continuous blocks. Each MPI
process calculates forces and torques only for the flakes in its block.
OpenMP can also divide that local block among threads.

Afterward, two `MPI_Allgatherv` operations combine the partial rotational
rates (`wdip`) and translational rates (`vdip`). Every process then has the
complete rates needed for the next simulation step.

The existing spatial cell list is retained. It limits calculations to flakes
within the selected cutoff instead of checking every possible pair.

## Current MPI design

This is a replicated-state MPI implementation:

- Every process stores all flake positions and orientations.
- Every process builds the same spatial cell list.
- Only the magnetic-force targets are divided between processes.
- Contact handling and time integration are still repeated on every process.
- `MPI_Allgatherv` synchronizes the magnetic rates after every force update.

This design preserves the original simulation behavior and provides a clear
MPI baseline. Its limitation is that memory and non-force work are not
distributed. A future domain-decomposed version could assign spatial cells
to processes and exchange boundary flakes.

## Build

From the `c` directory:

```bash
mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI \
    -o soft_mpi soft_mpi.c -lm
