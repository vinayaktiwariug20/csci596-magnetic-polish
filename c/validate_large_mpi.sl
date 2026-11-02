#!/bin/bash
#SBATCH --job-name=polish_largecheck
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:04:59
#SBATCH --output=validate_large_mpi.out
#SBATCH -A anakano_429

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/c" || exit 1

mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

export OMP_NUM_THREADS=1

mpirun -bind-to core -n 1 \
    ./soft_mpi field_center.txt 0.03 2000 free 1 \
    12 1 1.0 10 - 4 > large_mpi_1.txt

mpirun -bind-to core -n 8 \
    ./soft_mpi field_center.txt 0.03 2000 free 1 \
    12 1 1.0 10 - 4 > large_mpi_8.txt

sed -E '/^MPI_TIMING /d; s/wall [0-9.]+ s/wall TIME s/' \
    large_mpi_1.txt > large_mpi_1_normalized.txt

sed -E '/^MPI_TIMING /d; s/wall [0-9.]+ s/wall TIME s/' \
    large_mpi_8.txt > large_mpi_8_normalized.txt

if diff -u large_mpi_1_normalized.txt large_mpi_8_normalized.txt
then
    echo "LARGE VALIDATION PASSED: 1- and 8-process results match."
else
    echo "LARGE VALIDATION FAILED: results differ."
    exit 1
fi

echo
grep '^MPI_TIMING' large_mpi_1.txt
grep '^MPI_TIMING' large_mpi_8.txt
