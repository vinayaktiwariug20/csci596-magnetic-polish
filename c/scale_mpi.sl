#!/bin/bash
#SBATCH --job-name=polish_scale
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:09:59
#SBATCH --output=mpi_scaling.out
#SBATCH -A anakano_429

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/c" || exit 1

mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

export OMP_NUM_THREADS=1

echo "MPI magnetic-force scaling"
echo "N=2000, phi=0.03, static hold=1.0 s, cutoff=4a"
echo

for processes in 1 2 4 8
do
    for repeat in 1 2 3
    do
        echo "RUN processes=$processes repeat=$repeat"

        mpirun -bind-to core -n "$processes" \
            ./soft_mpi field_center.txt 0.03 2000 free 1 \
            12 1 1.0 10 - 4 |
            grep '^MPI_TIMING'
    done
done
