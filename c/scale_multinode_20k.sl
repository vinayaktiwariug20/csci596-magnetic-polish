#!/bin/bash
#SBATCH --job-name=polish_20k
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:19:59
#SBATCH --output=multinode_20k_scaling.out
#SBATCH -A anakano_429

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/c" || exit 1

mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=true
export OMP_PLACES=cores

echo "Large multi-node scaling"
echo "N=20000, phi=0.03, static hold=1.0 s"
echo "One MPI process and eight OpenMP threads per node"
echo

for nodes in 1 2 4
do
    for repeat in 1 2 3
    do
        echo "RUN nodes=$nodes processes=$nodes threads=8 repeat=$repeat"

         mpirun -n "$nodes" \
            --map-by ppr:1:node:PE=8 \
            --bind-to core \
            ./soft_mpi field_center.txt 0.03 20000 free 1 \
            12 1 1.0 10 - 4 |
            grep '^MPI_TIMING'
    done
done
