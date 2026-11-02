#!/bin/bash
#SBATCH --job-name=polish_2node
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:04:59
#SBATCH --output=validate_multinode.out
#SBATCH -A anakano_429

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/c" || exit 1

echo "Allocated compute nodes:"
scontrol show hostnames "$SLURM_JOB_NODELIST"
echo

gcc -O2 -fopenmp -DTIE_DCHI -o soft soft.c -lm
mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

export OMP_NUM_THREADS=1

./soft field_center.txt 0.03 500 free 1 12 1 0.1 10 - 4 \
    > multinode_openmp.txt

mpirun -bind-to none -n 2 \
    ./soft_mpi field_center.txt 0.03 500 free 1 12 1 0.1 10 - 4 \
    > multinode_mpi.txt

sed -E 's/wall [0-9.]+ s/wall TIME s/' multinode_openmp.txt \
    > multinode_openmp_normalized.txt

sed -E '/^MPI_TIMING /d; s/wall [0-9.]+ s/wall TIME s/' \
    multinode_mpi.txt > multinode_mpi_normalized.txt

if diff -u multinode_openmp_normalized.txt multinode_mpi_normalized.txt
then
    echo
    echo "MULTI-NODE VALIDATION PASSED: results match across two nodes."
else
    echo
    echo "MULTI-NODE VALIDATION FAILED: results differ."
    exit 1
fi
