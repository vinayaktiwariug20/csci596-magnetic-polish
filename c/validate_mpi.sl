#!/bin/bash
#SBATCH --job-name=polish_validate
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:02:59
#SBATCH --output=validate_mpi.out
#SBATCH -A anakano_429

cd "$SLURM_SUBMIT_DIR/c" || exit 1

gcc -O2 -fopenmp -DTIE_DCHI -o soft soft.c -lm
mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

export OMP_NUM_THREADS=1

./soft field_center.txt 0.03 500 free 1 12 1 0.1 10 - 4 \
    > validation_openmp.txt

mpirun -bind-to none -n 2 \
    ./soft_mpi field_center.txt 0.03 500 free 1 12 1 0.1 10 - 4 \
    > validation_mpi.txt

sed -E 's/wall [0-9.]+ s/wall TIME s/' validation_openmp.txt \
    > validation_openmp_normalized.txt

sed -E '/^MPI_TIMING /d; s/wall [0-9.]+ s/wall TIME s/' validation_mpi.txt \
    > validation_mpi_normalized.txt

if diff -u validation_openmp_normalized.txt validation_mpi_normalized.txt
then
    echo "VALIDATION PASSED: OpenMP and MPI scientific results match."
else
    echo "VALIDATION FAILED: results differ."
    exit 1
fi
