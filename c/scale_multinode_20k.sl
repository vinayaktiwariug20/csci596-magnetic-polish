#!/bin/bash
#SBATCH --job-name=polish_20k
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=8
#SBATCH --distribution=block:block
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

csv_files=()
for nodes in 1 2 4
do
    for repeat in 1 2 3
    do
        echo "RUN nodes=$nodes processes=$nodes threads=8 repeat=$repeat"
        echo "PLACEMENT:"
        mpirun --tag-output -n "$nodes" \
            --map-by ppr:1:node:PE=8 \
            --bind-to core \
            bash -c 'echo "HOST $(hostname)"; lscpu | grep "Model name"'
        output_csv="-"
        if [ "$repeat" -eq 1 ]
        then
            output_csv="multinode_20k_${nodes}_flakes.csv"
            csv_files+=("$output_csv")
        fi
         mpirun -n "$nodes" \
            --map-by ppr:1:node:PE=8 \
            --bind-to core \
            --report-bindings \
            ./soft_mpi field_center.txt 0.03 20000 free 1 \
            12 1 1.0 10 "$output_csv" 4 |
            grep '^MPI_TIMING'
    done
done

echo

if cmp -s "${csv_files[0]}" "${csv_files[1]}" &&
   cmp -s "${csv_files[0]}" "${csv_files[2]}"
then
    echo "LARGE CSV VALIDATION PASSED: 1-, 2-, and 4-node outputs are byte-identical."
else
    echo "LARGE CSV VALIDATION FAILED: output CSVs differ."
    exit 1
fi

rm -f "${csv_files[@]}"
