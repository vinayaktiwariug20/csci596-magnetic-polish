#!/bin/bash
#SBATCH --job-name=polish_hybrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:09:59
#SBATCH --output=hybrid_scaling.out
#SBATCH -A anakano_429

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/c" || exit 1

mpicc -O2 -Wall -Wextra -fopenmp -DTIE_DCHI -DUSE_MPI \
    -o soft_mpi soft.c -lm

echo "Hybrid MPI + OpenMP scaling"
echo "N=2000, phi=0.03, static hold=1.0 s, total cores=8"
echo

reference_created=0

for configuration in "1 8" "2 4" "4 2" "8 1"
do
    read -r processes threads <<< "$configuration"
    export OMP_NUM_THREADS="$threads"
    export OMP_PROC_BIND=true
    export OMP_PLACES=cores

    for repeat in 1 2 3
    do
        output="hybrid_${processes}x${threads}_${repeat}.txt"

        srun --exact \
           --ntasks="$processes" \
           --cpus-per-task="$threads" \
           --cpu-bind=cores \
            ./soft_mpi field_center.txt 0.03 2000 free 1 \
            12 1 1.0 10 - 4 > "$output"

        echo "RUN processes=$processes threads=$threads repeat=$repeat"
        grep '^MPI_TIMING' "$output"

        if [ "$repeat" -eq 1 ]
        then
            normalized="${output%.txt}_normalized.txt"

            sed -E '/^MPI_TIMING /d; s/wall [0-9.]+ s/wall TIME s/' \
                "$output" > "$normalized"

            if [ "$reference_created" -eq 0 ]
            then
                cp "$normalized" hybrid_reference.txt
                reference_created=1
            else
                diff -u hybrid_reference.txt "$normalized"
            fi
        fi
    done
done

echo
echo "HYBRID VALIDATION PASSED: all configurations match."
