#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10

module load gcc openmpi

make parallel_MPI

srun -n 20 ./test_parallel_MPI "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt"

