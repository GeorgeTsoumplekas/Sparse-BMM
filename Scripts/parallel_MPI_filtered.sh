#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10

module load gcc openmpi

make filtered_parallel_MPI

srun -n 20 ./test_parallel_filtered_MPI "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/F_coo.txt" "./matrices/C_filt_coo.txt" 100000
