#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

module load gcc openmpi

make filtered_combined

srun -n 2 ./test_filtered_combined "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/F_coo.txt" "./matrices/C_filt_coo.txt" 2 100000
