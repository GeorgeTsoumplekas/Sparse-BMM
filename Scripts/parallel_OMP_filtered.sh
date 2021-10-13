#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

module load gcc

make filtered_parallel

./test_parallel_filtered "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/F_coo.txt" "./matrices/C_filt_coo.txt" 10 100000