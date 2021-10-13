#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

module load gcc

make parallel_blocked

./test_blocked_parallel "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt" 10 100000