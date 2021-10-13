#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

module load gcc

make parallel

./test_parallel "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt" 10

