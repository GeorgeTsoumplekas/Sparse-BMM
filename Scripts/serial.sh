#!/bin/bash
#SBATCH --time=07:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load gcc

make serial

srun -N1 -n1 ./test_serial "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt"
