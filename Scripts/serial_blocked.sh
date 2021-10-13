#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load gcc

make blocked

srun -N1 -n1 ./test_blocked "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt" 150000
