#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

module load gcc openmpi

make combined_parallel

srun -n 2 ./test_combined_parallel "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/C_coo.txt" 2

