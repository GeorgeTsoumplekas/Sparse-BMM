#!/bin/bash
#SBATCH --time=00:02:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load gcc

make filtered_serial

srun -N1 -n1 ./test_serial_filtered "./matrices/A_coo.txt" "./matrices/B_coo.txt" "./matrices/F_coo.txt" "./matrices/C_filt_coo.txt" 50000