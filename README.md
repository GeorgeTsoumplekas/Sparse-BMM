# Sparse Boolean Matrix Multiplication with two levels of parallelization
# **Parallel and Distributed Computer Systems**  
## **Exercise 4**
### Tsoumplekas Georgios, gktsoump@ece.auth.gr | Arailopoulos Vasileios, varailop@ece.auth.gr

# Introductory

This repository contains the files created for the 4th assignment of Parallel and Distributed Computers Systems. The aim of the assignment was to create a computationaly efficient method to apply blocked multiplication between two sparse matrices (Sparse Blocked Matrix Multiplication). In order to achieve better performance 2 different methods of parallelism have been deployed: process-level parallelism using MPI and thread-level parallelism using OpenMP. We have included both blocked and non-blocked versions that use these libraries separately or combined together. Moreover we have included equivalent filtered versions were we the product is masked by a sparse matrix F. 

# Create the matrices

We have already included a file named ```testmatrices``` which contains some pre-made matrices with which you can test the programs. These matrices are of size n=10000. However, these programs have been made so to be applied to even bigger matrices. We have included a Matlab script named ```BMM_script.m``` which you can execute in Matlab to create the needed matrices with whatever size you want. After running the script the following .txt will be created: A_coo.txt, B_coo.txt, C_coo.txt, C_filt_coo.txt and F_coo.txt which contain the matrices you will need in COO format.

# Build

Use command line ```make all```  to build all programs included at once. It is necessary that MPI is installed so that you can use the ```mpicc``` and ```mpirun``` commands.

# Execution

Filtered versions include both blocked and non-blocked implementations. Note that it is necessary to put the paths to the matrices inside quotes ("").
As for the arguments that appear here ```b``` is the size of the block, ```thread_num``` is the number of threads to be deployed for each process and ```numtasks``` is the number of processes to be executed.

## Serial versions
  + Non-blocked: ```./test_serial "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt"```
  + Blocked: ```./test_blocked "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt b```
  + Filtered: ```./test_serial_filtered "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/F_coo.txt" "path/to/C_filt_coo.txt" b```

## Parallel versions with OpenMP
  + Non-blocked: ```./test_parallel "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt" thread_num```
  + Blocked: ```./test_blocked_parallel "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt" thread_num b```
  + Filtered: ```./test_parallel_filtered "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/F_coo.txt" "path/to/C_filt_coo.txt" thread_num b```

## Parallel versions with MPI
  + Non-blocked: ```mpirun -n numtasks ./test_parallel_MPI "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt"```
  + Blocked: ```mpirun -n numtasks ./test_blocked_parallel_MPI "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt" b```
  + Filtered: ```mpirun -n numtasks ./test_parallel_filtered_MPI "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/F_coo.txt" "path/to/C_filt_coo.txt" b```

## Parallel versions with both OpenMP and MPI
  + Non-blocked: ```mpirun -n numtasks ./test_combined_parallel "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt" thread_num```
  + Blocked: ```mpirun -n numtasks ./test_blocked_combined_parallel "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/C_coo.txt" thread_num b```
  + Filtered: ```mpirun -n numtasks test_filtered_combined "path/to/A_coo.txt" "path/to/B_coo.txt" "path/to/F_coo.txt" "path/to/C_filt_coo.txt" thread_num b```

All necessary information to execute the programs can also be found as comments in the ```Makefile```.

# Execution in HPC

All results presented in the report occured from testing our code in AUTH High Performance Computing (HPC) infrastructure. You can do so, too, by using the scripts provided in the ```Scripts``` file. Feel free to change the number of threads or processes deployed or tweak the rest of the parameters in the scripts to see how the code executes when provided different resources. 
