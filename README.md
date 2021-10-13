# Sparse Boolean Matrix Multiplication with two levels of parallelization
# **Parallel and Distributed Computer Systems**  
## **Exercise 4**
### Tsoumplekas Georgios, gktsoump@ece.auth.gr | Arailopoulos Vasileios, varailop@ece.auth.gr

# Introductory

This repository contains the files created for the 4th assignment of Parallel and Distributed Computers Systems. The aim of the assignment was to create a computationaly efficient method to apply blocked multiplication between two sparse matrices (Sparse Blocked Matrix Multiplication). In order to achieve better performance 2 different methods of parallelism have been deployed: process-level parallelism using MPI and thread-level parallelism using OpenMP. We have included both blocked and non-blocked versions that use these libraries separately or combined together. Moreover we have included equivalent filtered versions were we the product is masked by a sparse matrix F. 

# Create the matrices

We have already included a file named ```testmatrices``` which contains some pre-made matrices with which you can test the programs. These matrices are of size n=10^4. However, these programs have been made so to be applied to even bigger matrices. We have included a Matlab script named ```BMM_script.m``` which you can execute in Matlab to create the needed matrices with whatever size you want. After running the script the following .txt will be created: A_coo.txt, B_coo.txt, C_coo.txt, C_filt_coo.txt and F_coo.txt which contain the matrices you will need in COO format.

# Build

Use command line ```make all```  to build all programs included at once.

# Execution

