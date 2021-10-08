SHELL := /bin/bash


# ===============================
# COMMANDS

CFLAGS = -O3
CC = gcc
MPICC = mpicc
RM = rm -f

# ===============================
# TARGETS

EXECUTABLES = serial blocked parallel parallel_blocked \
			  filtered_serial filtered_parallel parallel_MPI \
			  blocked_parallel_MPI filtered_parallel_MPI combined_parallel \
			  blocked_combined_parallel filtered_combined

default:all

all: $(EXECUTABLES)

### The matrices that can be used as inputs in the following implementations
### are produced by the MATLAB script BBM_script.m. That file produces three
### random matrices A, B and F and the produce the result of the bmm and the
### filtered bmm matrices, C and C_filtered. The txt files that are produced
### serve as inputs in the tests below. In those tests we implement the BMM
### and use matrices C and C_filtered to evaluate the results.

# ======== Serial implementations ======== #

# Code that implements bmm of the input matrices
# USAGE: ./test_serial [A txt file] [B txt file] [C txt file]
serial: test_serial.c
	$(CC) $(CFLAGS) -o test_serial test_serial.c -lm

# Code that implements bmm by dividing the input matrices into blocks
# USAGE: ./test_blocked [A txt file] [B txt file] [C txt file] [no. cols/rows per block (b)]
blocked: test_blocked.c
	$(CC) $(CFLAGS) -o test_blocked test_blocked.c -lm

# Code that implements the filtered bmm by multiplying element-wise with a third matrix
# USAGE: ./test_serial_filtered [A txt file] [B txt file] [F txt file] [C txt file] [b]
filtered_serial: test_serial_filtered.c
	$(CC) $(CFLAGS) -o test_serial_filtered test_serial_filtered.c -lm


# ======== Parallelized implementations with OpenMP ======== #

# Code that implements the bmm of the input matrices using OpenMP
# USAGE: ./test_parallel [A txt file] [B txt file] [C txt file] [no. threads]
parallel: test_parallel.c
	$(CC) $(CFLAGS) -o test_parallel test_parallel.c -lm -fopenmp

# Code that implements bmm by dividing the input matrices into blocks using OpenMP
# USAGE: ./test_blocked_parallel [A txt file] [B txt file] [C txt file] [no. threads] [b]
parallel_blocked: test_blocked_parallel.c
	$(CC) $(CFLAGS) -o test_blocked_parallel test_blocked_parallel.c -lm -fopenmp

# Code that implements the filtered bmm by multiplying element-wise with a third matrix using OpenMP
# USAGE: ./test_parallel_filtered [A txt file] [B txt file] [F txt file] [C txt file] [no. threads] [b]
filtered_parallel: test_parallel_filtered.c
	$(CC) $(CFLAGS) -o test_parallel_filtered test_parallel_filtered.c -lm -fopenmp


# ======== Parallelized implementations with MPI ======== #

# Code that implements the bmm of the input matrices using MPI
# USAGE: mpirun -n [no. processes] ./test_parallel_MPI [A txt file] [B txt file] [C txt file]
parallel_MPI: test_parallel_MPI.c
	$(MPICC) $(CFLAGS) -o test_parallel_MPI test_parallel_MPI.c -lm

# Code that implements bmm by dividing the input matrices into blocks using MPI
# USAGE: mpirun -n [no. processes] ./test_blocked_parallel_MPI [A txt file] [B txt file] [C txt file] [b]
blocked_parallel_MPI: test_blocked_parallel_MPI.c
	$(MPICC) $(CFLAGS) -o test_blocked_parallel_MPI test_blocked_parallel_MPI.c -lm

# Code that implements the filtered bmm by multiplying element-wise with a third matrix using MPI
# USAGE: mpirun -n [no. processes] ./test_parallel_filtered_MPI [A txt file] [B txt file] [F txt file] [C txt file] [b]
filtered_parallel_MPI: test_parallel_filtered_MPI.c
	$(MPICC) $(CFLAGS) -o test_parallel_filtered_MPI test_parallel_filtered_MPI.c -lm


# ======== Parallelized implementations combining the two previous methods ======== #

# USAGE: mpirun -n [no. processes] ./test_combined_parallel [A txt file] [B txt file] [C txt file] [no. threads]
combined_parallel: test_combined_parallel.c
	$(MPICC) $(CFLAGS) -o test_combined_parallel test_combined_parallel.c -lm -fopenmp

# USAGE: mpirun -n [no. processes] ./test_blocked_combined_parallel [A txt file] [B txt file] [C txt file] [no. threads] [b]
blocked_combined_parallel: test_blocked_combined_parallel.c
	$(MPICC) $(CFLAGS) -o test_blocked_combined_parallel test_blocked_combined_parallel.c -lm -fopenmp

# USAGE: mpirun -n [no. processes] ./test_filtered_combined [A txt file] [B txt file] [F txt file] [C txt file] [no. threads] [b]
filtered_combined: test_filtered_combined.c
	$(MPICC) $(CFLAGS) -o test_filtered_combined test_filtered_combined.c -lm -fopenmp

.PHONY: clean

# =================================
# CLEAN

clean:
	$(RM) *.o *~ test_serial test_blocked \
				 test_parallel test_blocked_parallel \
				 test_serial_filtered test_parallel_filtered \
				 test_parallel_MPI test_blocked_parallel_MPI \
				 test_parallel_filtered_MPI test_combined_parallel \
				 test_blocked_combined_parallel test_filtered_combined
