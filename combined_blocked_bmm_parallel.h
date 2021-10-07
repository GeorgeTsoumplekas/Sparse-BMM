#include "blocked_bmm_parallel.h"
#include "blocked_bmm_parallel_MPI.h"


/*----------- Functions for the non-blocked parallel implementation -----------*/

/**
 * Function that implements non-blocked BMM using two levels of parallelism:
 * Each process handles a chunk of B and each thread of a process handles a chunk of A
**/
comp_matrix* combined_nonblocked_bmm_parallel(comp_matrix* A, comp_matrix* B, int rank, int thread_num){
    
    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    //Send matrix A to all processes
    A = transmit_A(A,rank);

    //Send a chunk of B to each process
    B = transmit_B_chunks(B,rank);

    comp_matrix* C_chunk;

    //No reason to calculate the product of A and B if B is zero
    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        C_chunk = nonblocked_bmm_parallel(A,B,thread_num);
    }

    //Concatenate the product chunks of each process to create the final product matrix
    comp_matrix* C = concat_C_chunks(C_chunk,rank,numtasks,A->n);

    free_comp_matrix(A);
    free_comp_matrix(B);

    return C;
}


/* ----------------- Functions for the blocked bmm parallel implementation -------------------- */


/**
 * Function that implements blocked BMM using two levels of parallelism:
 * Each process handles a chunk of B and each thread of a process handles a chunk of A
**/
block_comp_matrix* combined_blocked_bmm_parallel(block_comp_matrix* A, block_comp_matrix* B, int rank, int thread_num){
    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    //Send matrix A to all processes
    A = transmit_blocked_A(A,rank);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank == 0){
        printf("\nTime elapsed for A transmit: %.5f seconds.\n", elapsed);
    }
    
    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    //Send a chunk of B to each process
    B = transmit_blocked_B_chunks(B,rank);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank == 0){
        printf("\nTime elapsed for B transmit: %.5f seconds.\n", elapsed);
    }

    block_comp_matrix* C_chunk;

    //Proved to be important because otherwise multiplication may start without each process 
    //having the correct chunk of B
    MPI_Barrier(MPI_COMM_WORLD);

    //No reason to calculate the product of A and B if B is zero
    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        C_chunk = blocked_bmm_parallel(A,B,thread_num);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        if(rank == 0){
            printf("\nTime elapsed for blocked_bmm_parallel: %.5f seconds.\n", elapsed);
        }
    }

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    //Concatenate the product chunks of each process to create the final product matrix
    block_comp_matrix* C = concat_blocked_C_chunks(C_chunk,rank,numtasks,A->n_b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank == 0){
        printf("\nTime elapsed for concat_C: %.5f seconds.\n", elapsed);
    }

    free_block_comp_matrix(A);
    free_block_comp_matrix(B);

    return C;
}