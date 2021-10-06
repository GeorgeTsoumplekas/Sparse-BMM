#include "blocked_bmm_parallel.h"
#include "blocked_bmm_parallel_MPI.h"


/*----------- Functions for the non-blocked parallel implementation -----------*/


comp_matrix* combined_nonblocked_bmm_parallel(comp_matrix* A, comp_matrix* B, int rank, int thread_num){
    
    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    A = transmit_A(A,rank);

    B = transmit_B_chunks(B,rank);

    comp_matrix* C_chunk;

    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        C_chunk = nonblocked_bmm_parallel(A,B,thread_num);
    }

    comp_matrix* C = concat_C_chunks(C_chunk,rank,numtasks,A->n);

    free_comp_matrix(A);
    free_comp_matrix(B);

    return C;
}


/* ----------------- Functions for the blocked bmm parallel implementation -------------------- */


block_comp_matrix* combined_blocked_bmm_parallel(block_comp_matrix* A, block_comp_matrix* B, int rank, int thread_num){
    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

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

    MPI_Barrier(MPI_COMM_WORLD);

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