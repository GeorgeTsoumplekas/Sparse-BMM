#include "blocked_bmm_parallel.h"
#include "blocked_bmm_parallel_MPI.h"


/*----------- Function for the non-blocked parallel implementation -----------*/


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


/* ----------------- Function for the blocked bmm parallel implementation -------------------- */


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


/*----------- Function for the filtered non-blocked parallel implementation -----------*/


comp_matrix* combined_nonblocked_bmm_parallel_filtered(comp_matrix* A, comp_matrix* B, comp_matrix* F, int rank, int thread_num){

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    A = transmit_A(A,rank);

    uint32_t total_cols;
    if (rank == 0){
        total_cols = B->n;
    }
    MPI_Bcast(&total_cols,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    F = transmit_F(F,rank,total_cols);

    B = transmit_B_chunks(B,rank);

    comp_matrix* C_chunk = NULL;

    if(B == NULL){
        C_chunk = NULL;
        printf("C is null\n");
    }
    else{
        C_chunk = nonblocked_bmm_parallel_filtered(A,B,F,thread_num);
    }

    comp_matrix* C = concat_C_chunks(C_chunk,rank,numtasks,F->n);

    free_comp_matrix(A);
    free_comp_matrix(B);
    free_comp_matrix(F);    

    return C;
}

/* ----------------- Functions for the filtered blocked bmm parallel implementation -------------------- */

block_comp_matrix* combined_blocked_bmm_parallel_filtered(block_comp_matrix* A, block_comp_matrix* B, block_comp_matrix* F, int rank, int thread_num){

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    A = transmit_blocked_A(A,rank);

    uint32_t n_b;
    if (rank == 0){
        n_b = B->n_b;
    }
    MPI_Bcast(&n_b,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    F = transmit_blocked_F(F,rank,n_b);

    B = transmit_blocked_B_chunks(B,rank);

    block_comp_matrix* C_chunk;

    MPI_Barrier(MPI_COMM_WORLD);

    uint32_t offset;
    int remaining = n_b%numtasks;
    if (rank < numtasks-remaining){
        offset = rank * n_b/numtasks;
    }
    else{
        uint32_t previous_block_cols = (numtasks-remaining)*n_b/numtasks;
        offset = previous_block_cols + (rank-numtasks+remaining)*(n_b/numtasks+1);
    }

    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        C_chunk = blocked_bmm_parallel_filtered(A, B, F, thread_num, offset);
    }

    for (int i=0; i<C_chunk->nnz_blocks; i++){
        for(int j=0; j<C_chunk->blocks[i]->nnz; j++){
            C_chunk->blocks[i]->col[j] = C_chunk->blocks[i]->col[j] - rank*C_chunk->blocks[i]->n;
        }
    }

    block_comp_matrix* C = concat_blocked_C_chunks(C_chunk,rank,numtasks,F->n_b);

    free_block_comp_matrix(A);
    free_block_comp_matrix(B);
    free_block_comp_matrix(F);

    return C;
    
}