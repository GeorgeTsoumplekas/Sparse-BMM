#include "mpi.h"
#include "utilities.h"
#include "bmm_seq.h"


/*----------- Functions for the non-blocked parallel implementation -----------*/

//Transfers matrix A to all processes
comp_matrix* transmit_A(comp_matrix* A, int rank){

    if(rank!=0){
        A = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(A == NULL){
            printf("Couldn't allocate memory for A in transmit_A, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(&A->n,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&A->nnz,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    
    if(rank!=0){
        A->col = (uint32_t*)malloc(A->nnz*sizeof(uint32_t));
        if(A->col == NULL){
            printf("Couldn't allocate memory for A->col in transmit_A, process %d\n",rank);
            exit(-1);
        }

        A->row = (uint32_t*)malloc((A->n+1)*sizeof(uint32_t));
        if(A->row == NULL){
            printf("Couldn't allocate memory for A->row in transmit_A, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(A->col,A->nnz,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(A->row,A->n+1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    return A;
}

//Transfers a chunk of B to each process
comp_matrix* transmit_B_chunks(comp_matrix* B, int rank){

    if(rank!=0){
        B = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(B == NULL){
            printf("Couldn't allocate memory for B in transmit_B_chunks, process %d\n",rank);
            exit(-1);
        }
    }

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    MPI_Request reqs[numtasks-1];
    MPI_Request req_extra[numtasks-1];
    MPI_Status status;

    uint32_t n_total;

    if(rank==0){
        n_total = B->n;
    }

    MPI_Bcast(&n_total,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    int* col_sendcounts = (int*)malloc(numtasks*sizeof(int));
    if(col_sendcounts == NULL){
        printf("Couldn't allocate memory for col_sendcounts in transmit_B_chunks, process %d\n",rank);
        exit(-1);
    }
    for(int i=0;i<numtasks;++i){
        if(i<numtasks-(n_total%numtasks)){
            col_sendcounts[i] = n_total/numtasks + 1;
        }
        else{
            col_sendcounts[i] = n_total/numtasks + 2;
        }
    }

    int* col_displs = (int*)malloc(numtasks*sizeof(int));
    if(col_displs == NULL){
        printf("Couldn't allocate memory for col_displs in transmit_B_chunks, process %d\n",rank);
        exit(-1);
    }

    col_displs[0] = 0;
    for(int i=1;i<numtasks;++i){
        col_displs[i] = col_displs[i-1] + col_sendcounts[i-1] - 1;
    }

    int* row_sendcounts = (int*)malloc(numtasks*sizeof(int));
    if(row_sendcounts == NULL){
        printf("Couldn't allocate memory for row_sendcounts in transmit_B_chunks, process %d\n",rank);
        exit(-1);
    }

    if(rank == 0){
        for(int i=0;i<numtasks-1;++i){
            row_sendcounts[i] = B->col[col_displs[i+1]]-B->col[col_displs[i]];
        }
        row_sendcounts[numtasks-1] = B->nnz - B->col[col_displs[numtasks-1]];
    }

    MPI_Bcast(row_sendcounts,numtasks,MPI_INT,0,MPI_COMM_WORLD);

    int* row_displs = (int*)malloc(numtasks*sizeof(int));
    if(row_displs == NULL){
        printf("Couldn't allocate memory for col_displs in transmit_B_chunks, process %d\n",rank);
        exit(-1);
    }
    row_displs[0] = 0;
    for(int i=1;i<numtasks;++i){
        row_displs[i] = row_displs[i-1] + row_sendcounts[i-1];
    }

    B->nnz = row_sendcounts[rank];
    B->n = col_sendcounts[rank]-1;

    if(rank!=0){
        B->row = (uint32_t*)malloc(B->nnz*sizeof(uint32_t));
        if(B->row == NULL){
            printf("Couldn't allocate memory for B->row in transmit_B_chunks, process %d\n",rank);
            exit(-1);
        }

        B->col = (uint32_t*)malloc((B->n+1)*sizeof(uint32_t));
        if(B->col == NULL){
            printf("Couldn't allocate memory for B->col in transmit_B_chunks, process %d\n",rank);
            exit(-1);
        }
    }

    //Send row array
    if(rank==0){
        MPI_Scatterv(B->row,row_sendcounts,row_displs,MPI_UINT32_T,MPI_IN_PLACE,row_sendcounts[rank],MPI_UINT32_T,0,MPI_COMM_WORLD);
    }
    else{
        MPI_Scatterv(NULL,row_sendcounts,row_displs,MPI_UINT32_T,B->row,row_sendcounts[rank],MPI_UINT32_T,0,MPI_COMM_WORLD);
    }

    //Reallocate row array in root process to the new size
    if(rank==0){
        uint32_t* new_row = (uint32_t*)malloc(B->nnz*sizeof(uint32_t));
        if(new_row == NULL){
            printf("Couldn't allocate memory for new_row in transmit_B_chunks, process %d\n",rank);
            exit(-1);
        }
        memcpy(new_row,B->row,B->nnz*sizeof(uint32_t));
        free(B->row);
        //B->row = NULL; //isws xreastei?
        B->row = new_row;
    }
    
    //Send col array
    if(rank==0){
        for(int i=1;i<numtasks;++i){
            MPI_Isend(&B->col[col_displs[i]+1],col_sendcounts[i]-1,MPI_UINT32_T,i,0,MPI_COMM_WORLD,&reqs[i-1]);
        }

        MPI_Waitall(numtasks-1,reqs,MPI_STATUSES_IGNORE);

        //Reallocate col array in root process to the new size
        uint32_t* new_col = (uint32_t*)malloc((B->n+1)*sizeof(uint32_t));
        if(new_col == NULL){
            printf("Couldn't allocate memory for new_col in transmit_B_chunks, process %d\n",rank);
            exit(-1);
        }
        memcpy(new_col,B->col,(B->n+1)*sizeof(uint32_t));
        free(B->col);
        //B->col = NULL //isws na xreiastei?
        B->col = new_col;
    }
    else{
        MPI_Recv(&B->col[1],col_sendcounts[rank]-1,MPI_UINT32_T,0,0,MPI_COMM_WORLD,&status);
    }

    //Send the extra element missing from the col arrays (except the one in root process)
    if(rank!=numtasks-1){
        int next = rank+1;
        
        MPI_Isend(&B->col[B->n],1,MPI_UINT32_T,next,1,MPI_COMM_WORLD,&req_extra[1]);
    }

    if(rank!=0){
        int prev = rank-1;
        MPI_Irecv(&B->col[0],1,MPI_UINT32_T,prev,1,MPI_COMM_WORLD,&req_extra[0]);
        MPI_Wait(&req_extra[0],MPI_STATUSES_IGNORE);
    }

    free(col_sendcounts);
    free(col_displs);
    free(row_sendcounts);
    free(row_displs);

    return B;
}


comp_matrix* nonblocked_bmm_parallel_2(comp_matrix* A, comp_matrix* B){
    comp_matrix* result = NULL;

    int rank, numtasks;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    MPI_Status stats[2];
    MPI_Request reqs[2];

    

    return result;
}


/* ------------ Functions for the blocked parallel implementation ------------- */


block_comp_matrix* blocked_bmm_parallel_2(block_comp_matrix* A, block_comp_matrix* B){
    
}