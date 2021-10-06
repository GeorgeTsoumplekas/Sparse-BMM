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
        printf("Couldn't allocate memory for row_displs in transmit_B_chunks, process %d\n",rank);
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

        int offset = B->col[0];
        for(int i=0;i<B->n+1;++i){
            B->col[i] -= offset;
        }

    }

    free(col_sendcounts);
    free(col_displs);
    free(row_sendcounts);
    free(row_displs);

    return B;
}


comp_matrix* concat_C_chunks(comp_matrix* result, int rank, int numtasks, int n_total){

    MPI_Request rec_reqs[2*(numtasks-1)];
    MPI_Request send_reqs[2];

    uint32_t* chunks_nnz = NULL;
    uint32_t* chunks_n = NULL;
    uint32_t* offset = NULL;

    comp_matrix* C = NULL;

    if(rank==0){
        //n of each chunk
        chunks_n = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(chunks_n == NULL){
            printf("Couldn't allocate memory for chunks_n in concat_C_chunks.\n");
            exit(-1);
        }
        for(int i=0;i<numtasks;++i){
            if(i<numtasks-(n_total%numtasks)){
                chunks_n[i] = n_total/numtasks;
            }
            else{
                chunks_n[i] = n_total/numtasks + 1;
            }
        }

        chunks_nnz = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(chunks_nnz == NULL){
            printf("Couldn't allocate memory for chunks_nnz in concat_C_chunks.\n");
            exit(-1);
        }
    }

    if(result == NULL){
        result = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(result == NULL){
            printf("Couldn't allocate memory for result in concat_C_chunks process %d.\n",rank);
            exit(-1);
        }

        result->nnz = 0;
        result->n = n_total;

        result->col = NULL;
        result->row = (uint32_t*)calloc(n_total+1,sizeof(uint32_t));
        if(result->row == NULL){
            printf("Couldn't allocate memory for result->row in concat_C_chunks process %d.\n",rank);
            exit(-1);
        }
    }

    //nnz of each chunk
    MPI_Gather(&result->nnz,1,MPI_UINT32_T,chunks_nnz,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    //Send cols and rows
    if(rank!=0){
        MPI_Isend(result->row,n_total+1,MPI_UINT32_T,0,2*rank,MPI_COMM_WORLD,&send_reqs[0]);
        MPI_Isend(result->col,result->nnz,MPI_UINT32_T,0,2*rank+1,MPI_COMM_WORLD,&send_reqs[1]);

        MPI_Waitall(2,send_reqs,MPI_STATUSES_IGNORE);

        free_comp_matrix(result);
    }

    if(rank==0){

        comp_matrix** chunks = (comp_matrix**)malloc(numtasks*sizeof(comp_matrix*));
        if(chunks == NULL){
            printf("Couldn't allocate memory for chunks in concat_C_chunks.\n");
            exit(-1);
        }

        for(int i=0;i<numtasks;++i){
            chunks[i] = (comp_matrix*)malloc(sizeof(comp_matrix));
            if(chunks[i] == NULL){
                printf("Couldn't allocate memory for chunks[%d] in concat_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->row = (uint32_t*)malloc((n_total+1)*sizeof(uint32_t));
            if(chunks[i]->row == NULL){
                printf("Couldn't allocate memory for chunks[%d]->row in concat_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->col = (uint32_t*)malloc(chunks_nnz[i]*sizeof(uint32_t));
            if(chunks[i]->col == NULL){
                printf("Couldn't allocate memory for chunks[%d]->col in concat_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->nnz = chunks_nnz[i];
            chunks[i]->n = chunks_n[i];

            if(i==0){
                memcpy(chunks[i]->col,result->col,chunks[i]->nnz*sizeof(uint32_t));
                memcpy(chunks[i]->row,result->row,(n_total+1)*sizeof(uint32_t));
            }
            else{
                MPI_Irecv(chunks[i]->row,n_total+1,MPI_UINT32_T,i,2*i,MPI_COMM_WORLD,&rec_reqs[i-1]);
                MPI_Irecv(chunks[i]->col,chunks[i]->nnz,MPI_UINT32_T,i,2*i+1,MPI_COMM_WORLD,&rec_reqs[numtasks-1+i-1]);
            }
        }

        free_comp_matrix(result);
        free(chunks_nnz);
        free(chunks_n);

        MPI_Waitall(2*(numtasks-1),rec_reqs,MPI_STATUSES_IGNORE);

        //Create the final product csr matrix
        
        uint32_t total_nnz = 0;     //Total number of non-zero elements in the whole matrix
        for(int i=0;i<numtasks;++i){
            total_nnz += chunks[i]->nnz;
        }

        if(total_nnz==0){
            return NULL;
        }

        C = new_comp_matrix(total_nnz,n_total,"csr");

        uint32_t nnz_count=0;    //number of non-zero elements found up to this point
        uint32_t row_in_chunk_start, row_in_chunk_end;

        offset = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(offset == NULL){
            printf("Couldn't allocate memory for offset in concat_C_chunks.\n");
            exit(-1);
        }
        offset[0] = 0;
        for(int i=1;i<numtasks;++i){
            offset[i] = offset[i-1]+chunks[i-1]->n;
        }

        C->row[0] = 0;

        for(uint32_t i=0;i<n_total;++i){
            for(uint32_t j=0;j<numtasks;++j){
                row_in_chunk_start = chunks[j]->row[i];
                row_in_chunk_end = chunks[j]->row[i+1];

                for(uint32_t k=row_in_chunk_start;k<row_in_chunk_end;++k){
                    C->col[nnz_count] = chunks[j]->col[k] + offset[j];

                    nnz_count++;
                }

            }

            C->row[i+1] = nnz_count;
        }

        free(offset);

        for(int i=0;i<numtasks;++i){
            free_comp_matrix(chunks[i]);
        }
        free(chunks);
    }

    
    return C;
}


comp_matrix* nonblocked_bmm_parallel_2(comp_matrix* A, comp_matrix* B, int rank){

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    A = transmit_A(A,rank);

    B = transmit_B_chunks(B,rank);

    comp_matrix* C_chunk;

    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        C_chunk = bmm_seq(A,B,0);
    }

    comp_matrix* C = concat_C_chunks(C_chunk,rank,numtasks,A->n);

    free_comp_matrix(A);
    free_comp_matrix(B);

    return C;
}


/* ------------ Functions for the blocked parallel implementation ------------- */


block_comp_matrix* transmit_blocked_A(block_comp_matrix* A, int rank){

    if(rank!=0){
        A = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(A == NULL){
            printf("Couldn't allocate memory for A in transmit_blocked_A, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(&A->n_b,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&A->nnz_blocks,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&A->real_dim,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    if(rank!=0){
        A->block_col = (uint32_t*)malloc(A->nnz_blocks*sizeof(uint32_t));
        if(A->block_col == NULL){
            printf("Couldn't allocate memory for A->block_col in transmit_blocked_A, process %d\n",rank);
            exit(-1);
        }

        A->block_row = (uint32_t*)malloc((A->n_b+1)*sizeof(uint32_t));
        if(A->block_row == NULL){
            printf("Couldn't allocate memory for A->block_row in transmit_blocked_A, process %d\n",rank);
            exit(-1);
        }

        A->blocks = (comp_matrix**)malloc(A->nnz_blocks*sizeof(comp_matrix*));
        if(A->blocks == NULL){
            printf("Couldn't allocate memory for A->blocks in transmit_blocked_A, process %d\n",rank);
            exit(-1);
        }
        for(int i=0;i<A->nnz_blocks;++i){
            A->blocks[i] = NULL;
        }
    }

    MPI_Bcast(A->block_col,A->nnz_blocks,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(A->block_row,A->n_b+1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    for(int i=0;i<A->nnz_blocks;++i){
        A->blocks[i] = transmit_A(A->blocks[i],rank);
    }

    return A;
}


block_comp_matrix* transmit_blocked_B_chunks(block_comp_matrix* B, int rank){

    if(rank!=0){
        B = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(B == NULL){
            printf("Couldn't allocate memory for B in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }
    }

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    MPI_Request reqs[numtasks-1];
    MPI_Request req_extra[numtasks-1];
    MPI_Status status;

    uint32_t n_b_total;
    uint32_t real_dim;
    uint32_t total_nnz_blocks;

    if(rank==0){
        n_b_total = B->n_b;
        real_dim = B->real_dim;
        total_nnz_blocks = B->nnz_blocks;
    }

    MPI_Bcast(&n_b_total,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&real_dim,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    B->real_dim = real_dim;

    int* block_col_sendcounts = (int*)malloc(numtasks*sizeof(int));
    if(block_col_sendcounts == NULL){
        printf("Couldn't allocate memory for block_col_sendcounts in transmit_blocked_B_chunks, process %d\n",rank);
        exit(-1);
    }
    for(int i=0;i<numtasks;++i){
        if(i<numtasks-(n_b_total%numtasks)){
            block_col_sendcounts[i] = n_b_total/numtasks + 1;
        }
        else{
            block_col_sendcounts[i] = n_b_total/numtasks + 2;
        }
    }

    int* block_col_displs = (int*)malloc(numtasks*sizeof(int));
    if(block_col_displs == NULL){
        printf("Couldn't allocate memory for block_col_displs in transmit_blocked_B_chunks, process %d\n",rank);
        exit(-1);
    }

    block_col_displs[0] = 0;
    for(int i=1;i<numtasks;++i){
        block_col_displs[i] = block_col_displs[i-1] + block_col_sendcounts[i-1] - 1;
    }

    int* block_row_sendcounts = (int*)malloc(numtasks*sizeof(int));
    if(block_row_sendcounts == NULL){
        printf("Couldn't allocate memory for block_row_sendcounts in transmit_blocked_B_chunks, process %d\n",rank);
        exit(-1);
    }

    if(rank == 0){
        for(int i=0;i<numtasks-1;++i){
            block_row_sendcounts[i] = B->block_col[block_col_displs[i+1]]-B->block_col[block_col_displs[i]];
        }
        block_row_sendcounts[numtasks-1] = B->nnz_blocks - B->block_col[block_col_displs[numtasks-1]];
    }

    MPI_Bcast(block_row_sendcounts,numtasks,MPI_INT,0,MPI_COMM_WORLD);

    int* block_row_displs = (int*)malloc(numtasks*sizeof(int));
    if(block_row_displs == NULL){
        printf("Couldn't allocate memory for block_row_displs in transmit_blocked_B_chunks, process %d\n",rank);
        exit(-1);
    }
    block_row_displs[0] = 0;
    for(int i=1;i<numtasks;++i){
        block_row_displs[i] = block_row_displs[i-1] + block_row_sendcounts[i-1];
    }

    B->nnz_blocks = block_row_sendcounts[rank];
    B->n_b = block_col_sendcounts[rank]-1;
    
    if(rank!=0){
        B->block_row = (uint32_t*)malloc(B->nnz_blocks*sizeof(uint32_t));
        if(B->block_row == NULL){
            printf("Couldn't allocate memory for B->block_row in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }

        B->block_col = (uint32_t*)malloc((B->n_b+1)*sizeof(uint32_t));
        if(B->block_col == NULL){
            printf("Couldn't allocate memory for B->block_col in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }

        B->blocks = (comp_matrix**)malloc(B->nnz_blocks*sizeof(comp_matrix*));
        if(B->blocks == NULL){
            printf("Couldn't allocate memory for B->blocks in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }
        for(int i=0;i<B->nnz_blocks;++i){
            B->blocks[i] = (comp_matrix*)malloc(sizeof(comp_matrix));
            if(B->blocks[i] == NULL){
                printf("Couldn't allocate memory for B->blocks[%d] in transmit_blocked_B_chunks, process %d.\n",i,rank);
                exit(-1);
            }
        }
    }

    //Send block_row array
    if(rank==0){
        MPI_Scatterv(B->block_row,block_row_sendcounts,block_row_displs,MPI_UINT32_T,MPI_IN_PLACE,block_row_sendcounts[rank],MPI_UINT32_T,0,MPI_COMM_WORLD);
    }
    else{
        MPI_Scatterv(NULL,block_row_sendcounts,block_row_displs,MPI_UINT32_T,B->block_row,block_row_sendcounts[rank],MPI_UINT32_T,0,MPI_COMM_WORLD);
    }

    //Reallocate block_row array in root process to the new size
    if(rank==0){
        uint32_t* new_block_row = (uint32_t*)malloc(B->nnz_blocks*sizeof(uint32_t));
        if(new_block_row == NULL){
            printf("Couldn't allocate memory for new_block_row in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }
        memcpy(new_block_row,B->block_row,B->nnz_blocks*sizeof(uint32_t));
        free(B->block_row);
        //B->row = NULL; //isws xreastei?
        B->block_row = new_block_row;
    }

    //Send block_col array
    if(rank==0){
        for(int i=1;i<numtasks;++i){
            MPI_Isend(&B->block_col[block_col_displs[i]+1],block_col_sendcounts[i]-1,MPI_UINT32_T,i,0,MPI_COMM_WORLD,&reqs[i-1]);
        }

        MPI_Waitall(numtasks-1,reqs,MPI_STATUSES_IGNORE);

        //Reallocate block_col array in root process to the new size
        uint32_t* new_block_col = (uint32_t*)malloc((B->n_b+1)*sizeof(uint32_t));
        if(new_block_col == NULL){
            printf("Couldn't allocate memory for new_block_col in transmit_blocked_B_chunks, process %d\n",rank);
            exit(-1);
        }
        memcpy(new_block_col,B->block_col,(B->n_b+1)*sizeof(uint32_t));
        free(B->block_col);
        //B->col = NULL //isws na xreiastei?
        B->block_col = new_block_col;
    }
    else{
        MPI_Recv(&B->block_col[1],block_col_sendcounts[rank]-1,MPI_UINT32_T,0,0,MPI_COMM_WORLD,&status);
    }

    //Send the extra element missing from the block_col arrays (except the one in root process)
    if(rank!=numtasks-1){
        int next = rank+1;
        
        MPI_Isend(&B->block_col[B->n_b],1,MPI_UINT32_T,next,1,MPI_COMM_WORLD,&req_extra[1]);
    }

    if(rank!=0){
        int prev = rank-1;
        MPI_Irecv(&B->block_col[0],1,MPI_UINT32_T,prev,1,MPI_COMM_WORLD,&req_extra[0]);
        MPI_Wait(&req_extra[0],MPI_STATUSES_IGNORE);

        int offset = B->block_col[0];
        for(int i=0;i<B->n_b+1;++i){
            B->block_col[i] -= offset;
        }

    }

    MPI_Request blocks_req[2];
    MPI_Request blocks_req_2[2];

    //Send blocks
    if(rank!=0){
        for(int i=0;i<B->nnz_blocks;++i){
            MPI_Irecv(&B->blocks[i]->nnz,1,MPI_UINT32_T,0,4*i,MPI_COMM_WORLD,&blocks_req[0]);
            MPI_Irecv(&B->blocks[i]->n,1,MPI_UINT32_T,0,4*i+1,MPI_COMM_WORLD,&blocks_req[1]);

            MPI_Waitall(2,blocks_req,MPI_STATUSES_IGNORE);

            B->blocks[i]->col = (uint32_t*)malloc((B->blocks[i]->n+1)*sizeof(uint32_t));
            if(B->blocks[i]->col == NULL){
                printf("Couldn't allocate memory for B->blocks[%d]->col in transmit_blocked_B_chunks, process %d.\n",i,rank);
                exit(-1);
            }

            B->blocks[i]->row = (uint32_t*)malloc(B->blocks[i]->nnz*sizeof(uint32_t));
            if(B->blocks[i]->row == NULL){
                printf("Couldn't allocate memory for B->blocks[%d]->row in transmit_blocked_B_chunks, process %d.\n",i,rank);
                exit(-1);
            }

            if(i>0){
                MPI_Waitall(2,blocks_req_2,MPI_STATUSES_IGNORE);
            }

            MPI_Irecv(B->blocks[i]->col,B->blocks[i]->n+1,MPI_UINT32_T,0,4*i+2,MPI_COMM_WORLD,&blocks_req_2[0]);
            MPI_Irecv(B->blocks[i]->row,B->blocks[i]->nnz,MPI_UINT32_T,0,4*i+3,MPI_COMM_WORLD,&blocks_req_2[1]);

        }
    }

    MPI_Request block_send_reqs[4];
    
    if(rank==0){
        int nnz_blocks_count = B->nnz_blocks;

        for(int i=1;i<numtasks;++i){
            for(int j=0;j<block_row_sendcounts[i];++j){
                MPI_Isend(&B->blocks[nnz_blocks_count]->nnz,1,MPI_UINT32_T,i,4*j,MPI_COMM_WORLD,&block_send_reqs[0]);
                MPI_Isend(&B->blocks[nnz_blocks_count]->n,1,MPI_UINT32_T,i,4*j+1,MPI_COMM_WORLD,&block_send_reqs[1]);
                MPI_Isend(B->blocks[nnz_blocks_count]->col,B->blocks[nnz_blocks_count]->n+1,MPI_UINT32_T,i,4*j+2,MPI_COMM_WORLD,&block_send_reqs[2]);
                MPI_Isend(B->blocks[nnz_blocks_count]->row,B->blocks[nnz_blocks_count]->nnz,MPI_UINT32_T,i,4*j+3,MPI_COMM_WORLD,&block_send_reqs[3]);
                nnz_blocks_count++;
                MPI_Waitall(4,block_send_reqs,MPI_STATUSES_IGNORE);
            }
        }

        //Realloc blocks
        for(int i=B->nnz_blocks;i<total_nnz_blocks;++i){
            free_comp_matrix(B->blocks[i]);
        }
        B->blocks = (comp_matrix**)realloc(B->blocks,B->nnz_blocks*sizeof(comp_matrix*));
        if(B->blocks == NULL){
            printf("Couldn't reallocate mmemory for B->blocks in transmit_blocked_B_chunks.\n");
            exit(-1);
        }
    }

    if(B->nnz_blocks == 0){
        free_block_comp_matrix(B);
        B = NULL;
    }

    free(block_col_sendcounts);
    free(block_col_displs);
    free(block_row_sendcounts);
    free(block_row_displs);

    return B;
}


block_comp_matrix* concat_blocked_C_chunks(block_comp_matrix* result, int rank, int numtasks, int n_b_total){

    MPI_Request send_reqs[2];
    MPI_Request block_send_reqs[4];

    MPI_Request rec_reqs[2];
    MPI_Request block_rec_reqs[2];
    MPI_Request block_rec_reqs_2[2];

    uint32_t* chunks_nnz_blocks = NULL;
    uint32_t* chunks_n_b = NULL;
    uint32_t* offset = NULL;

    block_comp_matrix* C = NULL;

    if(rank==0){
        //n_b of each chunk
        chunks_n_b = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(chunks_n_b == NULL){
            printf("Couldn't allocate memory for chunks_n_b in concat_blocked_C_chunks.\n");
            exit(-1);
        }
        for(int i=0;i<numtasks;++i){
            if(i<numtasks-(n_b_total%numtasks)){
                chunks_n_b[i] = n_b_total/numtasks;
            }
            else{
                chunks_n_b[i] = n_b_total/numtasks + 1;
            }
        }

        chunks_nnz_blocks = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(chunks_nnz_blocks == NULL){
            printf("Couldn't allocate memory for chunks_nnz_blocks in concat_blocked_C_chunks.\n");
            exit(-1);
        }
    }

    if(result == NULL){
        result = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(result == NULL){
            printf("Couldn't allocate memory for result in concat_blocked_C_chunks process %d.\n",rank);
            exit(-1);
        }

        result->nnz_blocks = 0;

        if(rank<numtasks-(n_b_total%numtasks)){
            result->n_b = n_b_total/numtasks;
        }
        else{
            result->n_b = n_b_total/numtasks + 1;
        }

        result->block_col = NULL;
        result->blocks = NULL;
        result->block_row = (uint32_t*)calloc(n_b_total+1,sizeof(uint32_t));
        if(result->block_row == NULL){
            printf("Couldn't allocate memory for result->block_row in concat_blocked_C_chunks process %d.\n",rank);
            exit(-1);
        }
    }

    //nnz blocks of each chunk
    MPI_Gather(&result->nnz_blocks,1,MPI_UINT32_T,chunks_nnz_blocks,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    //Send block cols, block rows and blocks
    if(rank!=0){
        MPI_Isend(result->block_row,n_b_total+1,MPI_UINT32_T,0,2*rank,MPI_COMM_WORLD,&send_reqs[0]);
        MPI_Isend(result->block_col,result->nnz_blocks,MPI_UINT32_T,0,2*rank+1,MPI_COMM_WORLD,&send_reqs[1]);

        MPI_Waitall(2,send_reqs,MPI_STATUSES_IGNORE);

        for(int i=0;i<result->nnz_blocks;++i){
            MPI_Isend(&result->blocks[i]->n,1,MPI_UINT32_T,0,4*i+100,MPI_COMM_WORLD,&block_send_reqs[0]);
            MPI_Isend(&result->blocks[i]->nnz,1,MPI_UINT32_T,0,4*i+101,MPI_COMM_WORLD,&block_send_reqs[1]);
            MPI_Isend(result->blocks[i]->row,result->blocks[i]->n+1,MPI_UINT32_T,0,4*i+102,MPI_COMM_WORLD,&block_send_reqs[2]);
            MPI_Isend(result->blocks[i]->col,result->blocks[i]->nnz,MPI_UINT32_T,0,4*i+103,MPI_COMM_WORLD,&block_send_reqs[3]);
            MPI_Waitall(4,block_send_reqs,MPI_STATUSES_IGNORE);
        }

        free_block_comp_matrix(result);
    }

    if(rank==0){

        block_comp_matrix** chunks = (block_comp_matrix**)malloc(numtasks*sizeof(block_comp_matrix*));
        if(chunks == NULL){
            printf("Couldn't allocate memory for chunks in concat_blocked_C_chunks.\n");
            exit(-1);
        }

        for(int i=0;i<numtasks;++i){
            chunks[i] = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
            if(chunks[i] == NULL){
                printf("Couldn't allocate memory for chunks[%d] in concat_blocked_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->block_row = (uint32_t*)malloc((n_b_total+1)*sizeof(uint32_t));
            if(chunks[i]->block_row == NULL){
                printf("Couldn't allocate memory for chunks[%d]->block_row in concat_blocked_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->block_col = (uint32_t*)malloc(chunks_nnz_blocks[i]*sizeof(uint32_t));
            if(chunks[i]->block_col == NULL){
                printf("Couldn't allocate memory for chunks[%d]->block_col in concat_blocked_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->blocks = (comp_matrix**)malloc(chunks_nnz_blocks[i]*sizeof(comp_matrix*));
            if(chunks[i]->blocks == NULL){
                printf("Couldn't allocate memory for chunks[%d]->blocks in concat_blocked_C_chunks.\n",i);
                exit(-1);
            }

            chunks[i]->nnz_blocks = chunks_nnz_blocks[i];
            chunks[i]->n_b = chunks_n_b[i];
            chunks[i]->real_dim = result->real_dim;

            if(i==0){
                memcpy(chunks[i]->block_col,result->block_col,chunks[i]->nnz_blocks*sizeof(uint32_t));
                memcpy(chunks[i]->block_row,result->block_row,(n_b_total+1)*sizeof(uint32_t));

                for(int j=0;j<chunks[i]->nnz_blocks;++j){
                    chunks[i]->blocks[j] = (comp_matrix*)malloc(sizeof(comp_matrix));
                    if(chunks[i]->blocks[j] == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d] in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    chunks[i]->blocks[j]->nnz = result->blocks[j]->nnz;
                    chunks[i]->blocks[j]->n = result->blocks[j]->n;
                    
                    chunks[i]->blocks[j]->col = (uint32_t*)malloc(chunks[i]->blocks[j]->nnz*sizeof(uint32_t));
                    if(chunks[i]->blocks[j]->col == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d]->col in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    chunks[i]->blocks[j]->row = (uint32_t*)malloc((chunks[i]->blocks[j]->n+1)*sizeof(uint32_t));
                    if(chunks[i]->blocks[j]->row == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d]->row in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    memcpy(chunks[i]->blocks[j]->col,result->blocks[j]->col,chunks[i]->blocks[j]->nnz*sizeof(uint32_t));
                    memcpy(chunks[i]->blocks[j]->row,result->blocks[j]->row,(chunks[i]->blocks[j]->n+1)*sizeof(uint32_t));
                }
            }
            else{
                MPI_Irecv(chunks[i]->block_row,n_b_total+1,MPI_UINT32_T,i,2*i,MPI_COMM_WORLD,&rec_reqs[0]);
                MPI_Irecv(chunks[i]->block_col,chunks[i]->nnz_blocks,MPI_UINT32_T,i,2*i+1,MPI_COMM_WORLD,&rec_reqs[1]);

                MPI_Waitall(2,rec_reqs,MPI_STATUSES_IGNORE);

                for(int j=0;j<chunks[i]->nnz_blocks;++j){
                    chunks[i]->blocks[j] = (comp_matrix*)malloc(sizeof(comp_matrix));
                    if(chunks[i]->blocks[j] == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d] in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    MPI_Irecv(&chunks[i]->blocks[j]->n,1,MPI_UINT32_T,i,4*j+100,MPI_COMM_WORLD,&block_rec_reqs[0]);
                    MPI_Irecv(&chunks[i]->blocks[j]->nnz,1,MPI_UINT32_T,i,4*j+101,MPI_COMM_WORLD,&block_rec_reqs[1]);
                    
                    MPI_Waitall(2,block_rec_reqs,MPI_STATUSES_IGNORE);

                    chunks[i]->blocks[j]->col = (uint32_t*)malloc(chunks[i]->blocks[j]->nnz*sizeof(uint32_t));
                    if(chunks[i]->blocks[j]->col == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d]->col in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    chunks[i]->blocks[j]->row = (uint32_t*)malloc((chunks[i]->blocks[j]->n+1)*sizeof(uint32_t));
                    if(chunks[i]->blocks[j]->row == NULL){
                        printf("Couldn't allocate memory for chunks[%d]->blocks[%d]->row in concat_blocked_C_chunks.\n",i,j);
                        exit(-1);
                    }

                    MPI_Irecv(chunks[i]->blocks[j]->row,chunks[i]->blocks[j]->n+1,MPI_UINT32_T,i,4*j+102,MPI_COMM_WORLD,&block_rec_reqs_2[0]);
                    MPI_Irecv(chunks[i]->blocks[j]->col,chunks[i]->blocks[j]->nnz,MPI_UINT32_T,i,4*j+103,MPI_COMM_WORLD,&block_rec_reqs_2[1]);
                
                    MPI_Waitall(2,block_rec_reqs_2,MPI_STATUSES_IGNORE);
                }
            }
        }

        free_block_comp_matrix(result);
        free(chunks_nnz_blocks);
        free(chunks_n_b); 

        //Create the final matrix to be returned
        uint32_t total_nnz_blocks = 0;     //Total number of non-zero blocks in the whole matrix
        for(int i=0;i<numtasks;++i){
            total_nnz_blocks += chunks[i]->nnz_blocks;
        }

        if(total_nnz_blocks==0){
            return NULL;
        }

        C = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(C == NULL){
            printf("Couldn't allocate memory for C in concat_blocked_C_chunks.\n");
            exit(-1);
        }

        C->n_b = n_b_total;
        C->nnz_blocks = total_nnz_blocks;
        C->real_dim = chunks[rank]->real_dim;

        C->block_row = (uint32_t*)malloc((C->n_b+1)*sizeof(uint32_t));
        if(C->block_row == NULL){
            printf("Couldn't allocate memory for C->block_row in concat_blocked_C_chunks.\n");
            exit(-1);
        }

        C->block_col = (uint32_t*)malloc(C->nnz_blocks*sizeof(uint32_t));
        if(C->block_col == NULL){
            printf("Couldn't allocate memory for C->block_col in concat_blocked_C_chunks.\n");
            exit(-1);
        }

        C->blocks = (comp_matrix**)malloc(total_nnz_blocks*sizeof(comp_matrix*));
        if(C->blocks == NULL){
            printf("Couldn't allocate memory for blocks in concat_blocked_C_chunks.\n");
            exit(-1);
        }
        uint32_t nnz_block_count = 0;   //Number of non-zero blocks found up to this point
        
        uint32_t block_row_in_chunk_start, block_row_in_chunk_end;

        offset = (uint32_t*)malloc(numtasks*sizeof(uint32_t));
        if(offset == NULL){
            printf("Couldn't allocate memory for offset in concat_blocked_C_chunks.\n");
            exit(-1);
        }
        offset[0] = 0;
        for(int i=1;i<numtasks;++i){
            offset[i] = offset[i-1]+chunks[i-1]->n_b;
        }

        C->block_row[0] = 0;

        for(uint32_t i=0;i<n_b_total;++i){
            for(uint32_t j=0;j<numtasks;++j){

                if(chunks[j]->nnz_blocks == 0){
                    continue;
                }

                block_row_in_chunk_start = chunks[j]->block_row[i];
                block_row_in_chunk_end = chunks[j]->block_row[i+1];

                for(uint32_t k=block_row_in_chunk_start;k<block_row_in_chunk_end;++k){
                    C->block_col[nnz_block_count] = chunks[j]->block_col[k] + offset[j];

                    C->blocks[nnz_block_count] = (comp_matrix*)malloc(sizeof(comp_matrix));
                    if(C->blocks[nnz_block_count] == NULL){
                        printf("Couldn't allocate memory for C->blocks[%d] in concat_blocked_C_chunks.\n",nnz_block_count);
                        exit(-1);
                    }

                    C->blocks[nnz_block_count]->n = chunks[j]->blocks[k]->n;
                    C->blocks[nnz_block_count]->nnz = chunks[j]->blocks[k]->nnz;

                    C->blocks[nnz_block_count]->row = (uint32_t*)malloc((C->blocks[nnz_block_count]->n+1)*sizeof(uint32_t));
                    if(C->blocks[nnz_block_count]->row == NULL){
                        printf("Couldn't allocate memory for C->blocks[%d]->row in concat_blocked_C_chunks.\n",nnz_block_count);
                        exit(-1);
                    }

                    C->blocks[nnz_block_count]->col = (uint32_t*)malloc(C->blocks[nnz_block_count]->nnz*sizeof(uint32_t));
                    if(C->blocks[nnz_block_count]->col == NULL){
                        printf("Couldn't allocate memory for C->blocks[%d]->col in concat_blocked_C_chunks.\n",nnz_block_count);
                        exit(-1);
                    }

                    memcpy(C->blocks[nnz_block_count]->row,chunks[j]->blocks[k]->row,(C->blocks[nnz_block_count]->n+1)*sizeof(uint32_t));

                    for(uint32_t l=0;l<C->blocks[nnz_block_count]->nnz;++l){
                        C->blocks[nnz_block_count]->col[l] = chunks[j]->blocks[k]->col[l] + offset[j]*C->blocks[nnz_block_count]->n;
                    }
                    nnz_block_count++;
                }

            }

            C->block_row[i+1] = nnz_block_count;
        }

        free(offset);

        for(int i=0;i<numtasks;++i){
            free_block_comp_matrix(chunks[i]);
        }
        free(chunks);
    }

    return C;
}


block_comp_matrix* blocked_bmm_parallel_2(block_comp_matrix* A, block_comp_matrix* B, int rank){

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    A = transmit_blocked_A(A,rank);

    B = transmit_blocked_B_chunks(B,rank);

    block_comp_matrix* C_chunk;

    MPI_Barrier(MPI_COMM_WORLD);

    if(B == NULL){
        C_chunk = NULL;
    }
    else{
        C_chunk = blocked_bmm_seq(A,B);
    }

    block_comp_matrix* C = concat_blocked_C_chunks(C_chunk,rank,numtasks,A->n_b);

    free_block_comp_matrix(A);
    free_block_comp_matrix(B);

    return C;
    
}


/*----------- Functions for the non-blocked parallel filtered implementation -----------*/

comp_matrix* snip_F(comp_matrix* F, uint32_t total_cols){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    uint32_t min_col;
    uint32_t max_col;

    int remaining = total_cols%numtasks;
    
    if (rank < numtasks-remaining){
        min_col = rank * total_cols/numtasks;
        max_col = (rank + 1) * total_cols/numtasks;
    }
    else{
        min_col = (numtasks-remaining)*total_cols/numtasks + (rank-numtasks+remaining)*(total_cols/numtasks+1);
        max_col = (numtasks-remaining)*total_cols/numtasks + (rank+1-numtasks+remaining)*(total_cols/numtasks+1);
    }

    uint32_t* temp_row = (uint32_t*)malloc((F->n+1) * sizeof(uint32_t));
    if(temp_row == NULL){
        printf("Couldn't allocate memory for temp_row in snip_F, process %d\n",rank);
        exit(-1);
    }
    uint32_t* temp_col = (uint32_t*)malloc(F->nnz * sizeof(uint32_t));
    if(temp_col == NULL){
        printf("Couldn't allocate memory for temp_col in snip_F, process %d\n",rank);
        exit(-1);
    }

    uint32_t snipped_nnz = 0;

    for (int i=0; i<F->n; i++){
        temp_row[i] = snipped_nnz;
        for (int j=F->row[i]; j<F->row[i+1]; j++){
            if (F->col[j]<max_col && F->col[j]>=min_col){
                temp_col[snipped_nnz] = F->col[j]-min_col;
                snipped_nnz++;
            }
        }
    }
    temp_row[F->n] = snipped_nnz;

    free(F->col);
    free(F->row);

    temp_col = (uint32_t*)realloc(temp_col,snipped_nnz*sizeof(uint32_t));
    if(temp_col == NULL){
        printf("Couldn't reallocate memory for temp_col in snip_F.\n");
        exit(-1);
    }

    F->col = temp_col;
    F->row = temp_row;
    F->nnz = snipped_nnz;

    return F;
    
}

//Transfers matrix F to all processes
comp_matrix* transmit_F(comp_matrix* F, int rank, uint32_t total_cols){

    if(rank!=0){
        F = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(F == NULL){
            printf("Couldn't allocate memory for F in transmit_F, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(&F->n,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&F->nnz,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    
    if(rank!=0){
        F->col = (uint32_t*)malloc(F->nnz*sizeof(uint32_t));
        if(F->col == NULL){
            printf("Couldn't allocate memory for F->col in transmit_F, process %d\n",rank);
            exit(-1);
        }

        F->row = (uint32_t*)malloc((F->n+1)*sizeof(uint32_t));
        if(F->row == NULL){
            printf("Couldn't allocate memory for F->row in transmit_F, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(F->col,F->nnz,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(F->row,F->n+1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    F = snip_F(F,total_cols);

    return F;
}



comp_matrix* nonblocked_bmm_parallel_filtered_2(comp_matrix* A, comp_matrix* B, comp_matrix* F, int rank){

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
        C_chunk = bmm_filtered_seq(A,B,F,0);
    }

    comp_matrix* C = concat_C_chunks(C_chunk,rank,numtasks,F->n);

    free_comp_matrix(A);
    free_comp_matrix(B);
    free_comp_matrix(F);    

    return C;
}



/* ------------ Functions for the blocked parallel filtered implementation ------------- */

block_comp_matrix* snip_blocked_F(block_comp_matrix* F, uint32_t total_block_cols){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    uint32_t min_block_col;
    uint32_t max_block_col;

    int remaining = total_block_cols%numtasks;
    
    if (rank < numtasks-remaining){
        min_block_col = rank * total_block_cols/numtasks;
        max_block_col = (rank + 1) * total_block_cols/numtasks;
    }
    else{
        uint32_t previous_block_cols = (numtasks-remaining)*total_block_cols/numtasks;
        min_block_col = previous_block_cols + (rank-numtasks+remaining)*(total_block_cols/numtasks+1);
        max_block_col = previous_block_cols + (rank+1-numtasks+remaining)*(total_block_cols/numtasks+1);
    }

    uint32_t* temp_block_row = (uint32_t*)malloc((F->n_b+1) * sizeof(uint32_t));
    if(temp_block_row == NULL){
        printf("Couldn't allocate memory for temp_block_row in snip_blocked_F, process %d\n",rank);
        exit(-1);
    }
    uint32_t* temp_block_col = (uint32_t*)malloc(F->nnz_blocks * sizeof(uint32_t));
    if(temp_block_col == NULL){
        printf("Couldn't allocate memory for temp_block_col in snip_blocked_F, process %d\n",rank);
        exit(-1);
    }
    comp_matrix** temp_blocks = (comp_matrix**)malloc(F->nnz_blocks*sizeof(comp_matrix*));
    if(temp_blocks == NULL){
        printf("Couldn't allocate memory for temp_blocks in snip_blocked_F, process %d\n",rank);
        exit(-1);
    }

    uint32_t snipped_nnz_blocks = 0;

    for (int i=0; i<F->n_b; i++){
        temp_block_row[i] = snipped_nnz_blocks;

        for (int j=F->block_row[i]; j<F->block_row[i+1]; j++){
            if (F->block_col[j]<max_block_col && F->block_col[j]>=min_block_col){
                temp_block_col[snipped_nnz_blocks] = F->block_col[j];

                temp_blocks[snipped_nnz_blocks] = new_comp_matrix(F->blocks[j]->nnz, F->blocks[j]->n, "csr");

                memcpy(temp_blocks[snipped_nnz_blocks]->col, F->blocks[j]->col, F->blocks[j]->nnz*sizeof(uint32_t));
                memcpy(temp_blocks[snipped_nnz_blocks]->row, F->blocks[j]->row, (F->blocks[j]->n+1)*sizeof(uint32_t));

                snipped_nnz_blocks++;
            }
        }
    }
    temp_block_row[F->n_b] = snipped_nnz_blocks;

    free(F->block_col);
    free(F->block_row);
    for (int i=0; i<F->nnz_blocks; i++){
        free_comp_matrix(F->blocks[i]);
    }
    free(F->blocks);

    temp_block_col = (uint32_t*)realloc(temp_block_col,snipped_nnz_blocks*sizeof(uint32_t));
    if(temp_block_col == NULL){
        printf("Couldn't reallocate memory for temp_block_col in snip_blocked_F.\n");
        exit(-1);
    }
    temp_blocks = (comp_matrix**)realloc(temp_blocks,snipped_nnz_blocks*sizeof(comp_matrix*));
    if(temp_blocks == NULL){
        printf("Couldn't reallocate memory for temp_blocks in snip_blocked_F.\n");
        exit(-1);
    }

    F->block_col = temp_block_col;
    F->block_row = temp_block_row;
    F->blocks = temp_blocks;
    F->nnz_blocks = snipped_nnz_blocks;

    return F;
    
}


block_comp_matrix* transmit_blocked_F(block_comp_matrix* F, int rank, uint32_t total_block_cols){

    if(rank!=0){
        F = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(F == NULL){
            printf("Couldn't allocate memory for F in transmit_blocked_F, process %d\n",rank);
            exit(-1);
        }
    }

    MPI_Bcast(&F->n_b,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&F->nnz_blocks,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&F->real_dim,1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    if(rank!=0){
        F->block_col = (uint32_t*)malloc(F->nnz_blocks*sizeof(uint32_t));
        if(F->block_col == NULL){
            printf("Couldn't allocate memory for F->block_col in transmit_blocked_F, process %d\n",rank);
            exit(-1);
        }

        F->block_row = (uint32_t*)malloc((F->n_b+1)*sizeof(uint32_t));
        if(F->block_row == NULL){
            printf("Couldn't allocate memory for F->block_row in transmit_blocked_F, process %d\n",rank);
            exit(-1);
        }

        F->blocks = (comp_matrix**)malloc(F->nnz_blocks*sizeof(comp_matrix*));
        if(F->blocks == NULL){
            printf("Couldn't allocate memory for F->blocks in transmit_blocked_F, process %d\n",rank);
            exit(-1);
        }
        for(int i=0;i<F->nnz_blocks;++i){
            F->blocks[i] = NULL;
        }
    }

    MPI_Bcast(F->block_col,F->nnz_blocks,MPI_UINT32_T,0,MPI_COMM_WORLD);
    MPI_Bcast(F->block_row,F->n_b+1,MPI_UINT32_T,0,MPI_COMM_WORLD);

    for(int i=0;i<F->nnz_blocks;++i){
        F->blocks[i] = transmit_A(F->blocks[i],rank);
    }

    F = snip_blocked_F(F,total_block_cols);

    return F;
}



block_comp_matrix* blocked_bmm_parallel_filtered_2(block_comp_matrix* A, block_comp_matrix* B, block_comp_matrix* F, int rank){

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
        C_chunk = blocked_bmm_seq_filtered(A,B,F,offset);
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