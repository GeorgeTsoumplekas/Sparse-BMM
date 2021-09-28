#include <omp.h>

#include "utilities.h"
#include "bmm_seq.h"


/*----------- Functions for the non-blocked parallel implementation -----------*/


/**
 * This function concatenates the csr matrices held in chunks into one csr matrix. Each chunk
 * is a csr matrix that has been created by a different thread and represents a different part
 * of the final matrix. For example, if we have deployed 4 threads for the matrix multiplication
 * then chunks[0] has the 1st 1/4th of the product, chunks[1] the 2nd 1/4th and so forth. By
 * combining them all together we get the final product matrix.
**/
comp_matrix* concat_chunks(comp_matrix** chunks, int chunk_num, uint32_t n){

    uint32_t total_nnz = 0;     //Total number of non-zero elements in the whole matrix
    for(int i=0;i<chunk_num;++i){
        total_nnz += chunks[i]->nnz;
    }

    if(total_nnz==0){
        return NULL;
    }

    comp_matrix* total = new_comp_matrix(total_nnz,n,"csr");

    uint32_t nnz_count = 0;
    
    //Fill in the col array of the final matrix
    for(uint32_t i=0;i<chunk_num;++i){
        //If it is a zero elements only chunk
        if(chunks[i]==NULL){
            continue;
        }
        for(uint32_t j=0;j<chunks[i]->nnz;++j){
            total->col[nnz_count] = chunks[i]->col[j];
            nnz_count++;
        }
    }

    uint32_t row_count = 0; //Number of rows we have examined up to this point

    //Fill in the row array of the final matrix
    for(uint32_t i=0;i<chunk_num;++i){
        //If it is a zero elements only chunk
        if(chunks[i]->nnz==0){
            //If it is the first chunk then fill in these rows with 0
            if(i==0){
                for(uint32_t j=0;j<chunks[i]->n+1;++j){
                    total->row[j] = 0;
                }
                row_count = chunks[i]->n;
            }
            //If it's another chunk fill in these rows with the number of elements up to this point
            else{
                for(uint32_t j=0;j<chunks[i]->n+1;++j){
                    total->row[row_count+j] = total->row[row_count];
                }
                row_count += chunks[i]->n;
            }
        }
        else{
            for(uint32_t j=0;j<chunks[i]->n+1;++j){
                total->row[row_count+j] = total->row[row_count]+chunks[i]->row[j];
            }
            row_count += chunks[i]->n;
        }
    }

    return total;    
}



/**
 * Function that performs non-blocked BMM using thread parallelization with OpenMP.
 * Matrix A is split into chunks and each chunk is assigned to adifferent thread.
 * The split is done row-wise and each chunk has (almost) the same number of rows.
 * Then each chunk is multiplied with matrix B and a chunk of the final product matrix is created.
 * By cocatenating all these chunks of the final matrix we get the whole final matrix.
**/
comp_matrix* nonblocked_bmm_parallel(comp_matrix* A, comp_matrix* B, int n_threads){

    omp_set_num_threads(n_threads);

    uint32_t n = A->n;

    //Number of rows of A that each thread will examine
    int* chunk_size = (int*)malloc(n_threads*sizeof(int));
    if(chunk_size == NULL){
        printf("Couldn;t allocate memory for chunk_size in nonblocked_bmm_parallel.\n");
        exit(-1);
    }

    //If the size of the matrix is not exactly divisible by the number of threads
    //then give an extra row to an appropriate number of threads (this is to achieve better
    //load balance)
    uint32_t res = n % n_threads;
    for(int i=0;i<res;++i){
        chunk_size[i] = n/n_threads + 1;
    }
    for(int i=res;i<n_threads;++i){
        chunk_size[i] = n/n_threads;
    }

    //matrix containing the chunks (each chunk is represented as a csr matrix)
    comp_matrix** csr_chunks = (comp_matrix**)malloc(n_threads*sizeof(comp_matrix*));
    if(csr_chunks == NULL){
        printf("Couldn't allocate memory for csr_chunks in nonblocked_bmm_parallel.\n");
        exit(-1);
    }   

    //Beacause we don't know the exact number of non-zero elements that each chunk of the 
    //product will have, we give it an initial value proportional to the sum of non-zero
    //elements in A and B and if necessary we will resize it later
    int32_t nnz_initial = (A->nnz+B->nnz)/n_threads;
    printf("%d\n",nnz_initial);

    for(int i=0;i<n_threads;++i){
        csr_chunks[i] = new_comp_matrix(nnz_initial,chunk_size[i],"csr");
    }

    int i,j;    //Iterators used to run through the rows of A and columns of B
    uint32_t nnz_count = 0; //Non-zero elements found up to this point


    #pragma omp parallel shared(A,B,csr_chunks) private(i,j) firstprivate(nnz_count)
    {
        int thread_id = omp_get_thread_num();

        int i_start;    //Index of the first row that this thread will examine

        if(thread_id>=res){
            i_start = res*(n/n_threads+1) + (thread_id-res)*(n/n_threads);
        }
        else{
            i_start = thread_id*(n/n_threads+1);
        }

        uint32_t A_row_ptr, A_row_end;
        uint32_t B_col_ptr, B_col_end;

        //For every row in the chunk
        for(i=0;i<chunk_size[thread_id];++i){
            A_row_ptr = A->row[i_start+i];
            A_row_end = A->row[i_start+i+1];

            csr_chunks[thread_id]->row[i] = nnz_count;

            //For every column in B
            for(j=0;j<B->n;++j){
                B_col_ptr = B->col[j];
                B_col_end = B->col[j+1];

                //While there are still non-zero elements in the corresponding row and column
                while (A_row_ptr < A_row_end && B_col_ptr < B_col_end) {
                    //If A's current element's column index is smaller than B's current element's row index
                    if (A->col[A_row_ptr] < B->row[B_col_ptr]) {
                        A_row_ptr++;
                    }
                    //if A's current element's column index is bigger than B's current element's row index
                    else if (A->col[A_row_ptr] > B->row[B_col_ptr]) {
                        B_col_ptr++;
                    }
                    //If A's current element's column index is matching with B's current element's row index
                    //then we have a non-zero element in the product chunk in this position
                    else{
                        csr_chunks[thread_id]->col[nnz_count] = j;
                        nnz_count++;

                        //Extend the col array of the product chunk if there are not enough empty cells
                        if (csr_chunks[thread_id]->nnz <= nnz_count) {
                            csr_chunks[thread_id]->col = realloc(csr_chunks[thread_id]->col,2*csr_chunks[thread_id]->nnz*sizeof(uint32_t));
                            if(csr_chunks[thread_id]->col == NULL){
                                printf("Couldn't reallocate memory for csr_chunks[%d]->col in nonblocked_bmm_parallel.\n",thread_id);
                                exit(-1);
                            }
                            csr_chunks[thread_id]->nnz = 2*csr_chunks[thread_id]->nnz;
                        }

                        //If we have already found a non-zero element in this position,
                        //there is no reason to continue calculating for this row of A and column of B.
                        break;
                    }
                }

                A_row_ptr = A->row[i_start+i];
            }
        }

        //In case the multiplication gives a matrix with zero elements only, retun a NULL matrix
        if(nnz_count == 0){
            free_comp_matrix(csr_chunks[thread_id]);
            csr_chunks[thread_id] = NULL;
        }

        //Change to the true number of non zero values
        csr_chunks[thread_id]->nnz = nnz_count;

        csr_chunks[thread_id]->row[csr_chunks[thread_id]->n] = nnz_count;

        //Resize col index array to correct size
        csr_chunks[thread_id]->col = (uint32_t *)realloc(csr_chunks[thread_id]->col,csr_chunks[thread_id]->nnz*sizeof(uint32_t));
    }

    //Concatenate the csr matrices of each chunk into one for the whole product matrix
    comp_matrix* C = concat_chunks(csr_chunks,n_threads,A->n);

    for(uint32_t i=0;i<n_threads;++i){
        free_comp_matrix(csr_chunks[i]);
    }
    free(csr_chunks);
    free(chunk_size);

    return C;
}


/* ----------------- Functions for the blocked parallel implementation -------------------- */


/**
 * This function concatenates the csr matrices held in chunks into one csr matrix. 
 * It works the same way as the concat_chunks function above with the only difference that
 * the chunks are now blocked csr matrices and the final matrix is also in blocked csr format.
**/
block_comp_matrix* blocked_concat_chunks(block_comp_matrix** chunks, int chunk_num, uint32_t n_b, uint32_t real_dim){
    uint32_t total_nnz_blocks = 0;  //Total number of non-zero blocks in the whole matrix
    for(int i=0;i<chunk_num;++i){
        if(chunks[i]==NULL){
            continue;
        }
        total_nnz_blocks += chunks[i]->nnz_blocks ;
    }

    if(total_nnz_blocks==0){
        return NULL;
    }

    block_comp_matrix* total = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if(total == NULL){
        printf("Couldn't allocate memory for total in blocked_concat_chunks.\n");
        exit(-1);
    }

    total->block_col = (uint32_t*)malloc(total_nnz_blocks*sizeof(uint32_t));
    if(total->block_col == NULL){
        printf("Couldn't allocate memory for total->block_col in blocked_concat_chunks.\n");
        exit(-1);
    }

    total->block_row = (uint32_t*)malloc((n_b+1)*sizeof(uint32_t));
    if(total->block_row == NULL){
        printf("Couldn't allocate memory for total->block_col in blocked_concat_chunks.\n");
        exit(-1);
    }

    total->blocks = (comp_matrix**)malloc(total_nnz_blocks*sizeof(comp_matrix*));
    if(total->blocks == NULL){
        printf("Couldn't allocate memory for blocks in blocked_concat_chunks.\n");
        exit(-1);
    }

    total->n_b = n_b;
    total->nnz_blocks = total_nnz_blocks;
    total->real_dim = real_dim;

    uint32_t nnz_block_count = 0;   //Number of non-zero blocks found up to this point

    //Fill in the col and blocks array of the final matrix
    for(uint32_t i=0;i<chunk_num;++i){
        //If it is a zero elements only chunk
        if(chunks[i]==NULL){
            continue;
        }
        for(uint32_t j=0;j<chunks[i]->nnz_blocks;++j){
            total->block_col[nnz_block_count] = chunks[i]->block_col[j];
            
            total->blocks[nnz_block_count] = (comp_matrix*)malloc(sizeof(comp_matrix));
            if(total->blocks[nnz_block_count] == NULL){
                printf("Couldn't allocate memory for blocks[%d] in blocked_concat_chunks.\n",nnz_block_count);
                exit(-1);
            }

            total->blocks[nnz_block_count]->col = (uint32_t*)malloc(chunks[i]->blocks[j]->nnz*sizeof(uint32_t));
            if(total->blocks[nnz_block_count]->col == NULL){
                printf("Couldn't allocate memory for blocks[%d]->col in blocked_concat_chunks.\n",nnz_block_count);
                exit(-1);
            }

            total->blocks[nnz_block_count]->row = (uint32_t*)malloc((chunks[i]->blocks[j]->n+1)*sizeof(uint32_t));
            if(total->blocks[nnz_block_count]->row == NULL){
                printf("Couldn't allocate memory for blocks[%d]->row in blocked_concat_chunks.\n",nnz_block_count);
                exit(-1);
            }

            //Copy the examined non-zero block in the chunk to the blocks array of the whole matrix
            total->blocks[nnz_block_count]->nnz = chunks[i]->blocks[j]->nnz;
            total->blocks[nnz_block_count]->n = chunks[i]->blocks[j]->n;

            memcpy(total->blocks[nnz_block_count]->col,chunks[i]->blocks[j]->col,chunks[i]->blocks[j]->nnz*sizeof(uint32_t));
            memcpy(total->blocks[nnz_block_count]->row,chunks[i]->blocks[j]->row,(chunks[i]->blocks[j]->n+1)*sizeof(uint32_t));
            nnz_block_count++;
        }
    }

    uint32_t block_row_count = 0;   //Number of rows of blocks we have examined up to this point

    total->block_row[0] = 0;

    for(uint32_t i=0;i<chunk_num;++i){
        //If it is a zero elements only chunk
        if(chunks[i]->nnz_blocks==0){
            //If it is the first chunk then fill in these block rows with 0
            if(i==0){
                for(uint32_t j=0;j<chunks[i]->n_b+1;++j){
                    total->block_row[j] = 0;
                }
                block_row_count = chunks[i]->n_b;
            }
            //If it's another chunk fill in these block rows with the number of non-zero blocks up to this point
            else{
                for(uint32_t j=0;j<chunks[i]->n_b+1;++j){
                    total->block_row[block_row_count+j] = total->block_row[block_row_count];
                }
                block_row_count += chunks[i]->n_b;
            }
        }
        else{
            for(uint32_t j=0;j<chunks[i]->n_b+1;++j){
                total->block_row[block_row_count+j] = total->block_row[block_row_count]+chunks[i]->block_row[j];
            }
            block_row_count += chunks[i]->n_b;
        }
    }

    return total;    
}


/**
 * Function that performs non-blocked BMM using thread parallelization with OpenMP.
 * This function works the same way as the non-blocked parallel bmm with the difference that
 * we use blocked matrices to perform multiplication.
**/
block_comp_matrix* blocked_bmm_parallel(block_comp_matrix* A, block_comp_matrix* B, int n_threads){
    omp_set_num_threads(n_threads);

    uint32_t n_b = A->n_b;

    //Number of rows of blocks that each thread will examine
    int* chunk_size = (int*)malloc(n_threads*sizeof(int));
    if(chunk_size == NULL){
        printf("Couldn;t allocate memory for chunk_size in blocked_bmm_parallel.\n");
        exit(-1);
    }

    //If the number of block rows in the matrix is not exactly divisible by the number of threads
    //then give an extra block row to an appropriate number of threads (this is to achieve better
    //load balance)
    uint32_t res = n_b % n_threads;
    for(int i=0;i<res;++i){
        chunk_size[i] = n_b/n_threads + 1;
    }
    for(int i=res;i<n_threads;++i){
        chunk_size[i] = n_b/n_threads;
    }

    //matrix containing the chunks (each chunk is represented as a blocked csr matrix)
    block_comp_matrix** blocked_csr_chunks = (block_comp_matrix**)malloc(n_threads*sizeof(block_comp_matrix*));
    if(blocked_csr_chunks == NULL){
        printf("Couldn't allocate memory for csr_chunks in nonblocked_bmm_parallel.\n");
        exit(-1);
    }  

    //We assume A is a non-zero matrix
    uint32_t b = A->blocks[0]->n;

    for(int i=0;i<n_threads;++i){
        blocked_csr_chunks[i] = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
        if(blocked_csr_chunks[i]==NULL){
            printf("Couldn't allocate memory for blocked_csr_chunks[%d] in blocked_bmm_parallel.\n",i);
            exit(-1);
        }

        blocked_csr_chunks[i]->block_row = (uint32_t*)malloc((chunk_size[i]+1)*sizeof(uint32_t));
        if(blocked_csr_chunks[i]->block_row == NULL){
            printf("Couldn't allocate memory for block_row in blocked_bmm_parallel.\n");
            exit(-1);
        }
        
        blocked_csr_chunks[i]->block_col = (uint32_t*)malloc(chunk_size[i]*n_b*sizeof(uint32_t));
        if(blocked_csr_chunks[i]->block_col == NULL){
            printf("Couldn't allocate memory for block_col in blocked_bmm_parallel.\n");
            exit(-1);
        }

        blocked_csr_chunks[i]->blocks = (comp_matrix**)malloc(chunk_size[i]*n_b*sizeof(comp_matrix*));
        if(blocked_csr_chunks[i]->blocks == NULL){
            printf("Couldn't allocate memory for blocks in blocked_bmm_parallel.\n");
            exit(-1);
        }

        blocked_csr_chunks[i]->n_b = chunk_size[i];
        blocked_csr_chunks[i]->real_dim = A->real_dim;
    }

    //Number of the product chunks's non-zero blocks found up to this point
    uint32_t nnz_blocks_found = 0;

    //Used as a flag to see whether it is the first time we have matching blocks to apply multiplication to
    uint32_t first_match;       

    comp_matrix* prod = NULL;

    #pragma omp parallel shared(A,B,blocked_csr_chunks) private(first_match) firstprivate(nnz_blocks_found,prod)
    {
        int thread_id = omp_get_thread_num();

        int i_start;    //Index of the first block row that this thread will examine
        
        if(thread_id>=res){
            i_start = res*(n_b/n_threads+1) + (thread_id-res)*(n_b/n_threads);
        }
        else{
            i_start = thread_id*(n_b/n_threads+1);
        }

        uint32_t block_row_start;   //First non-zero block in this row of blocks
        uint32_t block_row_end;     //Last non-zero block in this row of blocks

        uint32_t block_col_start;   //First non-zero block in this column of blocks
        uint32_t block_col_end;     //Last non-zero block in this column of blocks

        uint32_t block_row_ptr;     //Iterator through the non-zero blocks of a row of blocks
        uint32_t block_col_ptr;     //Iterator through the non-zero blocks of a column of blocks

        blocked_csr_chunks[thread_id]->block_row[0] = 0;

        //For each row of blocks in the chunk
        for(uint32_t i=0;i<chunk_size[thread_id];++i){
            block_row_start = A->block_row[i_start+i];
            block_row_end = A->block_row[i_start+i+1];

            //If there aren't any non-zero blocks in this row of blocks in the chunk
            //then there aren't any non-zero blocks in this row of blocks in the
            //corresponding product chunk either
            if(block_row_start == block_row_end){
                blocked_csr_chunks[thread_id]->block_row[i+1] = blocked_csr_chunks[thread_id]->block_row[i];
                continue;
            }

            //For each column of blocks in B
            for(uint32_t j=0;j<n_b;++j){
                block_col_start = B->block_col[j];
                block_col_end = B->block_col[j+1];
                
                //If there aren't any non-zero blocks in this column go to the next one
                if(block_col_start == block_col_end){
                    continue;
                }
                
                block_row_ptr = block_row_start;
                block_col_ptr = block_col_start;

                first_match = 0;

                //While there are still non-zero blocks in the row of blocks in the chunk or the column of blocks in B
                while(block_row_ptr<block_row_end && block_col_ptr<block_col_end){
                    //Search for blocks in the chunk with the same col index as the row index of the blocks in B
                    if(A->block_col[block_row_ptr] > B->block_row[block_col_ptr]){
                        block_col_ptr++;
                    }
                    else if(A->block_col[block_row_ptr] < B->block_row[block_col_ptr]){
                        block_row_ptr++;
                    }
                    else{
                        //In this part we multiply the correspondind blocks
                        //prod is a csr matrix
                        prod = bmm_seq(A->blocks[block_row_ptr],B->blocks[block_col_ptr],j*b);

                        //If it is the first product we calculate for this block of the new blocked matrix
                        if(first_match==0){

                            //If the multiplication doesn't yield a non-zero matrix then do nothing
                            if(prod == NULL){
                                block_col_ptr++;
                                block_row_ptr++;
                                continue;
                            }

                            //Else create a new non-zero block for the new matrix
                            blocked_csr_chunks[thread_id]->blocks[nnz_blocks_found] = NULL;
                            blocked_csr_chunks[thread_id]->block_col[nnz_blocks_found] = j;
                            first_match = 1;
                            nnz_blocks_found++; 

                        }

                        //No reason to apply union if the product is NULL
                        if(prod != NULL){
                            blocked_csr_chunks[thread_id]->blocks[nnz_blocks_found-1] = block_union(blocked_csr_chunks[thread_id]->blocks[nnz_blocks_found-1],prod);
                            free_comp_matrix(prod); 
                            prod = NULL;
                        }

                        block_col_ptr++;
                        block_row_ptr++;    
            
                    }
                }
            }
            blocked_csr_chunks[thread_id]->block_row[i+1] = nnz_blocks_found;
        }

        //In case the multiplication gives a chunks with zero elements only, retun a NULL matrix
        if(nnz_blocks_found == 0){
            free_block_comp_matrix(blocked_csr_chunks[thread_id]);
            blocked_csr_chunks[thread_id] = NULL;
        }

        blocked_csr_chunks[thread_id]->nnz_blocks = nnz_blocks_found;

        //Resize the arrays whose size is proportional to the number of non-zero blocks in the product chunk
        blocked_csr_chunks[thread_id]->block_col = (uint32_t*)realloc(blocked_csr_chunks[thread_id]->block_col,nnz_blocks_found*sizeof(uint32_t));
        if(blocked_csr_chunks[thread_id]->block_col == NULL){
            printf("Couldn't reallocate memory for block_col in blocked_bmm_parallel.\n");
            exit(-1);
        }

        blocked_csr_chunks[thread_id]->blocks = (comp_matrix**)realloc(blocked_csr_chunks[thread_id]->blocks,nnz_blocks_found*sizeof(comp_matrix*));
        if(blocked_csr_chunks[thread_id]->blocks == NULL){
            printf("Couldn't reallocate memory for blocks in blocked_bmm_parallel.\n");
            exit(-1);
        }

    }

    //Concatenate the blocked csr matrices of each chunk into one for the whole blocked matrix C
    block_comp_matrix* C = blocked_concat_chunks(blocked_csr_chunks,n_threads,A->n_b,A->real_dim);

    for(uint32_t i=0;i<n_threads;++i){
        if(blocked_csr_chunks[i]==NULL){
            continue;
        }
        free_block_comp_matrix(blocked_csr_chunks[i]);
    }
    free(blocked_csr_chunks);

    free(chunk_size);

    return C;
}