#include "utilities.h"
#include <string.h>


/**
 * Function that takes two matrices A and B in compressed format and produces
 * the product of those matrices.
 * Matrix A is in CSR format, whereas matrix B is in CSC format.
 * The function returns the product in CSR format.
 * We also use here an argument named concept. The reason is that we use this function as a subroutine
 * in blocked multiplication to multiply blocks so it is necessary to add the blocks' offset
 * so that the column indexes are correct in regards to the whole matrix and not only the 
 * block locally.
 **/
comp_matrix *bmm_seq(comp_matrix *A, comp_matrix *B, uint32_t offset) {
    //The size of C->col is equal to the columns of B.
    //However, we don't know the value of nnz.
    //We assume it to be roughly A->nnz + B->nnz.
    //If this is not enough the array C->row will be reallocated.
    comp_matrix *C = new_comp_matrix(A->nnz + B->nnz, B->n, "csr");

    //Number of non-zero elements in C up to this point
    uint32_t nnz_count = 0;

    //Pointers used to go through A->col and B->row respectively
    uint32_t A_row_ptr, A_row_end;
    uint32_t B_col_ptr, B_col_end;

    //For every row in A
    for (uint32_t i = 0; i < A->n; i++) {
        A_row_ptr = A->row[i];
        A_row_end = A->row[i + 1];

        C->row[i] = nnz_count;

        //For every column in B
        for (uint32_t j = 0; j < B->n; j++) {
            B_col_ptr = B->col[j];
            B_col_end = B->col[j + 1];

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
                //then we have a non-zero element in C in this position
                else{
                    C->col[nnz_count] = j + offset;
                    nnz_count++;
                    // Extend the C->col array if there are not enough empty cells
                    if (C->nnz == nnz_count) {
                        C->col = (uint32_t*)realloc(C->col, 2 * C->nnz * sizeof(uint32_t));
                        if(C->col == NULL){
                            printf("Couldn't reallocate memory for col in bmm_seq.\n");
                            exit(-1);
                        }
                        C->nnz = 2 * C->nnz;
                    }

                    //If we have already found a non-zero element in this position,
                    //there is no reason to continue calculating for this row of A and column of B.
                    break;
                }
            }

            // Reset the pointer of A->row to the previous value
            A_row_ptr = A->row[i];
        }
    }

    //In case the multiplication gives a matrix with zero elements only, retun a NULL matrix
    if(nnz_count == 0){
        free_comp_matrix(C);
        C = NULL;
        return C;
    }

    //Change to the true number of non zero values
    C->nnz = nnz_count;

    C->row[C->n] = nnz_count;

    //Resize col index array to correct size
    C->col = (uint32_t *)realloc(C->col, C->nnz * sizeof(uint32_t));

    return C;
}


/**
 * Following the same principle as the previous function,
 * this function takes two matrices A and B in compressed format and produces
 * the filtered product of those matrices according to a third matrix F.
 * Matrices A and F are in CSR format, whereas matrix B is in CSC format.
 * The function returns the product in CSR format.
 **/
comp_matrix *bmm_filtered_seq(comp_matrix *A, comp_matrix *B, comp_matrix *F, uint32_t offset) {
    // The max non zero values of matrix C is equal to the nnz of F
    comp_matrix *C = new_comp_matrix(F->nnz, F->n, "csr");

    // Used to fill C->row
    uint32_t nnz_count = 0;

    // Pointers used to go through A->col and B->row respectively
    uint32_t A_row_ptr, A_row_end;
    uint32_t B_col_ptr, B_col_end;

    // Iterate for every row of F
    for (int i = 0; i < F->n; i++) {

        C->row[i] = nnz_count;

        // Go through every non zero element of the i-th row
        for (int j = F->row[i]; j < F->row[i + 1]; j++) {
            A_row_ptr = A->row[i];
            A_row_end = A->row[i + 1];

            B_col_ptr = B->col[F->col[j]-offset];
            B_col_end = B->col[F->col[j]-offset+1];

            while (A_row_ptr < A_row_end && B_col_ptr < B_col_end) {
                if (A->col[A_row_ptr] < B->row[B_col_ptr]) {
                    A_row_ptr++;
                }
                else if (A->col[A_row_ptr] > B->row[B_col_ptr]) {
                    B_col_ptr++;
                } 
                else {
                    C->col[nnz_count] = F->col[j];
                    nnz_count++;
                    break;
                }
            }
        }
    }

    // Change to the true number of non zero values
    C->nnz = nnz_count;

    C->row[C->n] = nnz_count;

    // Reduce array to correct size
    C->col = (uint32_t *)realloc(C->col, C->nnz * sizeof(uint32_t));

    return C;
}


/**
 * Test function for small matrices.
 * Performs the traditional matrix multiplication where the matrices are stored in the
 * classic format.
**/
matrix_2d *bmm_seq_2d(matrix_2d *A, matrix_2d *B) {
    int rows = A->rows;
    int cols = B->cols;

    // Initialize the matrix that holds the product
    uint32_t **C = (uint32_t **)malloc(rows * sizeof(uint32_t *));
    for (uint32_t i = 0; i < rows; i++) {
        C[i] = (uint32_t *)malloc(cols * sizeof(uint32_t));
    }

    // Run every element of the matrix C
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < A->cols; k++) {
                if (A->mat[i][k] == 0 || B->mat[k][j] == 0) {
                    C[i][j] = 0;
                } else {
                    C[i][j] = 1;
                    break;
                }
            }
        }
    }

    matrix_2d *c_mat = (matrix_2d *)malloc(sizeof(matrix_2d));
    if(c_mat == NULL){
        printf("Couldn't allocate memory for c_mat in bmm_seq_2d.\n");
        exit(-1);
    }
    c_mat->cols = cols;
    c_mat->rows = rows;
    c_mat->mat = C;

    return c_mat;
}


/**
 * Function that applies element-wise union in two matrices A and B.
 * Both A and B are in csr format.
 * This function works as the merge part of a mergesort for each row of the matrix, sorting
 * elements of corresponding rows in A and B based on their column index. In case of non-zero 
 * elements with the same row and column indexes this element is only counted once.
**/
comp_matrix* block_union(comp_matrix* A, comp_matrix* B){
    uint32_t nnz_count = 0; //number of non-zero elements of C found up to this point
    uint32_t b = B->n;


    uint32_t A_row_end;     //index of the last element in a specific row of A
    uint32_t A_row_ptr;     //used to iterate through the elements of a row in A

    uint32_t B_row_end;     //index of the last element in a specific row of B
    uint32_t B_row_ptr;     //used to iterate through the elements of a row in B

    //If A is NULL and B is not NULL, then returned matrix is same as B
    if(A == NULL){
        A = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(A==NULL){
            printf("Couldn't allocate memory for A in block_union.\n");
            exit(-1);
        }

        A->col = (uint32_t*)malloc(B->nnz*sizeof(uint32_t));
        if(A->col==NULL){
            printf("Couldn't allocate memory for A->col in block_union.\n");
            exit(-1);
        }

        A->row = (uint32_t*)malloc((b+1)*sizeof(uint32_t));
        if(A->row==NULL){
            printf("Couldn't allocate memory for A->row in block_union.\n");
            exit(-1);
        }

        memcpy(A->col,B->col,B->nnz*sizeof(uint32_t));
        memcpy(A->row,B->row,(b+1)*sizeof(uint32_t));

        A->n = b;
        A->nnz = B->nnz;
        return A;
    }

    //Number of non-zero elements in returned matrix is less or equal than the sum
    //of non-zero elements in A and B.
    uint32_t* new_col = (uint32_t*)malloc((A->nnz+B->nnz)*sizeof(uint32_t));
    if(new_col == NULL){
        printf("Couldn't allocate memory for new_col in block_union_csr.\n");
        exit(-1);
    }

    uint32_t* new_row = (uint32_t*)calloc(b+1,sizeof(uint32_t));
    if(new_row == NULL){
        printf("Couldn't allocate memory for new_row in block_union_csr.\n");
        exit(-1);
    }

    //For every row of the two matrices
    for(uint32_t i=0;i<b;++i){
        A_row_ptr = A->row[i];
        A_row_end = A->row[i+1];

        B_row_ptr = B->row[i];
        B_row_end = B->row[i+1];

        //If there are no elements in this row for both matrices
        if((A_row_ptr==A_row_end) && (B_row_ptr==B_row_end)){
            new_row[i+1] = nnz_count;
            continue;
        }

        //While there are still non-zero elements in both rows
        while(A_row_ptr<A_row_end && B_row_ptr<B_row_end){
            //If A's current element's column idex is smaller than that of B's current element's
            if(A->col[A_row_ptr] < B->col[B_row_ptr]){
                new_col[nnz_count] = A->col[A_row_ptr];
                A_row_ptr++;
            }
            //If B's current element's column idex is smaller than that of A's current element's
            else if(B->col[B_row_ptr] < A->col[A_row_ptr]){
                new_col[nnz_count] = B->col[B_row_ptr];
                B_row_ptr++;
            }
            //If the column index of the two elements is the same
            else{
                new_col[nnz_count] = A->col[A_row_ptr];
                A_row_ptr++;
                B_row_ptr++;
            }
            nnz_count++;
        }

        //If there are non-zero elements left only in the row of A just add them in the same order they appear
        while(A_row_ptr < A_row_end){
            new_col[nnz_count] = A->col[A_row_ptr];
            A_row_ptr++;
            nnz_count++;
        }

        //Same if there are non-zero elements left only in the row of B
        while(B_row_ptr < B_row_end){
            new_col[nnz_count] = B->col[B_row_ptr];
            B_row_ptr++;
            nnz_count++;
        }

        new_row[i+1] = nnz_count;
    }

    //Reallocate column index array of the returned matrix to the correct size
    if(nnz_count < (A->nnz+B->nnz)){
        new_col = (uint32_t*)realloc(new_col,nnz_count*sizeof(uint32_t));
        if(new_col == NULL){
            printf("Couldn;t reallocate memory for new_col in block_union_csr.\n");
            exit(-1);
        }
    }

    free(A->col);
    free(A->row);

    A->col = new_col;
    A->row = new_row;
    A->nnz = nnz_count;

    return A;
}


/**
 * Function that mulitplies two matrices A and B.
 * A is in blocked csr format and B is in blocked csc format.
 * Matrix C = A*B is in csr blocked format.
 * This function works as the simple bmm function but instead of non-zero elements
 * we use non-zero blocks. Then if we have a match we multiply those blocks using
 * the simple bmm function.
**/
block_comp_matrix* blocked_bmm_seq(block_comp_matrix* A, block_comp_matrix* B){

    //Check if the corresponding dimensions are correct
    if(A->n_b != B->n_b){
        printf("Dimensions of the 2 matrixes not matching.\n");
        return NULL;
    }

    uint32_t n_b = A->n_b;

    //We assume A is a non-zero matrix
    uint32_t b = A->blocks[0]->n;

    block_comp_matrix* C = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if(C==NULL){
        printf("Couldn't allocate memory for C in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->block_row = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(C->block_row == NULL){
        printf("Couldn't allocate memory for block_row in blocked_bmm_seq.\n");
        exit(-1);
    }
     
    C->block_col = (uint32_t*)malloc(pow(n_b,2)*sizeof(uint32_t));
    if(C->block_col == NULL){
        printf("Couldn't allocate memory for block_col in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->blocks = (comp_matrix**)malloc(pow(n_b,2)*sizeof(comp_matrix*));
    if(C->blocks == NULL){
        printf("Couldn't allocate memory for blocks in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->n_b = n_b;
    C->real_dim = A->real_dim;

    uint32_t block_row_start;   //First non-zero block in this row of blocks
    uint32_t block_row_end;     //Last non-zero block in this row of blocks

    uint32_t block_col_start;   //First non-zero block in this column of blocks
    uint32_t block_col_end;     //Last non-zero block in this column of blocks

    uint32_t block_row_ptr;     //Iterator through the non-zero blocks of a row of blocks
    uint32_t block_col_ptr;     //Iterator through the non-zero blocks of a column of blocks

    uint32_t nnz_blocks_found = 0;  //Number of C's non-zero blocks found up to this point

    uint32_t first_match;       //Used as a flag to see whether it is the first time we have matching blocks to apply multiplication to

    comp_matrix* prod = NULL;

    //For each row of blocks in A
    for(uint32_t i=0;i<n_b;++i){
        block_row_start = A->block_row[i];
        block_row_end = A->block_row[i+1];

        //If there aren't any non-zero blocks in this row of blocks in A
        //then there aren't any non-zero blocks in this row of blocks in C either
        if(block_row_start == block_row_end){
            C->block_row[i+1] = C->block_row[i];
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

            //While there are still non-zero blocks in the row of blocks in A or the column of blocks in B
            while(block_row_ptr<block_row_end && block_col_ptr<block_col_end){
                //Search for blocks in A with the same col index as the row index of the blocks in B
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
                        C->blocks[nnz_blocks_found] = NULL;
                        C->block_col[nnz_blocks_found] = j;
                        first_match = 1;
                        nnz_blocks_found++; 

                    }

                    //No reason to apply union if the product is NULL
                    if(prod != NULL){
                        C->blocks[nnz_blocks_found-1] = block_union(C->blocks[nnz_blocks_found-1],prod);
                        free_comp_matrix(prod); 
                        prod = NULL;
                    }

                    block_col_ptr++;
                    block_row_ptr++;    
          
                }
            }
        }
        C->block_row[i+1] = nnz_blocks_found;
    }

    //If we don't have any non-zero blocks then return a NULL matrix
    if(nnz_blocks_found==0){
        free_block_comp_matrix(C);
        C = NULL;
        return C;
    }

    C->nnz_blocks = nnz_blocks_found;

    //Resize the arrays whose size is proportional to the number of non-zero blocks in C
    C->block_col = (uint32_t*)realloc(C->block_col,nnz_blocks_found*sizeof(uint32_t));
    if(C->block_col == NULL){
        printf("Couldn't reallocate memory for block_col in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->blocks = (comp_matrix**)realloc(C->blocks,nnz_blocks_found*sizeof(comp_matrix*));
    if(C->blocks == NULL){
        printf("Couldn't reallocate memory for blocks in blocked_bmm_seq.\n");
        exit(-1);
    }

    return C;
}


block_comp_matrix *blocked_bmm_seq_filtered(block_comp_matrix *A, block_comp_matrix *B, block_comp_matrix *F){
    
    //Check if the corresponding dimensions are correct
    if(A->n_b != B->n_b || A->n_b != F->n_b){
        printf("Dimensions of the matrices not matching.\n");
        return NULL;
    }

    uint32_t n_b = F->n_b;

    //We assume F is a non-zero matrix
    uint32_t b = F->blocks[0]->n;

    block_comp_matrix* C = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if(C==NULL){
        printf("Couldn't allocate memory for C in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->block_row = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(C->block_row == NULL){
        printf("Couldn't allocate memory for block_row in blocked_bmm_seq.\n");
        exit(-1);
    }
     
    C->block_col = (uint32_t*)malloc(pow(n_b,2)*sizeof(uint32_t));
    if(C->block_col == NULL){
        printf("Couldn't allocate memory for block_col in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->blocks = (comp_matrix**)malloc(pow(n_b,2)*sizeof(comp_matrix*));
    if(C->blocks == NULL){
        printf("Couldn't allocate memory for blocks in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->n_b = n_b;
    C->real_dim = F->real_dim;

    uint32_t F_block_row_start; // First non-zero block in this row of blocks in matrix F
    uint32_t F_block_row_end;   // Last non-zero block in this row of blocks in matrix F

    uint32_t block_row_start;   //First non-zero block in this row of blocks
    uint32_t block_row_end;     //Last non-zero block in this row of blocks

    uint32_t block_col_start;   //First non-zero block in this column of blocks
    uint32_t block_col_end;     //Last non-zero block in this column of blocks

    uint32_t block_row_ptr;     //Iterator through the non-zero blocks of a row of blocks
    uint32_t block_col_ptr;     //Iterator through the non-zero blocks of a column of blocks

    uint32_t nnz_blocks_found = 0;  //Number of C's non-zero blocks found up to this point

    uint32_t first_match;       //Used as a flag to see whether it is the first time we have matching blocks to apply multiplication to

    comp_matrix* prod = NULL;

    //For each row of blocks in F
    for(uint32_t i=0;i<n_b;++i){
        F_block_row_start = F->block_row[i];
        F_block_row_end = F->block_row[i+1];

        block_row_start = A->block_row[i];
        block_row_end = A->block_row[i+1];

        //If there aren't any non-zero blocks in this row of blocks in F or A
        //then there aren't any non-zero blocks in this row of blocks in C either
        if(F_block_row_start == F_block_row_end || block_row_start == block_row_end){
            C->block_row[i+1] = C->block_row[i];
            continue;
        }

        //For each non zero block of F in the i-th block row
        for(uint32_t j=F_block_row_start;j<F_block_row_end;++j){

            block_col_start = B->block_col[F->block_col[j]];
            block_col_end = B->block_col[F->block_col[j]+1];

            //If there aren't any non-zero blocks in this column go to the next one
            if(block_col_start == block_col_end){
                continue;
            }
            
            block_row_ptr = block_row_start;
            block_col_ptr = block_col_start;

            first_match = 0;

            //While there are still non-zero blocks in the row of blocks in A or the column of blocks in B
            while(block_row_ptr<block_row_end && block_col_ptr<block_col_end){
                //Search for blocks in A with the same col index as the row index of the blocks in B
                if(A->block_col[block_row_ptr] > B->block_row[block_col_ptr]){
                    block_col_ptr++;
                }
                else if(A->block_col[block_row_ptr] < B->block_row[block_col_ptr]){
                    block_row_ptr++;
                }
                else{
                    //In this part we multiply the corresponding blocks
                    //prod is a csr matrix
                    prod = bmm_filtered_seq(A->blocks[block_row_ptr],B->blocks[block_col_ptr],F->blocks[j],F->block_col[j]*b);

                    //If it is the first product we calculate for this block of the new blocked matrix
                    if(first_match==0){

                        //If the multiplication doesn't yield a non-zero matrix then do nothing
                        if(prod == NULL){
                            block_col_ptr++;
                            block_row_ptr++;
                            continue;
                        }

                        //Else create a new non-zero block for the new matrix
                        C->blocks[nnz_blocks_found] = NULL;
                        C->block_col[nnz_blocks_found] = j;
                        first_match = 1;
                        nnz_blocks_found++; 

                    }

                    //No reason to apply union if the product is NULL
                    if(prod != NULL){
                        C->blocks[nnz_blocks_found-1] = block_union(C->blocks[nnz_blocks_found-1],prod);
                        free_comp_matrix(prod); 
                        prod = NULL;
                    }

                    block_col_ptr++;
                    block_row_ptr++;    
          
                }
            }
        }
        C->block_row[i+1] = nnz_blocks_found;
    }

    //If we don't have any non-zero blocks then return a NULL matrix
    if(nnz_blocks_found==0){
        free_block_comp_matrix(C);
        C = NULL;
        return C;
    }

    C->nnz_blocks = nnz_blocks_found;

    //Resize the arrays whose size is proportional to the number of non-zero blocks in C
    C->block_col = (uint32_t*)realloc(C->block_col,nnz_blocks_found*sizeof(uint32_t));
    if(C->block_col == NULL){
        printf("Couldn't reallocate memory for block_col in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->blocks = (comp_matrix**)realloc(C->blocks,nnz_blocks_found*sizeof(comp_matrix*));
    if(C->blocks == NULL){
        printf("Couldn't reallocate memory for blocks in blocked_bmm_seq.\n");
        exit(-1);
    }

    return C;
}