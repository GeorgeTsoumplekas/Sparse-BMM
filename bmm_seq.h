#include "utilities.h"
#include <string.h>

/**
 * Function that takes two matrices A and B in compressed format and produces
 * the product of those matrices.
 * Matrix A is in CSR format, whereas matrix B is in CSC format.
 * The function returns the product in CSR format.
 **/
comp_matrix *bmm_seq(comp_matrix *A, comp_matrix *B, uint32_t offset) {
    // The size of C->col is correct and equal to the columns of B.
    // However, we don't know the value of nnz.
    // We assume it at roughly A->nnz + B->nnz.
    // If this is not enough the array C->row will be reallocated.
    comp_matrix *C = new_comp_matrix(A->nnz + B->nnz, B->n, "csr");

    // Number of non-zero elements in C up to this point
    uint32_t nnz_count = 0;

    // Pointers used to go through A->col and B->row respectively
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

            while (A_row_ptr < A_row_end && B_col_ptr < B_col_end) {
                if (A->col[A_row_ptr] < B->row[B_col_ptr]) {
                    A_row_ptr++;
                }
                else if (A->col[A_row_ptr] > B->row[B_col_ptr]) {
                    B_col_ptr++;
                }
                else{
                    C->col[nnz_count] = j + offset;
                    nnz_count++;
                    // Extend the C->col array if there are not enough empty
                    // cells
                    if (C->nnz == nnz_count) {
                        C->col = realloc(C->col, 2 * C->nnz * sizeof(uint32_t));
                        C->nnz = 2 * C->nnz;
                    }

                    break;
                }
            }

            // Reset the pointer of A->row to the previous value
            A_row_ptr = A->row[i];
        }
    }

    //In case the multiplication gives a matrix with zero elements only
    if(nnz_count == 0){
        free(C->col);
        free(C->row);
        C = NULL;
        return C;
    }

    // Change to the true number of non zero values
    C->nnz = nnz_count;

    C->row[C->n] = nnz_count;

    // Reduce array to correct size
    C->col = (uint32_t *)realloc(C->col, C->nnz * sizeof(uint32_t));

    return C;
}

/**
 * Function that takes two matrices A and B in compressed format and produces
 * the filtered product of those matrices according to a third matrix F.
 * Matrix A is in CSR format, whereas matrices B and F are in CSC format.
 * The function returns the product in CSC format.
 **/
comp_matrix *bmm_filtered_seq(comp_matrix *A, comp_matrix *B, comp_matrix *F) {
    // TODO: check if line is entirely correct
    comp_matrix *C = new_comp_matrix(F->nnz, F->n, "csc");

    // Used to fill C->row and C->col
    uint32_t row_indx = 0;
    uint32_t col_indx = 0;

    // Used to go through a col in matrix F
    int start, end;

    // Pointers used to go through A->col and B->row respectively
    int ptr_in_col_A, end_in_col_A;
    int ptr_in_row_B, end_in_row_B;

    // Iterate for every column of F
    for (int i = 0; i < F->n; i++) {
        // Find the position to start and end the iteration in F->row
        start = F->col[i];
        end = F->col[i + 1];

        C->col[i] = col_indx;

        // Go through every non zero element of the i-th column
        for (int j = start; j < end; j++) {
            // Find the row index
            row_indx = F->row[j];

            ptr_in_col_A = A->row[row_indx];
            end_in_col_A = A->row[row_indx + 1];

            ptr_in_row_B = B->col[i];
            end_in_row_B = B->col[i + 1];

            while (ptr_in_col_A < end_in_col_A && ptr_in_row_B < end_in_row_B) {
                if (A->col[ptr_in_col_A] < B->row[ptr_in_row_B]) {
                    ptr_in_col_A++;
                } else if (A->col[ptr_in_col_A] > B->row[ptr_in_row_B]) {
                    ptr_in_row_B++;
                } else if (A->col[ptr_in_col_A] == B->row[ptr_in_row_B]) {
                    C->row[col_indx] = row_indx;
                    col_indx++;
                    break;
                }
            }
        }
    }

    // Change to the true number of non zero values
    C->nnz = col_indx;

    C->col[C->n] = col_indx;

    // Reduce array to correct size
    C->row = (uint32_t *)realloc(C->row, C->nnz * sizeof(uint32_t));

    return C;
}

// Test function for small matrices
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


//ousiastika tha leitourgei san merge-sort se kathe seira
//oi A kai B einai se csr morfh
comp_matrix* block_union(comp_matrix* A, comp_matrix* B){
    uint32_t nnz_count = 0;
    uint32_t b = B->n;

    uint32_t A_row_end;
    uint32_t A_row_ptr;

    uint32_t B_row_end;
    uint32_t B_row_ptr;

    //If A is NULL and B is not NULL
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

        while(A_row_ptr<A_row_end && B_row_ptr<B_row_end){
            if(A->col[A_row_ptr] < B->col[B_row_ptr]){
                new_col[nnz_count] = A->col[A_row_ptr];
                A_row_ptr++;
            }
            else if(B->col[B_row_ptr] < A->col[A_row_ptr]){
                new_col[nnz_count] = B->col[B_row_ptr];
                B_row_ptr++;
            }
            else{
                new_col[nnz_count] = A->col[A_row_ptr];
                A_row_ptr++;
                B_row_ptr++;
            }
            nnz_count++;
        }

        while(A_row_ptr < A_row_end){
            new_col[nnz_count] = A->col[A_row_ptr];
            A_row_ptr++;
            nnz_count++;
        }

        while(B_row_ptr < B_row_end){
            new_col[nnz_count] = B->col[B_row_ptr];
            B_row_ptr++;
            nnz_count++;
        }

        new_row[i+1] = nnz_count;
    }

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
 * A is in blocked csr format and B is in blocked csc format
 * Matrix C = A*B is in csr blocked format.
**/
block_comp_matrix_2* blocked_bmm_seq(block_comp_matrix_2* A, block_comp_matrix_2* B){

    //Check if the corresponding dimensions are correct
    if(A->n_b != B->n_b){
        printf("Dimensions of the 2 matrixes not matching.\n");
        return NULL;
    }

    uint32_t n_b = A->n_b;

    //We assume A is a non-zero matrix
    uint32_t b = A->blocks[0]->n;

    block_comp_matrix_2* C = (block_comp_matrix_2*)malloc(sizeof(block_comp_matrix_2));
    if(C==NULL){
        printf("Couldn't allocate memory for C in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->block_row = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(C->block_row == NULL){
        printf("Couldn;t allocate memory for block_row in blocked_bmm_seq.\n");
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

    uint32_t block_row_start;   //First non-zero block in this row of blocks
    uint32_t block_row_end;     //Last non-zero block in this row of blocks

    uint32_t block_col_start;   //First non-zero block in this column of blocks
    uint32_t block_col_end;     //Last non-zero block in this column of blocks

    uint32_t block_row_ptr;     //Iterator through the non-zero blocks of a row of blocks
    uint32_t block_col_ptr;     //Iterator through the non-zero blocks of a column of blocks

    uint32_t nnz_blocks_found = 0;  //number of C's non-zero blocks found up to this point

    uint32_t first_match;       //Used as a flag to see whether it is the first time we have matching blocks to apply multiplication to

    comp_matrix* prod = NULL;

    //For each row of blocks in A
    for(uint32_t i=0;i<n_b;++i){
        block_row_start = A->block_row[i];
        block_row_end = A->block_row[i+1];

        if(block_row_start == block_row_end){
            C->block_row[i+1] = C->block_row[i];
            continue;
        }

        //For each column of blocks in B
        for(uint32_t j=0;j<n_b;++j){
            block_col_start = B->block_col[j];
            block_col_end = B->block_col[j+1];
            
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

                    if(prod != NULL){
                        C->blocks[nnz_blocks_found-1] = block_union(C->blocks[nnz_blocks_found-1],prod);
                        free(prod); 
                        prod = NULL;
                    }

                    block_col_ptr++;
                    block_row_ptr++;    
          
                }
            }
        }
        C->block_row[i+1] = nnz_blocks_found;
    }

    C->nnz_blocks = nnz_blocks_found;

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

    //TODO: elegxos orthotitas molis oloklhrwthei h blocked_bmm_seq
    //TODO: blocked_bmm_seq_filtered