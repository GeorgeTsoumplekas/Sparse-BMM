#include "utilities.h"

/**
 * Function that takes two matrices A and B in compressed format and produces
 * the product of those matrices.
 * Matrix A is in CSR format, whereas matrix B is in CSC format.
 * The function returns the product in CSC format.
 **/
comp_matrix *bmm_seq(comp_matrix *A, comp_matrix *B) {
    // The size of C->col is correct and equal to the columns of B.
    // However, we don't know the value of nnz.
    // We assume it at roughly A->nnz + B->nnz.
    // If this is not enough the array C->row will be reallocated.
    comp_matrix *C = new_comp_matrix(A->nnz + B->nnz, B->n, "csc");

    // Used to fill C->row and C->col
    uint32_t col_indx = 0;

    // Pointers used to go through A->col and B->row respectively
    int ptr_in_col_A, end_in_col_A;
    int ptr_in_row_B, end_in_row_B;

    // Iterate for every col of B
    for (int i = 0; i < B->n; i++) {
        // Find the non zero elements position of the i-th column in the array
        // B->row
        ptr_in_row_B = B->col[i];
        end_in_row_B = B->col[i + 1];

        C->col[i] = col_indx;

        // And iterate for every row of A
        for (int j = 0; j < A->n; j++) {
            ptr_in_col_A = A->row[j];
            end_in_col_A = A->row[j + 1];

            while (ptr_in_col_A < end_in_col_A && ptr_in_row_B < end_in_row_B) {
                if (A->col[ptr_in_col_A] < B->row[ptr_in_row_B]) {
                    ptr_in_col_A++;
                } else if (A->col[ptr_in_col_A] > B->row[ptr_in_row_B]) {
                    ptr_in_row_B++;
                } else if (A->col[ptr_in_col_A] == B->row[ptr_in_row_B]) {
                    C->row[col_indx] = j;
                    col_indx++;
                    // Extend the C->row array if there are not enough empty
                    // cells
                    if (C->nnz == col_indx) {
                        C->row = realloc(C->row, 2 * C->nnz * sizeof(uint32_t));
                        C->nnz = 2 * C->nnz;
                    }

                    break;
                }
            }

            // Reset the pointer of B->row to the previous value
            ptr_in_row_B = B->col[i];
        }
    }

    // Change to the true number of non zero values
    C->nnz = col_indx;

    C->col[C->n] = col_indx;

    // Reduce array to correct size
    C->row = (uint32_t *)realloc(C->row, C->nnz * sizeof(uint32_t));

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

/**
 * Function that mulitplies two matrices A and B.
 * A is in blocked csr format and B is in blocked csc format
 * Matrix C = A*B is in csr blocked format.
**/
block_comp_matrix* blocked_bmm_seq(block_comp_matrix* A, block_comp_matrix* B){

    block_comp_matrix* C = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if(C==NULL){
        printf("Couldn't allocate memory for C in blocked_bmm_seq.\n");
        exit(-1);
    }

    C->nb = A->nb;
    C->b = A->b;

    uint32_t block_row_start;   //First non-zero block in this row of blocks
    uint32_t block_row_end;     //Last non-zero block in this row of blocks

    uint32_t block_col_start;   //First non-zero block in this column of blocks
    uint32_t block_col_end;     //Last non-zero block in this column of blocks

    uint32_t block_row_ptr;     //Iterator through the non-zero blocks of a row of blocks
    uint32_t block_col_ptr;     //Iterator through the non-zero blocks of a column of blocks

    //for each row of blocks in A
    for(uint32_t i=0;i<A->nb;++i){
        block_row_start = A->line_blocks[i];
        block_row_end = A->line_blocks[i+1];

        //for each column of blocks in B
        for(uint32_t j=0;j<B->nb;++j){
            block_col_start = B->line_blocks[j];
            block_col_end = B->line_blocks[j+1];

            block_row_ptr = block_row_start;
            block_col_ptr = block_col_start;

            //while there are still non-zero blocks in the row of blocks in A or the column of blocks in B
            while(block_row_ptr<block_row_end && block_col_ptr<block_col_end){
                //search for blocks in A with the same col index as the row index of the blocks in B
                if(block_row_ptr>block_col_ptr){
                    block_col_ptr++;
                }
                else if(block_row_ptr<block_col_ptr){
                    block_row_ptr++;
                }
                else{
                    //TODO: pollaplasiasmos twn block pinakwn gia tous opoious exoume antistoixish//
                    block_col_ptr++;
                    block_row_ptr++;                
                }
            }
        }
    }
    //TODO: elegxos orthotitas molis oloklhrwthei h blocked_bmm_seq
    //TODO: blocked_bmm_seq_filtered
    //TODO: sunartish pou metatrepei apo blocked se mh blocked morfh


    return C;
}