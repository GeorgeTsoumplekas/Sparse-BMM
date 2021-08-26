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
    int **C = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        C[i] = (int *)malloc(cols * sizeof(int));
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
    c_mat->cols = cols;
    c_mat->rows = rows;
    c_mat->mat = C;

    return c_mat;
}