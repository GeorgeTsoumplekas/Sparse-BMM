#include <omp.h>

#include "readmtx.h"
#include "utilities.h"

/*********   COO format functions   *********/
typedef struct {
    uint32_t *row;
    uint32_t *col;
    uint32_t nnz;
} coo;

coo *new_coo(uint32_t nnz) {
    coo *mat = (coo *)malloc(sizeof(coo));
    mat->col = (uint32_t *)calloc(nnz, sizeof(uint32_t));
    mat->row = (uint32_t *)calloc(nnz, sizeof(uint32_t));
    mat->nnz = nnz;
}

void free_coo(coo *mat) {
    free(mat->col);
    free(mat->row);
    free(mat);
}

void print_coo(coo *array) {
    printf("Rows: ");
    for (int i = 0; i < array->nnz; i++) {
        printf("%d ", array->row[i]);
    }
    printf("\nCols: ");
    for (int i = 0; i < array->nnz; i++) {
        printf("%d ", array->col[i]);
    }
    printf("\n");
}

// Function that takes two compressed matrices, one in csr format (A)
// and the other in csc format (B), finds the boolean multiplication
// on a specific point with coordinates (row_indx, col_indx) and stores
// the result in a matrix in coo format (C).
// The return value is 1 if a match was found, otherwise it is 0.
int mult_row_with_col(comp_matrix A, comp_matrix B, coo *C, int row_indx, int col_indx, int pos) {
    // Check if point coordinates are valid
    if (row_indx >= A.n || col_indx >= B.n) {
        return 0;
    }

    int ptr_in_col_A = A.row[row_indx];
    int end_in_col_A = A.row[row_indx + 1];

    int ptr_in_row_B = B.col[col_indx];
    int end_in_row_B = B.col[col_indx + 1];

    while (ptr_in_col_A < end_in_col_A && ptr_in_row_B < end_in_row_B) {
        if (A.col[ptr_in_col_A] < B.row[ptr_in_row_B]) {
            ptr_in_col_A++;
        } else if (A.col[ptr_in_col_A] > B.row[ptr_in_row_B]) {
            ptr_in_row_B++;
        } else if (A.col[ptr_in_col_A] == B.row[ptr_in_row_B]) {
            // Store the point (row_indx + 1, col_index + 1)
            // This becomes more clear in the following functions
            C->col[pos] = col_indx + 1;
            C->row[pos] = row_indx + 1;
            return 1;
        }
    }

    return 0;
}

comp_matrix *bmm_parallel(comp_matrix *A, comp_matrix *B) {
    // We devide the product matrix C into chunks to control the memory needed
    uint32_t chunk_size_rows = A->n;
    uint32_t chunk_size_cols = B->n;

    uint64_t chunk_size = chunk_size_cols * chunk_size_rows;

    // Don't let chunk_size surpass the size of 10^6
    if (chunk_size > 1000000) {
        chunk_size_rows = 1000000 / chunk_size_cols;
        chunk_size = chunk_size_cols * chunk_size_rows;
    }

    // We don't know the value of nnz and
    // we assume it at roughly A->nnz + B->nnz.
    // If this is not enough the arrays C->row and C->col will be reallocated.
    coo *C_reduced = new_coo(A->nnz + B->nnz);

    // Temporary matrix that holds the some of the non zero coordinates
    coo *C_temp = new_coo(chunk_size);

    int nnz = 0;
    int total_nnz = 0;

    // Find the iterations of the outer loop
    // The multiplication of the matrices happens gradually depending on the chunks
    int iterations = A->n / chunk_size_rows;

    int last_element_filled = 0;

    for (int i = 0; i < iterations; i++) {
        #pragma omp parallel for schedule(dynamic) reduction(+ : nnz) collapse(2)
        // Iterate every element in the chunk of rows
        for (int j = 0; j < chunk_size_cols; j++) {
            for (int k = 0; k < chunk_size_rows; k++) {
                nnz += mult_row_with_col(*A, *B, C_temp, i * chunk_size_rows + k, j, j * chunk_size_rows + k);
            }
        }

        // Holds the total number of non-zero elements
        total_nnz += nnz;
        // If the size of C_reduced isn't enough reallocate memory
        if (total_nnz > C_reduced->nnz) {
            C_reduced->col = realloc(C_reduced->col, 2 * C_reduced->nnz * sizeof(uint32_t));
            C_reduced->row = realloc(C_reduced->row, 2 * C_reduced->nnz * sizeof(uint32_t));
            C_reduced->nnz = 2 * C_reduced->nnz;
        }

        // Copy contents to the reduced C arrays and reset C_temp
        for (int j = 0; j < chunk_size; j++) {
            if (C_temp->col[j] != 0 && C_temp->row[j] != 0) {
                C_reduced->col[last_element_filled] = C_temp->col[j] - 1;
                C_reduced->row[last_element_filled] = C_temp->row[j] - 1;
                last_element_filled++;
            }
            C_temp->col[j] = 0;
            C_temp->row[j] = 0;
        }

        // Reset variable
        nnz = 0;
    }

    free_coo(C_temp);

    // Fix values of C_reduced
    C_reduced->col = realloc(C_reduced->col, total_nnz * sizeof(uint32_t));
    C_reduced->row = realloc(C_reduced->row, total_nnz * sizeof(uint32_t));
    C_reduced->nnz = total_nnz;

    // Transform the coo format of C to csc format
    comp_matrix *c_csc = new_comp_matrix(C_reduced->nnz, B->n, "csc");
    coo2csc(c_csc, C_reduced->row, C_reduced->col);

    free_coo(C_reduced);

    return c_csc;
}

comp_matrix *bmm_filtered_parallel(comp_matrix *A, comp_matrix *B, comp_matrix *F) {
    coo *C = new_coo(F->nnz);

    // Holds the number of the non zero elements of the product
    int nnz = 0;

    #pragma omp parallel for schedule(dynamic) reduction(+ : nnz)
    // Iterate for every column of F
    for (int i = 0; i < F->n; i++) {
        // Go through every non zero element of the i-th column
        for (int j = F->col[i]; j < F->col[i + 1]; j++) {
            nnz += mult_row_with_col(*A, *B, C, F->row[j], i, j);
        }
    }

    // Matrix that holds the correct size of C in coo format
    coo *C_reduced = new_coo(nnz);

    int last_element_filled = 0;
    // Copy the contents to the reduced matrix
    for (int i = 0; i < C->nnz; i++) {
        if (C->col[i] != 0 && C->row[i] != 0) {
            C_reduced->col[last_element_filled] = C->col[i] - 1;
            C_reduced->row[last_element_filled] = C->row[i] - 1;
            last_element_filled++;
        }
    }
    free_coo(C);

    // Transform the coo format of C to csc format
    comp_matrix *c_csc = new_comp_matrix(C_reduced->nnz, B->n, "csc");
    coo2csc(c_csc, C_reduced->row, C_reduced->col);

    free_coo(C_reduced);

    return c_csc;
}
