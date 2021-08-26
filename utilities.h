#ifndef UTILITIES_H
#define UTILITIES_H

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ......   Structs   ......*/
// This struct can either describe a Compressed Sparse Column (CSC) or a
// Compressed Sparse Row (CSR) matrix format.
// If it is CSC then the col array has n+1 values and the row array has nnz
// values.
// If it is CSR then the col array has nnz values and the row array has n+1
// values.
typedef struct {
    uint32_t* col;
    uint32_t* row;
    uint32_t nnz;
    uint32_t n;
} comp_matrix;

typedef struct {
    int** mat;
    int rows;
    int cols;
} matrix_2d;

/* ......   Functions   ......*/
// Functions for the comp_matrix struct
comp_matrix* new_comp_matrix(uint32_t nnz, uint32_t n, char* type) {
    comp_matrix* array = (comp_matrix*)malloc(sizeof(comp_matrix));
    if (array == NULL) exit(-1);

    if (!strcmp(type, "csc")) {
        array->col = (uint32_t*)malloc((n + 1) * sizeof(uint32_t));
        array->row = (uint32_t*)malloc(nnz * sizeof(uint32_t));
    } else if (!strcmp(type, "csr")) {
        array->col = (uint32_t*)malloc(nnz * sizeof(uint32_t));
        array->row = (uint32_t*)malloc((n + 1) * sizeof(uint32_t));
    } else {
        printf("Unknown type in new_comp_matrix function.\n");
    }

    if (array->col == NULL || array->row == NULL) {
        free(array);
        exit(-1);
    }

    array->nnz = nnz;
    array->n = n;
}

void free_compressed(comp_matrix* array) {
    free(array->col);
    free(array->row);
    free(array);
}

/* Functions for the matrix_2d struct */
// Creates a random 2d matrix
matrix_2d* rand_matrix_2d(int rows, int cols) {
    matrix_2d* array = (matrix_2d*)malloc(sizeof(matrix_2d));
    if (array == NULL) exit(-1);

    int** matrix = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (int*)malloc(cols * sizeof(int));
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand() % 2;
        }
    }

    array->mat = matrix;
    array->cols = cols;
    array->rows = rows;
}

void free_matrix_2d(matrix_2d* array) {
    for (int i = 0; i < array->rows; i++) {
        free(array->mat[i]);
    }
    free(array->mat);
    free(array);
}

// Prints an array of type matrix_2d
void print_matrix_2d(matrix_2d* array) {
    for (int i = 0; i < array->rows; i++) {
        for (int j = 0; j < array->cols; j++) {
            printf("%d ", array->mat[i][j]);
        }
        printf("\n");
    }
}

/* General functions */
// Function that turns a matrix_2d array to csc format
comp_matrix* matrix2csc(matrix_2d* array) {
    uint32_t nnz = 0;
    // Find the number of the non zero elements
    for (int i = 0; i < array->rows; i++) {
        for (int j = 0; j < array->cols; j++) {
            if (array->mat[i][j] != 0) nnz++;
        }
    }

    comp_matrix* csc = new_comp_matrix(nnz, (uint32_t)array->cols, "csc");

    int col_array_idx = 1;
    int row_array_idx = 0;

    csc->col[0] = 0;

    for (int j = 0; j < array->cols; j++) {
        for (int i = 0; i < array->rows; i++) {
            if (array->mat[i][j] != 0) {
                csc->row[row_array_idx] = i;
                row_array_idx++;
            }
        }
        csc->col[col_array_idx] = row_array_idx;
        col_array_idx++;
    }

    return csc;
}

// Function that turns a matrix_2d array to csr format
comp_matrix* matrix2csr(matrix_2d* array) {
    uint32_t nnz = 0;
    // Find the number of the non zero elements
    for (int i = 0; i < array->rows; i++) {
        for (int j = 0; j < array->cols; j++) {
            if (array->mat[i][j] != 0) nnz++;
        }
    }

    comp_matrix* csr = new_comp_matrix(nnz, (uint32_t)array->cols, "csr");

    int col_array_idx = 0;
    int row_array_idx = 1;

    csr->row[0] = 0;

    for (int i = 0; i < array->rows; i++) {
        for (int j = 0; j < array->cols; j++) {
            if (array->mat[i][j] != 0) {
                csr->col[col_array_idx] = j;
                col_array_idx++;
            }
        }
        csr->row[row_array_idx] = col_array_idx;
        row_array_idx++;
    }

    return csr;
}

// Prints a csc matrix
void print_csc(comp_matrix* array) {
    printf("Printing a csc matrix.\n");
    for (int i = 0; i < array->n + 1; i++) {
        printf("%d ", array->col[i]);
    }
    printf("\n");
    for (int i = 0; i < array->nnz; i++) {
        printf("%d ", array->row[i]);
    }
    printf("\n");
}

// Prints a csr matrix
void print_csr(comp_matrix* array) {
    printf("Printing a csr matrix.\n");
    for (int i = 0; i < array->n + 1; i++) {
        printf("%d ", array->row[i]);
    }
    printf("\n");
    for (int i = 0; i < array->nnz; i++) {
        printf("%d ", array->col[i]);
    }
    printf("\n");
}

#endif