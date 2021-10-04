#ifndef READMTX_H
#define READMTX_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mmio.h"
#include "utilities.h"

// Compare function used in qsort function
int cmpfunc(const void* a, const void* b) { return (*(int*)a - *(int*)b); }

/**
 *  \brief COO to CSC conversion
 *
 *  Converts a square matrix from COO to CSC format.
 *
 *  Note: The routine assumes the input COO and the output CSC matrix
 *  to be square.
 *
 */
void coo_to_csc(comp_matrix* csc, coo_matrix* coo) {
    // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < csc->n + 1; l++) {
        csc->col[l] = 0;
    }

    // ----- find the correct column sizes
    for (uint32_t l = 0; l < csc->nnz; l++) {
        csc->col[coo->col[l]]++;
    }

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < csc->n; i++) {
        uint32_t temp = csc->col[i];
        csc->col[i] = cumsum;
        cumsum += temp;
    }
    csc->col[csc->n] = csc->nnz;

    // ----- copy the row indices to the correct place
    for (uint32_t l = 0; l < csc->nnz; l++) {
        uint32_t col_l;
        col_l = coo->col[l];

        uint32_t dst = csc->col[col_l];
        csc->row[dst] = coo->row[l];

        csc->col[col_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < csc->n; i++) {
        uint32_t temp = csc->col[i];
        csc->col[i] = last;
        last = temp;
    }

    // ----- sort the elements of every col in csc->row
    for (uint32_t i = 0; i < csc->n; i++) {
        qsort(csc->row + csc->col[i], csc->col[i + 1] - csc->col[i],
              sizeof(int), cmpfunc);
    }
}

/**
 *  \brief COO to CSR conversion
 *
 *  Converts a square matrix from COO to CSR format.
 *
 *  Note: The routine assumes the input COO and the output CSR matrix
 *  to be square.
 *
 */
void coo_to_csr(comp_matrix* csr, coo_matrix* coo) {
    // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < csr->n + 1; l++) {
        csr->row[l] = 0;
    }

    // ----- find the correct row sizes
    for (uint32_t l = 0; l < csr->nnz; l++) {
        csr->row[coo->row[l]]++;
    }

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < csr->n; i++) {
        uint32_t temp = csr->row[i];
        csr->row[i] = cumsum;
        cumsum += temp;
    }
    csr->row[csr->n] = csr->nnz;

    // ----- copy the column indices to the correct place
    for (uint32_t l = 0; l < csr->nnz; l++) {
        uint32_t row_l;
        row_l = coo->row[l];

        uint32_t dst = csr->row[row_l];
        csr->col[dst] = coo->col[l];

        csr->row[row_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < csr->n; i++) {
        uint32_t temp = csr->row[i];
        csr->row[i] = last;
        last = temp;
    }

    // -----sort the elements of every row in csr->col
    for (uint32_t i = 0; i < csr->n; i++) {
        qsort(csr->col + csr->row[i], csr->row[i + 1] - csr->row[i],
              sizeof(int), cmpfunc);
    }
}

// Takes as input a mtx file and according to the type produces the appropriate
// compressed matrix.
// If type == "csc" it returns the matrix in csc format and if type == "csr" it
// returns the matrix in csr format.
comp_matrix* mtx2comp(char* filename, char* type) {
    int ret_code;
    MM_typecode matcode;
    FILE* f;
    int M, N, nz;
    int i;
    uint32_t *coo_row, *coo_col;
    double* val;

    if ((f = fopen(filename, "r")) == NULL) {
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0) {
        exit(1);
    }

    /* reserve memory for matrices */
    // We assume the matrices are symmetrical and allocate double the size
    coo_row = (uint32_t*)malloc(2 * nz * sizeof(uint32_t));
    coo_col = (uint32_t*)malloc(2 * nz * sizeof(uint32_t));

    if (coo_row == NULL){
        printf("Couldn't allocate memory for coo_row in mtx2comp.\n");
        exit(-1);
    }

    if (coo_col == NULL){
        printf("Couldn't allocate memory for coo_col in mtx2comp.\n");
        exit(-1);
    }

    // Acquire the points in the file
    for (i = 0; i < nz; i++) {
        int ign = fscanf(f, "%d %d\n", &coo_row[i], &coo_col[i]);
        coo_row[i]--; /* adjust from 1-based to 0-based */
        coo_col[i]--;

        // Add the symmetrical elements
        coo_row[i + nz] = coo_col[i];
        coo_col[i + nz] = coo_row[i];
    }

    if (f != stdin) {
        fclose(f);
    }

    coo_matrix* coo = (coo_matrix*)malloc(sizeof(coo_matrix));
    if(coo == NULL){
        printf("Couldn't allocate memory for coo in mtx2comp.\n");
        exit(-1);
    }

    coo->col = coo_col;
    coo->row = coo_row;

    if (!strcmp(type, "csc")) {
        comp_matrix* csc = new_comp_matrix(2 * nz, N, "csc");
        coo_to_csc(csc, coo);
        free_coo(coo);
        return csc;
    }
    else if (!strcmp(type, "csr")) {
        comp_matrix* csr = new_comp_matrix(2 * nz, M, "csr");
        coo_to_csr(csr, coo);
        free_coo(coo);
        return csr;
    }
    else {
        printf("Unknown type in function mtx2comp.\n");
        free_coo(coo);
        exit(-1);
    }
}

#endif