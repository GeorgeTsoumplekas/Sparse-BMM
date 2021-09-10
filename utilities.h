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
    uint32_t** mat;
    uint32_t rows;
    uint32_t cols;
} matrix_2d;

typedef struct {
    uint32_t* col;
    uint32_t* row;
    uint32_t* nnz;
    uint32_t* offsets;
    uint32_t* line_blocks;
    uint32_t b;
    uint32_t nb;
    uint32_t total_nnz;
    uint32_t nnz_blocks;
} block_comp_matrix;

/* ......   Functions   ......*/
// Functions for the comp_matrix struct

//Create a comp matrix struct
comp_matrix* new_comp_matrix(uint32_t nnz, uint32_t n, char* type) {

    comp_matrix* array = (comp_matrix*)malloc(sizeof(comp_matrix));
    if (array == NULL){
        printf("Couldn't allocate memory for array in new_comp_matrix.\n");
        exit(-1);
    } 

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
        printf("Couldn't allocate memory for array->col or array->row in new_comp_matrix.\n");
        free(array);
        exit(-1);
    }

    array->nnz = nnz;
    array->n = n;

    return array;
}

void free_comp_matrix(comp_matrix* array) {
    free(array->col);
    free(array->row);
    free(array);
}

/* Functions for the matrix_2d struct */
// Creates a random 2d matrix
matrix_2d* rand_matrix_2d(int rows, int cols) {
    matrix_2d* array = (matrix_2d*)malloc(sizeof(matrix_2d));
    if (array == NULL){
        printf("Couldn't allocate memory for array in rand_matrix_2d\n");
        exit(-1);
    } 

    array->mat = (uint32_t**)malloc(rows * sizeof(uint32_t*));
    if (array->mat == NULL){
        printf("Couldn't allocate memory for array->mat in rand_matrix_2d\n");
        exit(-1);
    }

    for (uint32_t i = 0; i < rows; i++) {
        array->mat[i] = (uint32_t*)malloc(cols * sizeof(uint32_t));
        if (array->mat[i] == NULL){
            printf("Couldn't allocate memory for array->mat[%d] in rand_matrix_2d\n",i);
            exit(-1);
        }   
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array->mat[i][j] = (uint32_t)(rand() % 2);
        }
    }

    array->cols = cols;
    array->rows = rows;

    return array;
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
            if (array->mat[i][j] != 0){
                nnz++;
            } 
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
            if (array->mat[i][j] != 0){
                nnz++;
            }
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

//Function that takes a matrix in csr format and creates its blocked csr format.
//b is the size of each block
block_comp_matrix* extract_blocked_csr(comp_matrix* array,uint32_t b){
    //number of blocks in each dimension(both zero and non-zero blocks included)
    uint32_t n_b = (array->n)/b;
    
    block_comp_matrix* blocked_csr = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if (blocked_csr == NULL){
        printf("Couldn't allocate memory for blocked_csr in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->nnz = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (blocked_csr->nnz == NULL){
        printf("Couldn't allocate memory for blocked_csr->nnz in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->offsets = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (blocked_csr->offsets == NULL){
        printf("Couldn't allocate memory for blocked_csr->offsets in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->col = (uint32_t*)calloc(array->nnz,sizeof(uint32_t));
    if (blocked_csr->col == NULL){
        printf("Couldn't allocate memory for blocked_csr->col in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->row = (uint32_t*)calloc((b+1)*pow(n_b,2),sizeof(uint32_t));
    if (blocked_csr->row == NULL){
        printf("Couldn't allocate memory for blocked_csr->row in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->line_blocks = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(blocked_csr->line_blocks==NULL){
        printf("Couldn't allocate memory for line_blocks in extract_blocked_csr.\n");
        exit(-1);
    }

    uint32_t nnz_in_block;      //number of non-zero elements in each block
    uint32_t nnz_in_row;        //number of non-zero blocks in each row of blocks  
    uint32_t nnz_found = 0;     //number of non-zero elements we have found until this point
    uint32_t nnz_blocks_count = 0; //number of blocks with non-zero elements

    //First 2 for loops to check all blocks
    for(uint32_t i=0;i<n_b;++i){
        nnz_in_row = 0;
        for(uint32_t j=0;j<n_b;++j){
            nnz_in_block = 0;
            
            //for each row included in the block
            for(uint32_t k=0;k<b;++k){
                //for each non-zero element of the row
                for(uint32_t l=array->row[i*b+k];l<array->row[i*b+k+1];++l){
                    //the non-zero element belongs to a following block in this row of blocks
                    if(array->col[l]>=(j+1)*b){
                        break;
                    }
                    //the non-zero element belongs to this block
                    else if(array->col[l]>=j*b){
                        nnz_in_block++;
                        blocked_csr->col[nnz_found] = array->col[l];
                        nnz_found++;
                    }
                    //the non-zero element belongs to a previous block in this row of blocks
                    else{
                        continue;
                    }
                }

                blocked_csr->row[nnz_blocks_count*(b+1)+k+1] = nnz_in_block;
            }
            
            //we have a non-zero block
            if(nnz_in_block>0){
                blocked_csr->offsets[nnz_blocks_count] =(uint32_t)(i*n_b+j);
                blocked_csr->nnz[nnz_blocks_count] = nnz_in_block;
                nnz_blocks_count++;
                nnz_in_row++;
            }
        }

        blocked_csr->line_blocks[i+1] = blocked_csr->line_blocks[i] + nnz_in_row;
    }

    blocked_csr->nnz_blocks = nnz_blocks_count;
    blocked_csr->total_nnz = array->nnz;
    blocked_csr->nb = n_b;
    blocked_csr->b = b;

    //readjust the length of the arrays that depend on the number of non-zero blocks
    blocked_csr->nnz = realloc(blocked_csr->nnz,blocked_csr->nnz_blocks * sizeof(uint32_t));
    if(blocked_csr->nnz == NULL){
        printf("Couldn't reallocate memory for nnz in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->offsets = realloc(blocked_csr->offsets, blocked_csr->nnz_blocks * sizeof(uint32_t));
    if(blocked_csr->offsets == NULL){
        printf("Couldn't reallocate memory for offsets in extract_blocked_csr.\n");
        exit(-1);
    }

    blocked_csr->row = realloc(blocked_csr->row,blocked_csr->nnz_blocks *(blocked_csr->b + 1)*sizeof(uint32_t));
    if(blocked_csr->row == NULL){
        printf("Couldn't reallocate memory for row in extract_blocked_csr.\n");
        exit(-1);
    }

    return blocked_csr;
}

//Function that takes a matrix in csc format and creates its blocked csr format.
//b is the size of each block
//Works in the same way as extract_blocked_csr with the role of col and row arrays inverted
block_comp_matrix* extract_blocked_csc(comp_matrix* array,uint32_t b){
    uint32_t n_b = (array->n)/b;
    
    block_comp_matrix* blocked_csc = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
    if (blocked_csc == NULL){
        printf("Couldn't allocate memory for blocked_csc in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->nnz = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (blocked_csc->nnz == NULL){
        printf("Couldn't allocate memory for blocked_csc->nnz in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->offsets = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (blocked_csc->offsets == NULL){
        printf("Couldn't allocate memory for blocked_csc->offsets in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->row = (uint32_t*)calloc(array->nnz,sizeof(uint32_t));
    if (blocked_csc->row == NULL){
        printf("Couldn't allocate memory for blocked_csc->row in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->col = (uint32_t*)calloc((b+1)*pow(n_b,2),sizeof(uint32_t));
    if (blocked_csc->col == NULL){
        printf("Couldn't allocate memory for blocked_csc->col in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->line_blocks = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(blocked_csc->line_blocks==NULL){
        printf("Couldn't allocate memory for line_blocks in extract_blocked_csc.\n");
        exit(-1);
    }

    uint32_t nnz_in_block;
    uint32_t nnz_in_col;    //number of non-zero blocks in each column of blocks
    uint32_t nnz_found = 0;
    uint32_t nnz_blocks_count = 0;

    for(uint32_t i=0;i<n_b;++i){
        nnz_in_col = 0;
        for(uint32_t j=0;j<n_b;++j){
            nnz_in_block = 0;

            for(uint32_t k=0;k<b;++k){
                for(uint32_t l=array->col[i*b+k];l<array->col[i*b+k+1];++l){

                    if(array->row[l]>=(j+1)*b){
                        break;
                    }
                    else if(array->row[l]>=j*b){
                        nnz_in_block++;
                        blocked_csc->row[nnz_found] = array->row[l];
                        nnz_found++;
                    }
                    else{
                        continue;
                    }
                }

                blocked_csc->col[nnz_blocks_count*(b+1)+k+1] = nnz_in_block;
            }
            
            if(nnz_in_block>0){
                blocked_csc->offsets[nnz_blocks_count] =(uint32_t)(i*n_b+j);
                blocked_csc->nnz[nnz_blocks_count] = nnz_in_block;
                nnz_blocks_count++;
                nnz_in_col++;
            }
        }
        blocked_csc->line_blocks[i+1] = blocked_csc->line_blocks[i] + nnz_in_col;
    }

    blocked_csc->nnz_blocks = nnz_blocks_count;
    blocked_csc->total_nnz = array->nnz;
    blocked_csc->nb = n_b;
    blocked_csc->b = b;

    blocked_csc->nnz = realloc(blocked_csc->nnz,blocked_csc->nnz_blocks * sizeof(uint32_t));
    if(blocked_csc->nnz == NULL){
        printf("Couldn't reallocate memory for nnz in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->offsets = realloc(blocked_csc->offsets, blocked_csc->nnz_blocks * sizeof(uint32_t));
    if(blocked_csc->offsets == NULL){
        printf("Couldn't reallocate memory for offsets in extract_blocked_csc.\n");
        exit(-1);
    }

    blocked_csc->col = realloc(blocked_csc->col,blocked_csc->nnz_blocks *(blocked_csc->b + 1)*sizeof(uint32_t));
    if(blocked_csc->col == NULL){
        printf("Couldn't reallocate memory for col in extract_blocked_csc.\n");
        exit(-1);
    }

    return blocked_csc;
}

//Prints a blocked csr matrix
void print_blocked_csr(block_comp_matrix* blocked_csr){
    uint32_t b = blocked_csr->b;
    uint32_t n_b = blocked_csr->nb;
    uint32_t nnz_blocks = blocked_csr->nnz_blocks;

    printf("Printing a blocked csr matrix.\n");
    printf("\nCol index: ");
    for(uint32_t i=0;i<blocked_csr->total_nnz;++i){
        printf("%d ", blocked_csr->col[i]);
    }
    printf("\n");

    printf("\nRow index: ");
    for(uint32_t i=0;i<nnz_blocks*(b+1);++i){
        printf("%d ", blocked_csr->row[i]);
    }
    printf("\n");

    printf("\nOffsets: ");
    for(uint32_t i=0;i<nnz_blocks;++i){
        printf("%d ", blocked_csr->offsets[i]);
    }
    printf("\n");

    printf("\nNNZ: ");
    for(uint32_t i=0;i<nnz_blocks;++i){
        printf("%d ", blocked_csr->nnz[i]);
    }
    printf("\n");

    printf("\nline_blocks: ");
    for(uint32_t i=0;i<n_b+1;++i){
        printf("%d ",blocked_csr->line_blocks[i]);
    }
    printf("\n");

    printf("\nb=%d, nb=%d, total_nnz=%d, nnz_blocks=%d\n",blocked_csr->b,blocked_csr->nb,blocked_csr->total_nnz,blocked_csr->nnz_blocks);
}

//Prints a blocked csc matrix
void print_blocked_csc(block_comp_matrix* blocked_csc){
    uint32_t b = blocked_csc->b;
    uint32_t n_b = blocked_csc->nb;
    uint32_t nnz_blocks = blocked_csc->nnz_blocks;

    printf("Printing a blocked csc matrix.\n");
    printf("\nRow index: ");
    for(uint32_t i=0;i<blocked_csc->total_nnz;++i){
        printf("%d ", blocked_csc->row[i]);
    }
    printf("\n");

    printf("\nCol index: ");
    for(uint32_t i=0;i<nnz_blocks*(b+1);++i){
        printf("%d ", blocked_csc->col[i]);
    }
    printf("\n");

    printf("\nOffsets: ");
    for(uint32_t i=0;i<nnz_blocks;++i){
        printf("%d ", blocked_csc->offsets[i]);
    }
    printf("\n");

    printf("\nNNZ: ");
    for(uint32_t i=0;i<nnz_blocks;++i){
        printf("%d ", blocked_csc->nnz[i]);
    }
    printf("\n");

    printf("\nline_blocks: ");
    for(uint32_t i=0;i<n_b+1;++i){
        printf("%d ",blocked_csc->line_blocks[i]);
    }
    printf("\n");

    printf("\nb=%d, nb=%d, total_nnz=%d, nnz_blocks=%d\n",blocked_csc->b,blocked_csc->nb,blocked_csc->total_nnz,blocked_csc->nnz_blocks);
}

//Frees memory used by a blocked matrix
void free_blocked_comp_matrix(block_comp_matrix* blocked_csr){
    free(blocked_csr->col);
    free(blocked_csr->row);
    free(blocked_csr->nnz);
    free(blocked_csr->offsets);
    free(blocked_csr->line_blocks);
    free(blocked_csr);
}

#endif