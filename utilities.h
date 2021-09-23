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


//an olapane kala auto tha parei podi
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

typedef struct {
    comp_matrix** blocks;
    uint32_t* block_row;
    uint32_t* block_col;
    uint32_t nnz_blocks;
    uint32_t n_b;
    uint32_t real_dim;
} block_comp_matrix_2;



/* -----------------------Functions for the 2d matrices----------------------------*/


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


//Prints an array of type matrix_2d
void print_matrix_2d(matrix_2d* array) {
    for (int i = 0; i < array->rows; i++) {
        for (int j = 0; j < array->cols; j++) {
            printf("%d ", array->mat[i][j]);
        }
        printf("\n");
    }
}


//Free memory allocated for a 2d matrix
void free_matrix_2d(matrix_2d* array) {
    for (int i = 0; i < array->rows; i++) {
        free(array->mat[i]);
    }
    free(array->mat);
    free(array);
}


/*------------------------Functions for the CSR/CSC matrices-------------------------------*/


//Create a comp matrix struct
comp_matrix* new_comp_matrix(uint32_t nnz, uint32_t n, char* type) {

    comp_matrix* array = (comp_matrix*)malloc(sizeof(comp_matrix));
    if (array == NULL){
        printf("Couldn't allocate memory for array in new_comp_matrix.\n");
        exit(-1);
    } 

    //If we want to create a csc matrix
    if (!strcmp(type, "csc")) {
        array->col = (uint32_t*)malloc((n + 1) * sizeof(uint32_t));
        array->row = (uint32_t*)malloc(nnz * sizeof(uint32_t));
    }
    //If we want to create a csr matrix
    else if (!strcmp(type, "csr")) {
        array->col = (uint32_t*)malloc(nnz * sizeof(uint32_t));
        array->row = (uint32_t*)malloc((n + 1) * sizeof(uint32_t));
    }
    else {
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


//Function that turns a matrix_2d array to csr format
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


//Prints a csc matrix
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


//Prints a csr matrix
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


//Free memory allocated for a csr/csc matrix
void free_comp_matrix(comp_matrix* array) {
    free(array->col);
    free(array->row);
    free(array);
}

//an ola pane kala me tous neous blocked auta tha paroun podi

/*
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
*/

/*------------------- Functions for the blocked CSR/CSC matrices ------------------------*/

/**
 * Function that transforms a csr matrix into a blocked csr matrix.
 * The basic idea is that we use a csr format to store the blocks
 *  while each block is also in csr format.
**/
block_comp_matrix_2* csr_to_blocked(comp_matrix* array,uint32_t b){
    //Number of blocks in each dimension(both zero and non-zero blocks included)
    uint32_t n_b = (array->n)/b;

    uint32_t remaining_rows = (array->n)%b;

    // If n_b is not an integer one more row of blocks should be added
    if (remaining_rows != 0){
        n_b = n_b + 1;

        // Temporarily change the array->row variable
        // Basically we pad the array with zeros to fit the dimensions of the blocks
        array->row = (uint32_t*)realloc(array->row, (n_b*b + 1)*sizeof(uint32_t));
        if (array->row == NULL) {
            printf("Couldn't reallocate memory for array->row in csr_to_blocked.\n");
            free_comp_matrix(array);
            exit(-1);
        }
        
        for (int i = array->n + 1; i < n_b*b + 1; i++){
            array->row[i] = array->nnz;
        }
    }

    uint32_t* nnz = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (nnz == NULL){
        printf("Couldn't allocate memory for nnz in csr_to_blocked.\n");
        exit(-1);
    }

    uint32_t* col = (uint32_t*)calloc(array->nnz,sizeof(uint32_t));
    if (col == NULL){
        printf("Couldn't allocate memory for col in csr_to_blocked.\n");
        exit(-1);
    }

    uint32_t* row = (uint32_t*)calloc((b+1)*pow(n_b,2),sizeof(uint32_t));
    if (row == NULL){
        printf("Couldn't allocate memory for row in csr_to_blocked.\n");
        exit(-1);
    }

    uint32_t* block_row = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(block_row == NULL){
        printf("Couldn't allocate memory for block_row in csr_to_blocked.\n");
        exit(-1);
    }

    uint32_t* block_col = (uint32_t*)malloc(pow(n_b,2)*sizeof(uint32_t));
    if(block_row == NULL){
        printf("Couldn't allocate memory for block_col in csr_to_blocked.\n");
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
            
            //For each row included in the block
            for(uint32_t k=0;k<b;++k){
                //For each non-zero element of the row
                for(uint32_t l=array->row[i*b+k];l<array->row[i*b+k+1];++l){
                    //The non-zero element belongs to a following block in this row of blocks
                    if(array->col[l]>=(j+1)*b){
                        break;
                    }
                    //The non-zero element belongs to this block
                    else if(array->col[l]>=j*b){
                        nnz_in_block++;
                        col[nnz_found] = array->col[l];
                        nnz_found++;
                    }
                    //The non-zero element belongs to a previous block in this row of blocks
                    else{
                        continue;
                    }
                }

                row[nnz_blocks_count*(b+1)+k+1] = nnz_in_block;
            }
            
            //We have a non-zero block
            if(nnz_in_block>0){
                block_col[nnz_blocks_count] = j;
                nnz[nnz_blocks_count] = nnz_in_block;
                nnz_blocks_count++;
                nnz_in_row++;
            }
        }

        block_row[i+1] = block_row[i] + nnz_in_row;
    }

    //Readjust the length of the arrays that depend on the number of non-zero blocks
    nnz = (uint32_t*)realloc(nnz,nnz_blocks_count * sizeof(uint32_t));
    if(nnz == NULL){
        printf("Couldn't reallocate memory for nnz in csr_to_blocked.\n");
        exit(-1);
    }

    block_col = (uint32_t*)realloc(block_col, nnz_blocks_count * sizeof(uint32_t));
    if(block_col == NULL){
        printf("Couldn't reallocate memory for block_col in csr_to_blocked.\n");
        exit(-1);
    }

    row = (uint32_t*)realloc(row,nnz_blocks_count *(b + 1)*sizeof(uint32_t));
    if(row == NULL){
        printf("Couldn't reallocate memory for row in csr_to_blocked.\n");
        exit(-1);
    } 

    block_comp_matrix_2* blocked_matrix = (block_comp_matrix_2*)malloc(sizeof(block_comp_matrix_2));
    if(blocked_matrix == NULL){
        printf("Couldn't allocate memory for blocked_matrix in csr_to_blocked.\n");
        exit(-1);
    }   

    blocked_matrix->n_b = n_b;
    blocked_matrix->nnz_blocks = nnz_blocks_count;
    blocked_matrix->block_col = block_col;
    blocked_matrix->block_row = block_row;
    blocked_matrix->real_dim = array->n;

    nnz_found = 0;  //Non-zero elements found up to this point

    //Create the csr format of each non-zero block
    blocked_matrix->blocks = (comp_matrix**)malloc(nnz_blocks_count*sizeof(comp_matrix*));
    if(blocked_matrix->blocks == NULL){
        printf("Couldn't allocate memory for blocks in csr_to_blocked.\n");
        exit(-1);
    }

    for(int i=0;i<nnz_blocks_count;++i){
        blocked_matrix->blocks[i] = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(blocked_matrix->blocks[i] == NULL){
            printf("Couldn't allocate memory for blocks[%d] in csr_to_blocked.\n", i);
            exit(-1);
        }
        blocked_matrix->blocks[i]->n = b;
        blocked_matrix->blocks[i]->nnz = nnz[i];
        
        blocked_matrix->blocks[i]->col = (uint32_t*)malloc(nnz[i]*sizeof(uint32_t));
        if(blocked_matrix->blocks[i]->col == NULL){
            printf("Couldn't allocate memory for blocks[%d]->col in csr_to_blocked.\n", i);
            exit(-1);    
        }
        for(int j=0;j<nnz[i];++j){
            blocked_matrix->blocks[i]->col[j] = col[nnz_found+j];
        }
        nnz_found += nnz[i];

        blocked_matrix->blocks[i]->row = (uint32_t*)malloc((b+1)*sizeof(uint32_t));
        if(blocked_matrix->blocks[i]->row == NULL){
            printf("Couldn't allocate memory for blocks[%d]->row in csr_to_blocked.\n", i);
            exit(-1);    
        }
        for(int j=0;j<b+1;++j){
            blocked_matrix->blocks[i]->row[j] = row[i*(b+1)+j];
        }
    }

    // Return array->row array to its original size
    if (remaining_rows != 0){
        array->row = (uint32_t*)realloc(array->row, (array->n + 1)*sizeof(uint32_t));
        if (array->row == NULL) {
            printf("Couldn't reallocate memory for array->row in csr_to_blocked.\n");
            free_comp_matrix(array);
            exit(-1);
        }
    }

    free(row);
    free(col);
    free(nnz);

    return blocked_matrix;
}


/**
 * Similar to csr_to_blocked but for csc format.
 * Blocks are stored in csc format while each block is also in csc format.
**/
block_comp_matrix_2* csc_to_blocked(comp_matrix* array,uint32_t b){
    //Number of blocks in each dimension(both zero and non-zero blocks included)
    uint32_t n_b = (array->n)/b;
    uint32_t remaining_cols = (array->n) % b;

    // If n_b is not an integer one more row of blocks should be added
    if (remaining_cols != 0) {
        n_b = n_b + 1;
        // Temporarily change the array->col variable
        // Basically we pad the array with zeros to fit the dimensions of the blocks
        array->col = (uint32_t*)realloc(array->col, (n_b * b + 1) * sizeof(uint32_t));
        if (array->col == NULL) {
            printf("Couldn't reallocate memory for array->col in csc_to_blocked.\n");
            free_comp_matrix(array);
            exit(-1);
        }
        for (int i = array->n + 1; i < (n_b * b + 1); i++) {
            array->col[i] = array->nnz;
        }
    }

    uint32_t* nnz = (uint32_t*)calloc(pow(n_b,2),sizeof(uint32_t));
    if (nnz == NULL){
        printf("Couldn't allocate memory for nnz in csc_to_blocked.\n");
        exit(-1);
    }

    uint32_t* row = (uint32_t*)calloc(array->nnz,sizeof(uint32_t));
    if (row == NULL){
        printf("Couldn't allocate memory for row in csc_to_blocked.\n");
        exit(-1);
    }

    uint32_t* col = (uint32_t*)calloc((b+1)*pow(n_b,2),sizeof(uint32_t));
    if (col == NULL){
        printf("Couldn't allocate memory for col in csc_to_blocked.\n");
        exit(-1);
    }

    uint32_t* block_col = (uint32_t*)calloc(n_b+1,sizeof(uint32_t));
    if(block_col == NULL){
        printf("Couldn't allocate memory for block_col in csc_to_blocked.\n");
        exit(-1);
    }

    uint32_t* block_row = (uint32_t*)malloc(pow(n_b,2)*sizeof(uint32_t));
    if(block_row == NULL){
        printf("Couldn't allocate memory for block_row in csc_to_blocked.\n");
        exit(-1);
    }

    uint32_t nnz_in_block;          //Number of non-zero elements in each block
    uint32_t nnz_in_col;            //Number of non-zero blocks in each row of blocks  
    uint32_t nnz_found = 0;         //Number of non-zero elements found until this point
    uint32_t nnz_blocks_count = 0;  //Number of blocks with non-zero elements

    //First 2 for loops to check all blocks
    for(uint32_t i=0;i<n_b;++i){
        nnz_in_col = 0;
        for(uint32_t j=0;j<n_b;++j){
            nnz_in_block = 0;
            
            //For each row included in the block
            for(uint32_t k=0;k<b;++k){
                //For each non-zero element of the column
                for(uint32_t l=array->col[i*b+k];l<array->col[i*b+k+1];++l){
                    //The non-zero element belongs to a following block in this column of blocks
                    if(array->row[l]>=(j+1)*b){
                        break;
                    }
                    //The non-zero element belongs to this block
                    else if(array->row[l]>=j*b){
                        nnz_in_block++;
                        row[nnz_found] = array->row[l];
                        nnz_found++;
                    }
                    //The non-zero element belongs to a previous block in this column of blocks
                    else{
                        continue;
                    }
                }

                col[nnz_blocks_count*(b+1)+k+1] = nnz_in_block;
            }
            
            //We have a non-zero block
            if(nnz_in_block>0){
                block_row[nnz_blocks_count] = j;
                nnz[nnz_blocks_count] = nnz_in_block;
                nnz_blocks_count++;
                nnz_in_col++;
            }
        }

        block_col[i+1] = block_col[i] + nnz_in_col;
    }
    
    //readjust the length of the arrays that depend on the number of non-zero blocks
    nnz = (uint32_t*)realloc(nnz,nnz_blocks_count * sizeof(uint32_t));
    if(nnz == NULL){
        printf("Couldn't reallocate memory for nnz in csc_to_blocked.\n");
        exit(-1);
    }

    block_row = (uint32_t*)realloc(block_row, nnz_blocks_count * sizeof(uint32_t));
    if(block_row == NULL){
        printf("Couldn't reallocate memory for block_row in csc_to_blocked.\n");
        exit(-1);
    }

    col = (uint32_t*)realloc(col,nnz_blocks_count *(b + 1)*sizeof(uint32_t));
    if(col == NULL){
        printf("Couldn't reallocate memory for col in csc_to_blocked.\n");
        exit(-1);
    } 

    block_comp_matrix_2* blocked_matrix = (block_comp_matrix_2*)malloc(sizeof(block_comp_matrix_2));
    if(blocked_matrix == NULL){
        printf("Couldn't allocate memory for blocked_matrix in csc_to_blocked.\n");
        exit(-1);
    }   

    blocked_matrix->n_b = n_b;
    blocked_matrix->nnz_blocks = nnz_blocks_count;
    blocked_matrix->block_col = block_col;
    blocked_matrix->block_row = block_row;
    blocked_matrix->real_dim = array->n;

    nnz_found = 0;

    blocked_matrix->blocks = (comp_matrix**)malloc(nnz_blocks_count*sizeof(comp_matrix*));
    if(blocked_matrix->blocks == NULL){
        printf("Couldn't allocate memory for blocks in csc_to_blocked.\n");
        exit(-1);
    }

    for(int i=0;i<nnz_blocks_count;++i){
        blocked_matrix->blocks[i] = (comp_matrix*)malloc(sizeof(comp_matrix));
        if(blocked_matrix->blocks[i] == NULL){
            printf("Couldn't allocate memory for blocks[%d] in csc_to_blocked.\n", i);
            exit(-1);
        }
        blocked_matrix->blocks[i]->n = b;
        blocked_matrix->blocks[i]->nnz = nnz[i];
        
        blocked_matrix->blocks[i]->row = (uint32_t*)malloc(nnz[i]*sizeof(uint32_t));
        if(blocked_matrix->blocks[i]->row == NULL){
            printf("Couldn't allocate memory for blocks[%d]->row in csc_to_blocked.\n", i);
            exit(-1);    
        }
        for(int j=0;j<nnz[i];++j){
            blocked_matrix->blocks[i]->row[j] = row[nnz_found+j];
        }
        nnz_found += nnz[i];

        blocked_matrix->blocks[i]->col = (uint32_t*)malloc((b+1)*sizeof(uint32_t));
        if(blocked_matrix->blocks[i]->col == NULL){
            printf("Couldn't allocate memory for blocks[%d]->col in csc_to_blocked.\n", i);
            exit(-1);    
        }
        for(int j=0;j<b+1;++j){
            blocked_matrix->blocks[i]->col[j] = col[i*(b+1)+j];
        }
    }

    // Return array->col array to its original size
    if (remaining_cols != 0) {
        array->col = (uint32_t*)realloc(array->col, (array->n + 1) * sizeof(uint32_t));
        if (array->col == NULL) {
            printf("Couldn't reallocate memory for array->col in csc_to_blocked.\n");
            free_comp_matrix(array);
            exit(-1);
        }
    }

    free(col);
    free(row);
    free(nnz);

    return blocked_matrix;
}


//Function that transforms a matrix in blocked csr format to csr format
comp_matrix* blocked_to_csr(block_comp_matrix_2* blocked){

    uint32_t b = blocked->blocks[0]->n; //We assume that the matrix has at least one non-zero block
    uint32_t n_b = blocked->n_b;
    uint32_t n = b*n_b;             
    uint32_t blocks_in_row;         //Non-zero blocks in a row of blocks   
    uint32_t nnz_in_row_of_block;   //Non-zero elements in a row of elements in a block
    uint32_t nnz_count = 0;         //Non-zero elements found up to this point
    uint32_t block_row_start;       //Index of first non-zero block in a row of blocks
    uint32_t row_inside_block_start; //Index of first non-zero element in a row of a block

    uint32_t total_nnz = 0;         //Total non-zero elements
    for(uint32_t i=0;i<blocked->nnz_blocks;++i){
        total_nnz += blocked->blocks[i]->nnz;
    }

    comp_matrix* csr = (comp_matrix*)malloc(sizeof(comp_matrix));
    if(csr == NULL){
        printf("Couldn't allocate memory for csr in blocked_to_csr.\n");
        exit(-1);
    }

    csr->n = blocked->real_dim;
    csr->nnz = total_nnz;

    csr->row = (uint32_t*)calloc(n+1,sizeof(uint32_t));
    if(csr->row == NULL){
       printf("Couldn't allocate memory for row in blocked_to_csr.\n");
        exit(-1);
    } 

    csr->col = (uint32_t*)malloc(total_nnz*sizeof(uint32_t));
    if(csr->col == NULL){
       printf("Couldn't allocate memory for col in blocked_to_csr.\n");
        exit(-1);
    }

    //For every row of blocks
    for(uint32_t i=0;i<n_b;++i){

        //Index of first block of this row of blocks
        block_row_start = blocked->block_row[i];

        //Number of blocks in this row of blocks
        blocks_in_row = blocked->block_row[i+1]-blocked->block_row[i];

        //If there are no blocks in this row of blocks
        if(blocks_in_row == 0){
            for(uint32_t m=0;m<b;++m){
                //There are no non-zero elements in these rows of elements
                csr->row[i*b+m+1] = csr->row[i*b+m];
            }
            continue;
        }

        //For every row of elements in a row of blocks
        for(uint32_t j=0;j<b;++j){

            //For every block in this row of blocks
            for(uint32_t k=0;k<blocks_in_row;++k){

                //Index of first element in the j-th row of the k-th block in this row of blocks
                row_inside_block_start = blocked->blocks[block_row_start+k]->row[j];

                //Number of non-zero elements in the j-th row of the k-th block in this row of blocks
                nnz_in_row_of_block = blocked->blocks[block_row_start+k]->row[j+1]-blocked->blocks[block_row_start+k]->row[j];
                
                //If there are no elements in this row of elements of this block
                if(nnz_in_row_of_block == 0){
                    continue;
                }
                
                //For every element in this row of elements in this block
                for(uint32_t l=0;l<nnz_in_row_of_block;++l){
                    csr->col[nnz_count] = blocked->blocks[block_row_start+k]->col[row_inside_block_start+l];
                    nnz_count++;
                }

            }   
            csr->row[i*b+j+1] = nnz_count;   
        }
    }

    // If the original dimensions mismatch reallocate csr->row to its real size
    if(blocked->real_dim != n){
        csr->row = (uint32_t*)realloc(csr->row, (csr->n + 1)*sizeof(uint32_t));
        if (csr->row == NULL) {
            printf("Couldn't reallocate memory for csr->row in blocked_to_csr.\n");
            free_comp_matrix(csr);
            exit(-1);
        }
    }
    
    return csr;
}


/**
 * Function that transforms a matrix in blocked csc format to csc format
 * Works similarly as the blocked_to_csr function
**/
comp_matrix* blocked_to_csc(block_comp_matrix_2* blocked){

    uint32_t b = blocked->blocks[0]->n;
    uint32_t n_b = blocked->n_b;
    uint32_t n = b*n_b;
    uint32_t blocks_in_col;           
    uint32_t nnz_in_col_of_block;   
    uint32_t nnz_count = 0;         
    uint32_t block_col_start;
    uint32_t col_inside_block_start;

    uint32_t total_nnz = 0;
    for(uint32_t i=0;i<blocked->nnz_blocks;++i){
        total_nnz += blocked->blocks[i]->nnz;
    }

    comp_matrix* csc = (comp_matrix*)malloc(sizeof(comp_matrix));
    if(csc == NULL){
        printf("Couldn't allocate memory for csc in blocked_to_csc.\n");
        exit(-1);
    }

    csc->n = blocked->real_dim;
    csc->nnz = total_nnz;

    csc->col = (uint32_t*)calloc(n+1,sizeof(uint32_t));
    if(csc->col == NULL){
       printf("Couldn't allocate memory for col in blocked_to_csc.\n");
        exit(-1);
    } 

    csc->row = (uint32_t*)malloc(total_nnz*sizeof(uint32_t));
    if(csc->row == NULL){
       printf("Couldn't allocate memory for row in blocked_to_csc.\n");
        exit(-1);
    }

    //For every column of blocks
    for(uint32_t i=0;i<n_b;++i){

        //Index of first block of this column of blocks
        block_col_start = blocked->block_col[i];

        //Number of blocks in this column of blocks
        blocks_in_col = blocked->block_col[i+1]-blocked->block_col[i];

        //If there are no blocks in this column of blocks
        if(blocks_in_col == 0){
            //There are no non-zero elements in these columns of elements 
            for(uint32_t m=0;m<b;++m){
                csc->col[i*b+m+1] = csc->col[i*b+m];
            }
            continue;
        }

        //For every column of elements in a column of blocks
        for(uint32_t j=0;j<b;++j){

            //For every block in this column of blocks
            for(uint32_t k=0;k<blocks_in_col;++k){

                //Index of first element in the j-th column of the k-th block in this column of blocks
                col_inside_block_start = blocked->blocks[block_col_start+k]->col[j];

                //Number of non-zero elements in the j-th column of the k-th block in this column of blocks
                nnz_in_col_of_block = blocked->blocks[block_col_start+k]->col[j+1]-blocked->blocks[block_col_start+k]->col[j];
                
                //If there are no elements in this column of elements of this block
                if(nnz_in_col_of_block == 0){
                    continue;
                }
                
                //For every element in this column of elements in this block
                for(uint32_t l=0;l<nnz_in_col_of_block;++l){
                    csc->row[nnz_count] = blocked->blocks[block_col_start+k]->row[col_inside_block_start+l];
                    nnz_count++;
                }

            }   
            csc->col[i*b+j+1] = nnz_count;   
        }
    }

    // If the original dimensions mismatch reallocate csc->col to its real size
    if (blocked->real_dim != n) {
        csc->col = (uint32_t*)realloc(csc->col, (csc->n + 1) * sizeof(uint32_t));
        if (csc->col == NULL) {
            printf("Couldn't reallocate memory for csc->col in blocked_to_csc.\n");
            free_comp_matrix(csc);
            exit(-1);
        }
    }
    
    return csc;
}

    /**
     * Function that prints a blocked csr matrix.
     * First it prints the csr arrays we have in the block-level
     * and then prints the csr arrays of each non-zero block
     **/
    void print_blocked_csr(block_comp_matrix_2 * blocked_matrix) {
        printf("Printing a blocked csr matrix.\n");

        printf("Block col: ");
        for (int i = 0; i < blocked_matrix->nnz_blocks; ++i) {
            printf("%d ", blocked_matrix->block_col[i]);
        }
        printf("\n");

        printf("Block row: ");
        for (int i = 0; i < blocked_matrix->n_b + 1; ++i) {
            printf("%d ", blocked_matrix->block_row[i]);
        }
        printf("\n");

        for (int i = 0; i < blocked_matrix->nnz_blocks; ++i) {
            printf("Block %d:\n", i);
            print_csr(blocked_matrix->blocks[i]);
            printf("\n");
        }
        printf("\n");
}


/**
 * Function that prints a blocked csc matrix.
 * First it prints the csc arrays we have in the block-level
 * and then prints the csc arrays of each non-zero block
**/
void print_blocked_csc(block_comp_matrix_2* blocked_matrix){
    printf("Printing a blocked csc matrix.\n");

    printf("Block row: ");
    for(int i=0;i<blocked_matrix->nnz_blocks;++i){
        printf("%d ",blocked_matrix->block_row[i]);
    }
    printf("\n");

    printf("Block col: ");
    for(int i=0;i<blocked_matrix->n_b+1;++i){
        printf("%d ",blocked_matrix->block_col[i]);
    }
    printf("\n");

    for(int i=0;i<blocked_matrix->nnz_blocks;++i){
        printf("Block %d:\n",i);
        print_csc(blocked_matrix->blocks[i]);
        printf("\n");
    }
    printf("\n");
}


//Free allocated memory used for a matrix in blocked csr/csc format
void free_blocked_comp_matrix_2(block_comp_matrix_2* blocked_matrix){
    free(blocked_matrix->block_col);
    free(blocked_matrix->block_row);
    
    for(int i=0;i<blocked_matrix->nnz_blocks;++i){
        free_comp_matrix(blocked_matrix->blocks[i]);
    }
    free(blocked_matrix->blocks);

    free(blocked_matrix);
}

//TODO an ola leitourgoun kala, na dw an mporw kapws na sumpth3w tis sunarthseis
//wste na mhn einai 3exwrista gia tous csc kai tous csr

#endif