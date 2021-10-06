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


//Struct that holds a dense 2-dimensional matrix
typedef struct {
    uint32_t** mat;
    uint32_t rows;
    uint32_t cols;
} matrix_2d;


//Struct that holds a matrix in coo format
typedef struct {
    uint32_t* row;
    uint32_t* col;
    uint32_t nnz;
    uint32_t n;
} coo_matrix;


/**
 * Struct that holds a blocked csc/csr matrix
 * We hold only the non-zero blocks. Their position is stored is csc(csr} format but
 * on a block level this time. Moreover, each block is also stored in csc(csr) format.
**/
typedef struct {
    comp_matrix** blocks;   //array containing the non-zero blocks in csc(csr) format
    uint32_t* block_row;
    uint32_t* block_col;
    uint32_t nnz_blocks;
    uint32_t n_b;
    uint32_t real_dim;
} block_comp_matrix;



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


/*------------------------Functions for COO matrices--------------------------------*/


/**
 * Reads a .txt file which contains a sparse square matrix in COO format
 * and creates a coo matrix struct
**/
coo_matrix* readCOO(char* filename){
    uint32_t n;
    uint32_t nnz;

    coo_matrix* coo = (coo_matrix*)malloc(sizeof(coo_matrix));
    if(coo == NULL){
        printf("Couldn't allocate memory for coo in readCOO.\n");
        exit(-1);
    }

    FILE *stream = fopen(filename,"r");
    if(stream == NULL){
        printf("Error opening %s\n",filename);
        exit(-1);
    }

    //First line of the file contains the size of the square matrix and its number of non-zero elements
    int ign = fscanf(stream, "%u,%u\n", &n, &nnz);

    coo->row = (uint32_t*)malloc(nnz*sizeof(uint32_t));
    if(coo->row == NULL){
        printf("Couldn't allocate memory for coo->row in readCOO.\n");
        exit(-1);
    }

    coo->col = (uint32_t*)malloc(nnz*sizeof(uint32_t));
    if(coo->col == NULL){
        printf("Couldn't allocate memory for coo->col in readCOO.\n");
        exit(-1);
    }

    for(uint32_t i=0;i<nnz;++i){
        int ign = fscanf(stream, "%u,%u\n", &coo->row[i], &coo->col[i]);
    }

    coo->n = n;
    coo->nnz = nnz;

    fclose(stream);

    return coo;
}


//Prints a COO matrix
void print_coo(coo_matrix* coo){
    printf("Printing a COO matrix.\n");
    for(uint32_t i=0;i<coo->nnz;++i){
        printf("(%u, %u)\n",coo->row[i],coo->col[i]);
    }
}


//Frees memory allocated for a COO matrix
void free_coo(coo_matrix* coo){
    free(coo->col);
    free(coo->row);
    free(coo);
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
        array->col = (uint32_t*)calloc((n + 1), sizeof(uint32_t));
        array->row = (uint32_t*)calloc(nnz, sizeof(uint32_t));
    }
    //If we want to create a csr matrix
    else if (!strcmp(type, "csr")) {
        array->col = (uint32_t*)calloc(nnz, sizeof(uint32_t));
        array->row = (uint32_t*)calloc((n + 1), sizeof(uint32_t));
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


/**
 * Function that turns a matrix from coo to csr format
 * Note for this function to be correct, the non-zero elements should be sorted first by row
 * and then by column
**/
comp_matrix* coo2csr(coo_matrix* coo){
    comp_matrix* csr = (comp_matrix*)malloc(sizeof(comp_matrix));
    if(csr == NULL){
        printf("Couldn't allocate memory for csr in coo2csr./n");
        exit(-1);
    }

    csr->n = coo->n;
    csr->nnz = coo->nnz;

    csr->col = (uint32_t*)malloc(csr->nnz*sizeof(uint32_t));
    if(csr->col == NULL){
        printf("Couldn't allocate memory for csr->col in coo2csr./n");
        exit(-1);
    }

    csr->row = (uint32_t*)calloc(csr->n+1,sizeof(uint32_t));
    if(csr->row == NULL){
        printf("Couldn't allocate memory for csr->row in coo2csr./n");
        exit(-1);
    }

    uint32_t prev_elem_row = 0;     //Row in which the previous non-zero element belongs to
    uint32_t row_count = 1;         //Index of current row

    //For each non-zero element
    for(uint32_t i=0;i<csr->nnz;++i){
        csr->col[i] = coo->col[i]-1; //Because csr matrices are 0-based

        //If current element is the first element of a new row
        if(coo->row[i]-1 > prev_elem_row){
            //Checking if one or more consecutive rows have zero elements only
            if(coo->row[i]-1-prev_elem_row > 1){
                //Then the row index of these rows is the same as the last one's with a non-zero element
                for(uint32_t j=0;j<(coo->row[i]-1-prev_elem_row);++j){
                    csr->row[row_count] = i;
                    row_count++;
                }
            }
            //If there are no zero-elements-only rows in between
            else{
                csr->row[row_count] = i;
                row_count++;
            }

            csr->row[row_count] = i+1;
        }
        //If current non-zero element belongs to the same row as its previous one
        else{
            csr->row[row_count] = i+1;
        }


        prev_elem_row = coo->row[i]-1;
    }

    //If there are rows with zero elements only after the last non-zero element
    while(row_count<coo->n+1){
        csr->row[row_count] = coo->nnz;
        row_count++;
    }

    return csr;
}


/**
 * Function that turns a matrix from coo to csc format
 * Note for this function to be correct, the non-zero elements should be sorted first by column
 * and then by row
**/
comp_matrix* coo2csc(coo_matrix* coo){
    comp_matrix* csc = (comp_matrix*)malloc(sizeof(comp_matrix));
    if(csc == NULL){
        printf("Couldn't allocate memory for csc in coo2csc./n");
        exit(-1);
    }

    csc->n = coo->n;
    csc->nnz = coo->nnz;

    csc->row = (uint32_t*)malloc(csc->nnz*sizeof(uint32_t));
    if(csc->row == NULL){
        printf("Couldn't allocate memory for csc->row in coo2csc./n");
        exit(-1);
    }

    csc->col = (uint32_t*)calloc(csc->n+1,sizeof(uint32_t));
    if(csc->col == NULL){
        printf("Couldn't allocate memory for csc->col in coo2csc./n");
        exit(-1);
    }

    uint32_t prev_elem_col = 0;     //Column in which the previous non-zero element belongs
    uint32_t col_count = 1;         //Index of current column

    //For each non-zero element
    for(uint32_t i=0;i<csc->nnz;++i){
        csc->row[i] = coo->row[i]-1; //Because csc matrices are 0-based

        //If current element is the first element of a new column
        if(coo->col[i]-1 > prev_elem_col){
            //Checking if one or more consecutive columns have zero elements only
            if(coo->col[i]-1-prev_elem_col > 1){
                //Then the col index of these columns is the same as the last one's with a non-zero element
                for(uint32_t j=0;j<(coo->col[i]-1-prev_elem_col);++j){
                    csc->col[col_count] = i;
                    col_count++;
                }
            }
            //If there are no zero elements only columns in between
            else{
                csc->col[col_count] = i;
                col_count++;
            }

            csc->col[col_count] = i+1;
        }
        //If current non-zero element belongs to the same column as its previous one
        else{
            csc->col[col_count] = i+1;
        }


        prev_elem_col = coo->col[i]-1;
    }

    //If there are rows with zero elements only after the last non-zero element
    while(col_count<coo->n+1){
        csc->col[col_count] = coo->nnz;
        col_count++;
    }

    return csc;
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


/*------------------- Functions for the blocked CSR/CSC matrices ------------------------*/


/**
 * Function that transforms a csr matrix into a blocked csr matrix.
 * The basic idea is that we use a csr format to store the blocks
 *  while each block is also in csr format.
**/
block_comp_matrix* csr2blocked(comp_matrix* array,uint32_t b){
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

    block_comp_matrix* blocked_matrix = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
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
block_comp_matrix* csc2blocked(comp_matrix* array,uint32_t b){
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

    block_comp_matrix* blocked_matrix = (block_comp_matrix*)malloc(sizeof(block_comp_matrix));
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
comp_matrix* blocked2csr(block_comp_matrix* blocked){

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
comp_matrix* blocked2csc(block_comp_matrix* blocked){

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
void print_blocked_csr(block_comp_matrix * blocked_matrix) {
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
void print_blocked_csc(block_comp_matrix* blocked_matrix){
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
void free_block_comp_matrix(block_comp_matrix* blocked_matrix){
    free(blocked_matrix->block_col);
    free(blocked_matrix->block_row);
    
    for(int i=0;i<blocked_matrix->nnz_blocks;++i){
        free_comp_matrix(blocked_matrix->blocks[i]);
    }
    free(blocked_matrix->blocks);

    free(blocked_matrix);
}


/*------------------------ Other functions ------------------------------*/


/**
 * Function that checks if the result created by the multiplication of the two matrices is correct.
 * This is done by transforming the true result in CSR format and comparing it with the one that
 * we have calculated.
 * Returns 1 if result is correct, else 0.
**/
uint32_t check_result(char* filename_C,comp_matrix* C){
    
    //Transform real matrix C in COO format
    coo_matrix* C_coo_real = readCOO(filename_C);

    //Transform real matrix C from COO to CSR format
    comp_matrix* C_csr_real = coo2csr(C_coo_real);

    //Compare the two structs element by element
    if(C->nnz != C_csr_real->nnz){
        free_coo(C_coo_real);
        free_comp_matrix(C_csr_real);
        return 0;
    }

    if(C->n != C_csr_real->n){
        free_coo(C_coo_real);
        free_comp_matrix(C_csr_real);
        return 0;
    }

    for(uint32_t i=0;i<C->nnz;++i){
        if(C->col[i] != C_csr_real->col[i]){
            free_coo(C_coo_real);
            free_comp_matrix(C_csr_real);
            return 0;
        }
    }

    for(uint32_t i=0;i<C->n+1;++i){
        if(C->row[i] != C_csr_real->row[i]){
            free_coo(C_coo_real);
            free_comp_matrix(C_csr_real);
            return 0;
        }
    }

    free_coo(C_coo_real);
    free_comp_matrix(C_csr_real);

    return 1;
}

#endif