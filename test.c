#include "bmm_seq.h"

int main(int argc, char* argv[]){

    uint32_t rows = 4;
    uint32_t cols = 4;

    matrix_2d* array1 = (matrix_2d*)malloc(sizeof(matrix_2d));
    if (array1 == NULL){
        printf("Couldn't allocate memory for array in rand_matrix_2d\n");
        exit(-1);
    } 

    array1->mat = (uint32_t**)malloc(rows * sizeof(uint32_t*));
    if (array1->mat == NULL){
        printf("Couldn't allocate memory for array->mat in rand_matrix_2d\n");
        exit(-1);
    }

    for (uint32_t i = 0; i < rows; i++) {
        array1->mat[i] = (uint32_t*)malloc(cols * sizeof(uint32_t));
        if (array1->mat[i] == NULL){
            printf("Couldn't allocate memory for array->mat[%d] in rand_matrix_2d\n",i);
            exit(-1);
        }   
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array1->mat[i][j] = 0;
        }
    }

    array1->mat[0][2] = 1;
    array1->mat[1][3] = 1;
    array1->mat[2][1] = 1;
    array1->mat[2][2] = 1;
    array1->mat[3][3] = 1;

    array1->cols = cols;
    array1->rows = rows;

    print_matrix_2d(array1);
    printf("\n");

    comp_matrix* test_matrix_1 = matrix2csr(array1);
    print_csr(test_matrix_1);
    printf("\n");

    uint32_t b = 2;

    block_comp_matrix* test_matrix_2 = extract_blocked_csr(test_matrix_1,b);
    print_blocked_csr(test_matrix_2);
    printf("\n");

    matrix_2d* array2 = (matrix_2d*)malloc(sizeof(matrix_2d));
    if (array2 == NULL){
        printf("Couldn't allocate memory for array in rand_matrix_2d\n");
        exit(-1);
    } 

    array2->mat = (uint32_t**)malloc(rows * sizeof(uint32_t*));
    if (array2->mat == NULL){
        printf("Couldn't allocate memory for array->mat in rand_matrix_2d\n");
        exit(-1);
    }

    for (uint32_t i = 0; i < rows; i++) {
        array2->mat[i] = (uint32_t*)malloc(cols * sizeof(uint32_t));
        if (array2->mat[i] == NULL){
            printf("Couldn't allocate memory for array->mat[%d] in rand_matrix_2d\n",i);
            exit(-1);
        }   
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array2->mat[i][j] = 0;
        }
    }

    array2->mat[0][2] = 1;
    array2->mat[3][1] = 1;

    array2->cols = cols;
    array2->rows = rows;

    print_matrix_2d(array2);
    printf("\n");

    comp_matrix* test_matrix_3 = matrix2csc(array2);
    print_csc(test_matrix_3);
    printf("\n");

    block_comp_matrix* test_matrix_4 = extract_blocked_csc(test_matrix_3,b);
    print_blocked_csc(test_matrix_4);
    printf("\n");

/*
    matrix_2d* final_2d = bmm_seq_2d(array1,array2);
    print_matrix_2d(final_2d);
    printf("\n");

    comp_matrix* transf = matrix2csc(final_2d);
    print_csc(transf);
    printf("\n");

    comp_matrix* final_csc = bmm_seq(test_matrix_1,test_matrix_2);
    print_csc(final_csc);
    printf("\n");
*/    
    free_comp_matrix(test_matrix_1);
    free_comp_matrix(test_matrix_3);
//    free_comp_matrix(final_csc);
//    free_comp_matrix(transf);

//    free_matrix_2d(final_2d);
    free_matrix_2d(array1);
    free_matrix_2d(array2);

    free_blocked_comp_matrix(test_matrix_2);
    free_blocked_comp_matrix(test_matrix_4);
}