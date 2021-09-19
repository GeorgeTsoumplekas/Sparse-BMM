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

    printf("A:\n");
    print_matrix_2d(array1);
    printf("\n");

    //comp_matrix* test_matrix_1 = NULL;

    comp_matrix* test_matrix_1 = matrix2csr(array1);
    printf("A in csr:\n");
    print_csr(test_matrix_1);
    printf("\n");

    uint32_t b = 2;
/*
    block_comp_matrix* test_matrix_2 = extract_blocked_csr(test_matrix_1,b);
    print_blocked_csr(test_matrix_2);
    printf("\n");
*/
    block_comp_matrix_2* big_test = csr_to_blocked(test_matrix_1,b);
    printf("A in blocked csr:\n");
    print_blocked_csr(big_test);
    printf("\n");

    comp_matrix* return_1 = blocked_to_csr(big_test);
    printf("A transformed back to csr:\n");
    print_csr(return_1);
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

    printf("B:\n");
    print_matrix_2d(array2);
    printf("\n");

    comp_matrix* test_matrix_3 = matrix2csc(array2);
    printf("B in csc:\n");
    print_csc(test_matrix_3);
    printf("\n");
/*
    comp_matrix* test_matrix_3 = matrix2csr(array2);
    print_csr(test_matrix_3);
    printf("\n");
*/
/*
    block_comp_matrix* test_matrix_4 = extract_blocked_csc(test_matrix_3,b);
    print_blocked_csc(test_matrix_4);
    printf("\n");
*/

    block_comp_matrix_2* big_test_2 = csc_to_blocked(test_matrix_3,b);
    printf("B in blocked csc:\n");
    print_blocked_csc(big_test_2);
    printf("\n");

    comp_matrix* return_2 = blocked_to_csc(big_test_2);
    printf("B transormed back to csc:\n");
    print_csc(return_2);
    printf("\n");    


    matrix_2d* final_2d = bmm_seq_2d(array1,array2);
    printf("C:\n");
    print_matrix_2d(final_2d);
    printf("\n");

    comp_matrix* transf = matrix2csr(final_2d);
    printf("C (2d -> csr):\n");
    print_csr(transf);
    printf("\n");

    comp_matrix* final_csr = bmm_seq(test_matrix_1,test_matrix_3,0);
    if(final_csr == NULL){
        printf("Found NULL correctly.\n");
        return 0;
    }
    printf("C (csr multiplication):\n");
    print_csr(final_csr);
    printf("\n");

    block_comp_matrix_2* final_blocked = blocked_bmm_seq(big_test,big_test_2);
    printf("C blocked (blocked multiplication):\n");
    print_blocked_csr(final_blocked);
    printf("\n");

    comp_matrix* final_back_csr = blocked_to_csr(final_blocked);
    printf("C (blocked multiplication) transformed from blocked to csr:\n");
    print_csr(final_back_csr);
    printf("\n");
  
/*
    test_matrix_1 = block_union(test_matrix_1,test_matrix_3);
    printf("After union: \n");
    print_csr(test_matrix_1);
*/
    free_comp_matrix(test_matrix_1);
    free_comp_matrix(test_matrix_3);
    free_comp_matrix(return_1);
    free_comp_matrix(return_2);
    free_comp_matrix(transf);
    free_comp_matrix(final_csr);
    free_comp_matrix(final_back_csr);

    free_matrix_2d(final_2d);
    free_matrix_2d(array1);
    free_matrix_2d(array2);

    free_blocked_comp_matrix_2(big_test);
    free_blocked_comp_matrix_2(big_test_2);
    free_blocked_comp_matrix_2(final_blocked);

}