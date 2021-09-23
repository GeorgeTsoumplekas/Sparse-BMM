#include "utilities.h"
#include "bmm_seq.h"

void main(int argc, char* argv[]){
    if(argc<4){
        printf("Not enough arguments.\n");
        exit(-1);
    }

    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    char* filename_A = argv[1];
    char* filename_B = argv[2];
    char* filename_C = argv[3];

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    coo_matrix* A_coo = readCOO(filename_A);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A txt->COO: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    coo_matrix* B_coo = readCOO(filename_B);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B txt->COO: %.5f seconds.\n", elapsed);
    
/*
    print_coo(A_coo);
    printf("n=%u, nnz=%u\n",A_coo->n,A_coo->nnz);
    printf("\n");
*/

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* A_csr = coo2csr(A_coo);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A COO->CSR: %.5f seconds.\n", elapsed);
/*
    print_csr(A_csr);
    printf("\n");
*/
/*
    print_coo(B_coo);
    printf("n=%u, nnz=%u\n",B_coo->n,B_coo->nnz);
    printf("\n");
*/
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* B_csc = coo2csc(B_coo);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B COO->CSC: %.5f seconds.\n", elapsed);

/*
    print_csc(B_csc);
    printf("\n");
*/

    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* C_csr = bmm_seq(A_csr,B_csc,0);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for bmm: %.5f seconds.\n", elapsed);

/*
    print_csr(C_csr);
    printf("\n");
*/

    uint32_t check = check_result(filename_C,C_csr);
    if(check==0){
        printf("Wrong result.\n");
    }
    else{
        printf("Correct result.\n");
    }

    free_coo(A_coo);
    free_coo(B_coo);

    free_comp_matrix(A_csr);
    free_comp_matrix(B_csc);
    free_comp_matrix(C_csr);

}