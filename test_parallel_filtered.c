#include "blocked_bmm_parallel.h"
#include "utilities.h"

int main(int argc, char* argv[]) {
    printf("\nThis file completes the filtered bmm with OpenMP parallelization with simple and blocked matrices.\n");

    if (argc < 7) {
        printf("Not enough arguments.\n");
        printf("Usage: %s [txt file of matrix A] [txt file of matrix B] ", argv[0]);
        printf("[txt file of matrix F] [txt file of matrix C] [no. threads] [b]\n");
        exit(-1);
    }

    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    char* filename_A = argv[1];
    char* filename_B = argv[2];
    char* filename_F = argv[3];
    char* filename_C = argv[4];

    int thread_num = atoi(argv[5]);
    printf("You have chosen %d threads.\n", thread_num);

    uint32_t b = atoi(argv[6]);
    printf("b=%d\n", b);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    coo_matrix* A_coo = readCOO(filename_A);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A txt->COO: %.5f seconds.\n", elapsed);

    //It is necessary that b * # of threads less or equal than n
    if(b*thread_num>A_coo->n){
        printf("ERROR: Condition that b * # of threads <= n is not held. Please try again with differend parameters.\n");
        exit(-1);
    }

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    coo_matrix* B_coo = readCOO(filename_B);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B txt->COO: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    coo_matrix* F_coo = readCOO(filename_F);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for F txt->COO: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* A_csr = coo2csr(A_coo);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A COO->CSR: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* B_csc = coo2csc(B_coo);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B COO->CSC: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* F_csr = coo2csr(F_coo);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for F COO->CSR: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* A_blocked = csr2blocked(A_csr, b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A CSR->blocked: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* B_blocked = csc2blocked(B_csc, b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B CSC->blocked: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* F_blocked = csr2blocked(F_csr, b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for F CSR->blocked: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* C_csr = nonblocked_bmm_parallel_filtered(A_csr, B_csc, F_csr, thread_num);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for parallel filtered bmm: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* C_blocked = blocked_bmm_parallel_filtered(A_blocked, B_blocked, F_blocked,thread_num,0);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for blocked parallel filtered bmm: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* C_csr_from_blocked = blocked2csr(C_blocked);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for C blocked->CSR: %.5f seconds.\n\n", elapsed);

    // Check result
    uint32_t check = check_result(filename_C, C_csr);
    if (check == 0) {
        printf("Wrong result of parallel bmm filtered.\n");
    } 
    else {
        printf("Correct result of parallel bmm filtered.\n");
    }

    check = check_result(filename_C, C_csr_from_blocked);
    if (check == 0) {
        printf("Wrong result of parallel bmm blocked filtered.\n");
    }
    else {
        printf("Correct result of parallel bmm blocked filtered.\n");
    }

    free_coo(A_coo);
    free_coo(B_coo);
    free_coo(F_coo);

    free_comp_matrix(A_csr);
    free_comp_matrix(B_csc);
    free_comp_matrix(F_csr);
    free_comp_matrix(C_csr);
    free_comp_matrix(C_csr_from_blocked);

    free_block_comp_matrix(A_blocked);
    free_block_comp_matrix(B_blocked);
    free_block_comp_matrix(F_blocked);
    free_block_comp_matrix(C_blocked);

    return 0;
}