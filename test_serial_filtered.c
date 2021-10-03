#include "bmm_seq.h"
#include "utilities.h"

int main(int argc, char* argv[]) {
    printf(
        "This file completes the filtered bmm with the serial way and the "
        "serial blocked way\n");

    if (argc < 6) {
        printf("Not enough arguments.\n");
        printf(
            "Usage: %s [txt file of matrix A] [txt file of matrix B] [txt file "
            "of matrix F] [txt file of matrix C] [b]\n",
            argv[0]);
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

    uint32_t b = atoi(argv[5]);
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

    printf("\nTime elapsed for B txt->COO: %.5f seconds.\n", elapsed);

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

    comp_matrix* C_csr = bmm_filtered_seq(A_csr, B_csc, F_csr, 0);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for filtered bmm: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* C_blocked = blocked_bmm_seq_filtered(A_blocked, B_blocked, F_blocked);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for blocked filtered bmm: %.5f seconds.\n", elapsed);

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
        printf("Wrong result of bmm filtered.\n");
    } 
    else {
        printf("Correct result of bmm filtered.\n");
    }

    check = check_result(filename_C, C_csr_from_blocked);
    if (check == 0) {
        printf("Wrong result of bmm blocked filtered.\n");
    }
    else {
        printf("Correct result of bmm blocked filtered.\n");
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