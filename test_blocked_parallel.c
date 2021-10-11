#include "utilities.h"
#include "blocked_bmm_parallel.h"

int main(int argc, char* argv[]){
    if(argc<6){
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

    int thread_num = atoi(argv[4]);
    printf("You have chosen %d threads.\n",thread_num);

    uint32_t b = atoi(argv[5]);
    printf("b=%d\n",b);

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
    
    block_comp_matrix* A_blocked = csr2blocked(A_csr,b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for A CSR->blocked: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* B_blocked = csc2blocked(B_csc,b);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for B CSC->blocked: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* C_blocked = blocked_bmm_parallel(A_blocked,B_blocked,thread_num);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for blocked bmm: %.5f seconds.\n", elapsed);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);
    
    comp_matrix* C_csr = blocked2csr(C_blocked);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf("\nTime elapsed for C blocked->CSR: %.5f seconds.\n", elapsed);

    //Check result
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

    free_block_comp_matrix(A_blocked);
    free_block_comp_matrix(B_blocked);
    free_block_comp_matrix(C_blocked);

    return 0;
}