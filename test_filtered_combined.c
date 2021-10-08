#include "combined_blocked_bmm_parallel.h"

int main(int argc, char* argv[]){

    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    int initialized, finalized;
    int rank;

    char* filename_C = NULL;

    coo_matrix* A_coo = NULL;
    coo_matrix* B_coo = NULL;
    coo_matrix* F_coo = NULL;

    comp_matrix* A_csr = NULL;
    comp_matrix* B_csc = NULL;
    comp_matrix* F_csr = NULL;

    block_comp_matrix* A_blocked = NULL;
    block_comp_matrix* B_blocked = NULL;
    block_comp_matrix* F_blocked = NULL;

    int thread_num;

    MPI_Initialized(&initialized);
    if(!initialized){
        MPI_Init(&argc, &argv);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank==0){
        printf("\nThis file completes the filtered bmm with both MPI and OpenMP parallelization with simple and blocked matrices.\n");

        if (argc < 7) {
            printf("Not enough arguments.\n");
            printf("Usage: mpirun -n [no. processes] %s [txt file of matrix A] [txt file of matrix B] ",argv[0]);
            printf("[txt file of matrix F] [txt file of matrix C] [b] [no. threads for OpenMP]\n");
                    
            exit(-1);
        }

    }

    thread_num = atoi(argv[5]);
    if(rank==0){
        printf("You have chosen %d threads.\n",thread_num);
    }

    if(rank==0){
        

        char* filename_A = argv[1];
        char* filename_B = argv[2];
        char* filename_F = argv[3];
        filename_C = argv[4];
        
        uint32_t b = atoi(argv[6]);
        printf("b=%d\n", b);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        A_coo = readCOO(filename_A);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for A txt->COO: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        B_coo = readCOO(filename_B);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for B txt->COO: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        F_coo = readCOO(filename_F);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for F txt->COO: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        A_csr = coo2csr(A_coo);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for A COO->CSR: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        B_csc = coo2csc(B_coo);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for B COO->CSC: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        F_csr = coo2csr(F_coo);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for F COO->CSR: %.5f seconds.\n", elapsed);

            // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        A_blocked = csr2blocked(A_csr, b);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for A CSR->blocked: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        B_blocked = csc2blocked(B_csc, b);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for B CSC->blocked: %.5f seconds.\n", elapsed);

        // Start timer
        clock_gettime(CLOCK_MONOTONIC, &begin);

        F_blocked = csr2blocked(F_csr, b);

        // End timer
        clock_gettime(CLOCK_MONOTONIC, &end);
        seconds = end.tv_sec - begin.tv_sec;
        nanoseconds = end.tv_nsec - begin.tv_nsec;
        elapsed = seconds + nanoseconds * 1e-9;

        printf("\nTime elapsed for F CSR->blocked: %.5f seconds.\n", elapsed);
    }

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* C_csr = combined_nonblocked_bmm_parallel_filtered(A_csr,B_csc,F_csr,rank,10);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank==0){

        printf("\nTime elapsed for bmm filtered: %.5f seconds.\n", elapsed);

        // Check result
        uint32_t check = check_result(filename_C,C_csr);
        if(check==0){
            printf("Wrong result.\n");
        }
        else{
            printf("Correct result.\n");
        }
    }

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    block_comp_matrix* C_blocked = combined_blocked_bmm_parallel_filtered(A_blocked,B_blocked,F_blocked,rank,10);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank==0){

        printf("\nTime elapsed for blocked bmm filtered: %.5f seconds.\n", elapsed);

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
        uint32_t check = check_result(filename_C,C_csr_from_blocked);
        if(check==0){
            printf("Wrong result.\n");
        }
        else{
            printf("Correct result.\n");
        }

        free_coo(A_coo);
        free_coo(B_coo);
        free_coo(F_coo);

        free_comp_matrix(C_csr);
        free_block_comp_matrix(C_blocked);
        free_comp_matrix(C_csr_from_blocked);
    }


    MPI_Finalized(&finalized);
    if(!finalized){
        MPI_Finalize();
    }
    
    return 0;
}