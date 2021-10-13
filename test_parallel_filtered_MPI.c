#include "blocked_bmm_parallel_MPI.h"

int main(int argc, char* argv[]){

    if(argc<5){
        printf("Not enough arguments.\n");
        exit(-1);
    }

    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    int initialized, finalized;
    int rank, numtasks;

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

    MPI_Initialized(&initialized);
    if(!initialized){
        MPI_Init(&argc, &argv);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    if(rank==0){
        printf("\nThis file completes the filtered bmm with MPI parallelization with simple and blocked matrices.\n");

        if (argc < 6) {
            printf("Not enough arguments.\n");
            printf("Usage: mpirun -n [no. processes] %s [txt file of matrix A] [txt file of matrix B] ",argv[0]);
            printf("[txt file of matrix F] [txt file of matrix C] [b]\n");
                
            exit(-1);
        }

        printf("Number of processes: %d\n",numtasks);

        char* filename_A = argv[1];
        char* filename_B = argv[2];
        char* filename_F = argv[3];
        filename_C = argv[4];
        
        uint32_t b = atoi(argv[5]);
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

        //It is necessary that b * # of processes less or equal than n
        if(b*numtasks>A_coo->n){
            printf("ERROR: Condition that b * # of processes <= n is not held. Please try again with differend parameters.\n");
            exit(-1);
        }

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

    comp_matrix* C_csr = nonblocked_bmm_parallel_filtered_2(A_csr,B_csc,F_csr,rank);

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

    block_comp_matrix* C_blocked = blocked_bmm_parallel_filtered_2(A_blocked,B_blocked,F_blocked,rank);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank==0){

        printf("\nTime elapsed for blocked bmm filtered: %.5f seconds.\n", elapsed);

        free_coo(A_coo);
        free_coo(B_coo);
        free_coo(F_coo);

        free_comp_matrix(C_csr);
        free_block_comp_matrix(C_blocked);
    }


    MPI_Finalized(&finalized);
    if(!finalized){
        MPI_Finalize();
    }
    
    return 0;
}