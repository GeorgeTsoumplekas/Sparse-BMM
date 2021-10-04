#include "blocked_bmm_parallel_MPI.h"

int main(int argc, char* argv[]){

    if(argc<4){
        printf("Not enough arguments.\n");
        exit(-1);
    }

    struct timespec begin, end;
    long seconds;
    long nanoseconds;
    double elapsed;

    int initialized, finalized;
    int rank;

    char* filename_C = NULL;

    coo_matrix* A_coo = NULL;
    coo_matrix* B_coo = NULL;

    comp_matrix* A_csr = NULL;
    comp_matrix* B_csc = NULL;

    MPI_Initialized(&initialized);
    if(!initialized){
        MPI_Init(&argc, &argv);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank==0){
        char* filename_A = argv[1];
        char* filename_B = argv[2];
        filename_C = argv[3];

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
    }

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    comp_matrix* C_csr = nonblocked_bmm_parallel_2(A_csr,B_csc,rank);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    if(rank==0){

        printf("\nTime elapsed for bmm: %.5f seconds.\n", elapsed);

        // Check result
        uint32_t check = check_result(filename_C,C_csr);
        if(check==0){
            printf("Wrong result.\n");
        }
        else{
            printf("Correct result.\n");
        }

        free_coo(A_coo);
        free_coo(B_coo);

        free_comp_matrix(C_csr);
    }


    MPI_Finalized(&finalized);
    if(!finalized){
        MPI_Finalize();
    }
    
    return 0;
}