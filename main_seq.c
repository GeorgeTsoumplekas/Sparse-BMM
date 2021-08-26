#include "bmm_seq.h"

int main(void) {
    struct timespec begin, end;
    srand(time(NULL));
    matrix_2d *a = rand_matrix_2d(5000, 5000);
    matrix_2d *b = rand_matrix_2d(5000, 5000);
    matrix_2d *f = rand_matrix_2d(3, 3);

    // print_matrix_2d(a);
    // printf("\n\n");
    // print_matrix_2d(b);
    // printf("\n\n");
    // print_matrix_2d(f);
    // printf("\n\n");

    comp_matrix *A = matrix2csr(a);
    comp_matrix *B = matrix2csc(b);
    comp_matrix *F = matrix2csc(f);

    // print_csr(A);
    // print_csc(B);
    // print_csc(F);

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &begin);

    // comp_matrix *C = bmm_filtered_seq(A, B, F);
    comp_matrix *C = bmm_seq(A, B);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds * 1e-9;

    // printf("I'm C\n");
    // print_csc(C);

    // matrix_2d *c = bmm_seq_2d(a, b);
    // comp_matrix *c_ = matrix2csc(c);

    // printf("I'm c_\n");
    // print_csc(c_);

    printf("\nTime elapsed: %.5f seconds.\n\n", elapsed);

    free_matrix_2d(a);
    free_matrix_2d(b);
    free_matrix_2d(f);
    // free_matrix_2d(c);

    free_compressed(A);
    free_compressed(B);
    free_compressed(F);
    free_compressed(C);
    // free_compressed(c_);

    return 0;
}