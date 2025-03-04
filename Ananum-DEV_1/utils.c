#include "utils.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))


void print_sym_band(double *A, int n, int b, char *name) {
    if (n > MAX_SIZE_PRINT)
        return;
    printf("\nSymmetric band matrix %s\n", name);
    int idx;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i - b; j++) {
            printf("%6s ", "");
        }
        for (int j = MAX(0, b - i); j < b + 1; j++) {
            idx = i * (b + 1) + j;
            // printf("%6d ", idx);
            if (fabs(A[idx]) < TOL_PRINT)
                printf("%6s ", "");
            else
                printf("%6.2lf ", A[idx]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_band(double *A, int n, int k, char *name) {
    if (n > MAX_SIZE_PRINT)
        return;
    printf("\nBand matrix %s\n", name);
    int idx;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i - k; j++) {
            printf("%6s ", "");
        }
        for (int j = MAX(0, k - i); j < MIN(2 * k + 1, k + n - i); j++) {
            idx = i * (2 * k + 1) + j;
            // printf("%6d ", idx);
            if (fabs(A[idx]) < TOL_PRINT)
                printf("%6s ", "");
            else
                printf("%6.2lf ", A[idx]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_full(double *A, int n, char *name) {
    if (n > MAX_SIZE_PRINT)
        return;
    printf("\nFull matrix %s\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // printf("%6d ", i * n + j);
            if (fabs(A[i * n + j]) < TOL_PRINT)
                printf("%6s ", "");
            else
                printf("%6.2lf ", A[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}