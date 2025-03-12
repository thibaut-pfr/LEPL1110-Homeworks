#include "utils.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define get(A, i, j) (A[((j)-(i)) + (i) * ((k) + 1) + (k)])
#define get_safe(A, i, j) (((i) < 0 || (i) >= n || (j) < 0 || (j) >= n || abs((i) - (j)) > k) ? 0 : ((j) < (i) ? A[((j)-(i)) + ((i)) * ((k) + 1) + (k)] : A[((i)-(j)) + (j) * ((k) + 1) + (k)]))

#define set(A, i, j, val) (A[((j)-(i))+ (i) * ((k) + 1) + (k)] = (val))
#define set_safe(A, i, j, val) \
    if (!((i) < 0 || (i) >= n || (j) < 0 || (j) >= n || abs((i) - (j)) > k || (j) > (i))) { \
        A[((j)-(i)) + (i) * ((k) + 1) + (k)] = (val); \
    }

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

void compute_eigenvalues(double *A, int n, int k) {
    // Convert band symmetric matrix to full format
    double *full_A = (double *)malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            full_A[i * n + j] = get_safe(A, i, j);
        }
    }
    // Compute the eigenvalues using LAPACKE_dsyev
    char jobz = 'N'; // Compute eigenvalues only
    char uplo = 'U'; // Upper triangle of A is stored
    int lda = n;
    double *w = (double *)malloc(n * sizeof(double));
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, full_A, lda, w);
    // Sort the eigenvalues
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (w[i] > w[j]) {
                double temp = w[i];
                w[i] = w[j];
                w[j] = temp;
            }
        }
    }
    if (info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
    } else {
        printf("Eigenvalues = [");
        for (int i = 0; i < n; i++) {
            printf("%f", w[i]);
            if (i < n - 1) printf(", ");
        }
        printf("]\n");
    }

    free(full_A);
    free(w);
}