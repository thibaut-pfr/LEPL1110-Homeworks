#include "code_4.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define TOL_PRINT 0e-15
#define MAX_SIZE_PRINT 25

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

void save_eigv(
    char *name, double *eigv, double eigw, double lx, double ly, int nx, int ny
) {
    FILE *f = fopen(name, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", name);
        exit(1);
    }
    fprintf(f, "# %le %le %le %d %d\n", eigw, lx, ly, nx, ny);
    for (int i = 0; i < nx * ny; i++) {
        fprintf(f, "%20.15le\n", eigv[i]);
    }
    fclose(f);
}


/**
 * @brief Computes eigenvalues and eigenvectors of a discrete Laplacian matrix
 * @param lx Length of the domain in the x direction
 * @param ly Length of the domain in the y direction
 * @param nx Number of inner nodes in the x direction
 * @param ny Number of inner nodes in the y direction
 */
int main(int argc, char *argv[]) {

    double lx = 5.0;
    double ly = 5.0;
    int nx = 2;
    int ny = 2;
    int n = nx * ny;

    double *A;

    A = create_matrix(nx, ny, lx, ly, 0);
//    print_full(A, n, "A");
    free(A);

    A = create_matrix(nx, ny, lx, ly, 1);
//    print_band(A, n, nx, "L");
    free(A);

    A = create_matrix(nx, ny, lx, ly, 2);
//   print_sym_band(A, n, nx, "L");
    free(A);

    double *x = (double *)calloc(n, sizeof(double));
    for(int i = 0; i < n; i++) {
        x[i] = 1.0;
    }
    save_eigv("./eigv.txt", x, 1.0, lx, ly, nx, ny);
    free(x);

    // Create a positive definite diagonal matrix
//    double *D = (double *)calloc(n * n, sizeof(double));
/*    for (int i = 0; i < n; i++) {
        D[i * n + i] = (double)(5.5 - i);
    }
    D[(n-1)*n + n-1] = 5.0;
    D[1] = 1.0; D[n] = 1.0;
    D[6] = 3.0; D[9] = 3.0;
    print_full(D, n, "D");
*/

// Create a positive definite tridiagonal matrix stored in band format
double *D = (double *)calloc(n * 3, sizeof(double));
for (int i = 0; i < n; i++) {
    D[i * 3 + 1] = 2.0; // Main diagonal
    if (i > 0) {
        D[i * 3] = -1.0; // Sub-diagonal
    }
    if (i < n - 1) {
        D[i * 3 + 2] = -1.0; // Super-diagonal
    }
    D[1] = 3.0;
}
print_band(D, n, 1, "D");

    // Test the inverse_interation from code_4.c
    double mu = 0.0;
    double *v = (double *)calloc(n, sizeof(double));
    // inverse_power(D, n, 0, &mu, v);
    // printf("Eigenvalue: %lf\n", mu);
    free(v);

    // Test the deflation function from code_4.c
    double *eigvs = (double *)calloc(n * n, sizeof(double));
    double *eigws = (double *)calloc(n, sizeof(double));

    deflation(D, n, 1, n, eigws, eigvs);
    // print eigws
    printf("Eigenvalues\n");
    for (int i = 0; i < n; i++) {
        printf("%6.2lf ", eigws[i]);
    }
    free(eigvs);
    free(eigws);
    free(D);
    return 0;
}
