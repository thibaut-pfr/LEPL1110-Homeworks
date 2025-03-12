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

void print_tridiag(double *d, double *e, int n, char *name){
    printf("\nTridiagonal matrix %s\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                printf("%6.2lf ", d[i]);
            } else if (i == j - 1) {
                printf("%6.2lf ", e[i+1]);
            } else if (i == j + 1) {
                printf("%6.2lf ", e[j+1]);
            } else {
                printf("%6s ", "");
            }
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


// ===================================================== //
// ================== TEST FUNCTIONS ================== //
// ===================================================== //

void test_tridiagonalize_band_large() {
    // Define a 12x12 symmetric band matrix with band width 4
    int k = 4;
    int n = 12;
    double A[72] = {
        0, 0, 0, 0, 12,
           0, 0, 0, 12, 11,
              0, 0, 12, 11, 10,
                 0, 12, 11, 10, 9,
                    12, 11, 10, 9, 8,
                       11, 10, 9, 8, 7,
                          10, 9, 8, 7, 6,
                             9, 8, 7, 6, 5,
                                8, 7, 6, 5, 4,
                                   7, 6, 5, 4, 3,
                                      6, 5, 4, 3, 2,
                                         5, 4, 3, 2, 1
    };
    double d[12] = {0}, e[12] = {0};

    // Print the original matrix
    print_sym_band(A, 12, 4, "Original A");

    compute_eigenvalues(A, 12, 4); // True eigenvalues

    // Call the tridiagonalize_band function
    tridiagonalize_band(A, 12, 4, d, e);

    // Print the tridiagonalized matrix
    print_sym_band(A, 12, 4, "Tridiagonalized A");

    compute_eigenvalues(A, 12, 4); // Eigenvalues of the tridiagonalized matrix

    // Print the diagonal and sub-diagonal elements
    printf("d = [");
    for (int i = 0; i < 12; i++) {
        printf("%f", d[i]);
        if (i < 11) printf(", ");
    }
    printf("]\n");

    printf("e = [");
    for (int i = 0; i < 12; i++) {
        printf("%f", e[i]);
        if (i < 11) printf(", ");
    }
    printf("]\n");

    // Write the Hessenberg matrix in full format to a txt file
    FILE *file = fopen("hessenberg_matrix_large.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (abs(i - j) <= 1) {
                fprintf(file, "%f ", get_safe(A, i, j));
            } else {
                fprintf(file, "0.000000 ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}


void test_step_qr_tridiag() {
    // Define a symmetric tridiagonal matrix using vectors d and e
    int m = 5;
    double d[5] = {4.0, 3.0, 2.0, 1.0, 0.5};
    double e[5] = {0.0, 1.0, 1.0, 1.0, 1.0}; // e[0] is not used

    double eps = 1e-6;

    // Print the original tridiagonal matrix
    print_tridiag(d, e, m, "Original tridiagonal matrix");

    double A[12] = {
              0, 4,
                 1, 3,
                    1, 2,
                       1, 1,
                          1, 0.5
                        };   

    // Call the compute_eigenvalues function
    compute_eigenvalues(A, m, 1);

    int new_m = 5;
    // Perform one QR step with Wilkinson shift
    for(int i = 0; i < 10; i++){
        new_m = step_qr_tridiag(d, e, new_m, eps);
    }

    // Print the updated tridiagonal matrix
    print_tridiag(d, e, new_m, "Updated tridiagonal matrix after one QR step");

    // Print the new active size of the matrix
    printf("New active size of the matrix: %d\n", new_m);
}

void test_qr_band(){
    // Define a 12x12 symmetric band matrix with band width 4
    int k = 4;
    int n = 12;
    double A[72] = {
        0, 0, 0, 0, 12,
           0, 0, 0, 12, 11,
              0, 0, 12, 11, 10,
                 0, 12, 11, 10, 9,
                    12, 11, 10, 9, 8,
                       11, 10, 9, 8, 7,
                          10, 9, 8, 7, 6,
                             9, 8, 7, 6, 5,
                                8, 7, 6, 5, 4,
                                   7, 6, 5, 4, 3,
                                      6, 5, 4, 3, 2,
                                         5, 4, 3, 2, 1
    };

    // Print the original matrix
    print_sym_band(A, n, k, "Original A");

    // Call the qr_eigs_band function
    int max_iter = 100;
    double eps = 1e-6;
    int iter = qr_eigs_band(A, n, k, eps, max_iter);

    // Print the matrix
    print_sym_band(A, n, k, "Final A");

    // Print the number of iterations
    printf("Number of iterations: %d\n", iter);
}

void test_qr_band_small(){
    // Define a 12x12 symmetric band matrix with band width 4
    int k = 2;
    int n = 3;
    double A[9] = {0, 0, 1, 
                       0, 2, 1, 
                          3, 4, 1};

    // Print the original matrix
    print_sym_band(A, n, k, "Original A");

    // Call the qr_eigs_band function
    int max_iter = 100;
    double eps = 1e-6;
    int iter = qr_eigs_band(A, n, k, eps, max_iter);

    // Print the matrix
    print_sym_band(A, n, k, "Final A");

    // Print the number of iterations
    printf("Number of iterations: %d\n", iter);
}


void test_tridiagonalize_band() {
    // Define a 10x10 symmetric band matrix with band width 7
    int k = 6;
    int n = 10;
    double A[70] = {
        0, 0, 0, 0, 0, 0, 10,
           0, 0, 0, 0, 0, 10, 9,
              0, 0, 0, 0, 10, 9, 8,
                 0, 0, 0, 10, 9, 8, 7,
                    0, 0, 10, 9, 8, 7, 6,
                       0, 10, 9, 8, 7, 6, 5,
                          10, 9, 8, 7, 6, 5, 4,
                              9, 8, 7, 6, 5, 4, 3,
                                 8, 7, 6, 5, 4, 3, 2,
                                    7, 6, 5, 4, 3, 2, 1
    };
    double d[10] = {0}, e[10] = {0};

    // Print the original matrix
    print_sym_band(A, 10, 6, "Original A");

    // print eigenvalues
    compute_eigenvalues(A, 10, 6);

    // Call the tridiagonalize_band function
    //tridiagonalize_band(A, 10, 6, d, e);

    // Print the tridiagonalized matrix
    print_sym_band(A, 10, 6, "Tridiagonalized A");

    compute_eigenvalues(A, 10, 6);

    // Print the diagonal and sub-diagonal elements
    printf("d = [");
    for (int i = 0; i < 10; i++) {
        printf("%f", d[i]);
        if (i < 9) printf(", ");
    }
    printf("]\n");

    printf("e = [");
    for (int i = 0; i < 10; i++) {
        printf("%f", e[i]);
        if (i < 9) printf(", ");
    }
    printf("]\n");


    // Write the Hessenberg matrix in full format to a txt file
    FILE *file = fopen("hessenberg_matrix.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            if (abs(i - j) <= 1) {
                fprintf(file, "%f ", get_safe(A, i, j));
            } else {
                fprintf(file, "0.000000 ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);

    compute_eigenvalues(A, 10, 6);
}

void test_tridiagonalize_band_small() {
    // Define a 4x4 symmetric band matrix with band width 2
    int k = 2;
    int n = 4;
    double A[12] = {
        0, 0, 8,
           0, 1, 4,
              4, 3, 1,
                 1, 2, 9
    };
    double d[4] = {0}, e[4] = {0};

    // Print the original matrix
    print_sym_band(A, 4, 2, "Original A");

    // Call the tridiagonalize_band function
    tridiagonalize_band(A, 4, 2, d, e);

    // Print the tridiagonalized matrix
    print_sym_band(A, 4, 2, "Tridiagonalized A");

    // Print the diagonal and sub-diagonal elements
    printf("d = [");
    for (int i = 0; i < 4; i++) {
        printf("%f", d[i]);
        if (i < 3) printf(", ");
    }
    printf("]\n");

    printf("e = [");
    for (int i = 0; i < 4; i++) {
        printf("%f", e[i]);
        if (i < 3) printf(", ");
    }
    printf("]\n");

    // Write the Hessenberg matrix in full format to a txt file
    FILE *file = fopen("hessenberg_matrix_small.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (abs(i - j) <= 1) {
                fprintf(file, "%f ", get_safe(A, i, j));
            } else {
                fprintf(file, "0.000000 ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

