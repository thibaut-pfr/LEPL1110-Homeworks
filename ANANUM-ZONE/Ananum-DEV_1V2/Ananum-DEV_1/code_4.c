#include "code_4.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQUARE(a) ((a) * (a))

/**
 * @brief Fills the matrix with the (opposite of the) Laplacian operator
 * @param nx Number of inner nodes in the x direction
 * @param ny Number of inner nodes in the y direction
 * @param lx Length of the domain in the x direction
 * @param ly Length of the domain in the y direction
 * @param storage 0 if full, 1 if band, 2 if sym_band
 * @return the pointer to the discrete Laplacian matrix
 */
double *create_matrix(int nx, int ny, double lx, double ly, int storage) {

    int lda, k;
    int size = nx * ny; // Number of nodes/unknowns
    double dx2 = SQUARE(lx / (nx + 1));
    double dy2 = SQUARE(ly / (ny + 1));
    double alpha, beta, gamma;
    double *L;

    // Choice of node numbering, here nx=4, ny=5
    // j\i    0   1   2   3
    //     .  .   .   .   .  .
    // 0   .  0   1   2   3  .
    // 1   .  4   5   6   7  .
    // 2   .  8   9  10  11  .
    // 3   . 12  13  14  15  .
    // 4   . 16  17  18  19  .
    //     .  .   .   .   .  .

    k = nx;
    alpha = 1. / dx2;
    beta = 1. / dy2;
    gamma = 2 * (alpha + beta);

    if (storage == 2) { // Symmetric band storage
        lda = k + 1;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int l = 0; l < size; l++) {
            L[l * lda + k - k] = -beta; // (i,j)->(i-1,j)
            if (l % k != 0)
                L[l * lda + k - 1] = -alpha; // (i,j)->(i,j-1)
            L[l * lda + k - 0] = +gamma;     // (i,j)->(i  ,j)
        }
    } else if (storage == 1) { // Band storage
        lda = 2 * k + 1;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int l = 0; l < size; l++) {
            L[l * lda + k - k] = -beta; // (i,j)->(i-1,j)
            if (l % k != 0)
                L[l * lda + k - 1] = -alpha; // (i,j)->(i,j-1)
            L[l * lda + k + 0] = +gamma;     // (i,j)->(i  ,j)
            if (l % k != k - 1)
                L[l * lda + k + 1] = -alpha; // (i,j)->(i,j+1)
            L[l * lda + k + k] = -beta;      // (i,j)->(i+1,j)
        }
    } else { // Full storage
        lda = size;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int idx, i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                idx = i * k + j;
                L[idx * lda + idx] = gamma; // (i,j)->(i  ,j)
                if (0 < i)
                    L[idx * lda + idx - k] = -beta; // (i,j)->(i-1,j)
                if (i < ny - 1)
                    L[idx * lda + idx + k] = -beta; // (i,j)->(i+1,j)
                if (0 < j)
                    L[idx * lda + idx - 1] = -alpha; // (i,j)->(i,j-1)
                if (j < nx - 1)
                    L[idx * lda + idx + 1] = -alpha; // (i,j)->(i,j+1)
            }
        }
    }
    return L;
}

int LU_f_band(double *A, int n, int k) {
    // Perform the LU factorization in-place in A
    for (int i = 0; i < n; i++) {
        // Compute the multipliers and update the lower part of the matrix
        for (int j = i+1; j < MIN(n, i + k + 1); j++) {
            A[j * (2 * k + 1) + k - (j - i)] /= A[i * (2 * k + 1) + k];
            double c = A[j * (2 * k + 1) + k - (j - i)];
            for (int l = i + 1; l < MIN(n, i + k + 1); l++) {
                A[j * (2 * k + 1) + k - (j - l)] -= c * A[i * (2 * k + 1) + k - (i - l)];
            }
        }
    }
    return 0;
}

/**
 * @brief Computes the eigenvector/eigenvalue pair of a symmetric band matrix
 * using the inverse iteration/Rayleigh quotien
 * @param A Laplacian matrix of desired storage
 * @param n Size of the square matrix A
 * @param k Number of non-zero off-diagonals of A
 * @param mu Pointer to the shift value, contains the eigenvalue at the end
 * @param v Pointer to the eigenvector (size n)
 * @return 0 if the algorithm converged, 1 otherwise
 */
int inverse_power(double *A, int n, int k, double *mu, double *v) {

    int max_iter = 1000;
    double tol = 1e-12;
    double norm, lambda_new;
    double lambda_old = *mu;
    double *temp = (double *)malloc(n * sizeof(double));

    // compute LU factorisation of A-mu*I
    double *LU = (double *)malloc(n * (2*k+1) * sizeof(double));
    memcpy(LU, A, n * (2*k+1) * sizeof(double));

    for (int i = 0; i < n; i++) {
        LU[i * (2*k+1) + k] -= *mu;
    }

    // LU FACTO
    LU_f_band(LU, n, k);

    // Initialize eigenvector with random values
    for (int i = 0; i < n; i++) {
        v[i] = (double) rand() / RAND_MAX;
    }

    // Compute the norm of v
    norm = cblas_dnrm2(n, v, 1);
    cblas_dscal(n, 1.0 / norm, v, 1);

    for (int iter = 0; iter < max_iter; iter++) {
        // Solve Ay = z using LU factorisation knowing that A = LU is stored in band format with k non-zero off-diagonals
        // forward substitution
        for (int i = 0; i < n; i++){
            double sum = v[i];
            for (int j = MAX(0, i - k); j < i; j++){
                sum -= LU[(j - i) + i * (2 * k + 1) + k] * v[j];
            }
            v[i] = sum;
        }

        // backward substitution
        for (int i = n-1; i >= 0; i--){
            double sum = v[i];
            for (int j = i+1; j < MIN(n, i + k + 1); j++){
                sum -= LU[(j - i) + i * (2 * k + 1) + k] * v[j];
            }
            v[i] = sum/LU[(i) * (2 * k + 1) + k];
        }

        // Compute the norm of v
        norm = cblas_dnrm2(n, v, 1);
        cblas_dscal(n, 1.0 / norm, v, 1);
        
        // Rayleigh quotient
        // Perform the matrix-vector multiplication A * v
        cblas_dgbmv(CblasRowMajor, CblasNoTrans, n, n, k, k, 1.0, A, 2 * k + 1, v, 1, 0.0, temp, 1);
        // Compute the dot product v^T * (A * v)
        lambda_new = cblas_ddot(n, v, 1, temp, 1);

        // Check for convergence
        if (fabs(lambda_new - lambda_old) < tol) {
            *mu = lambda_new;
            free(LU); free(temp);
            return 0;
        }
        lambda_old = lambda_new;
    }
    free(LU);
    free(temp);
    return 1;
}

int inverse_power_steroid(double *A, double *LU, int n, int k, double *mu, double *v, double *eigws, double *eigvs, int idx) {
    int max_iter = 1000;
    double tol = 1e-12;
    double norm, lambda_new;
    double lambda_old = *mu;
    double *temp = (double *)malloc(n * sizeof(double));

    // Initialize eigenvector with random values
    for (int i = 0; i < n; i++) {
        v[i] = (double) rand() / RAND_MAX;
    }

    // Compute the norm of v
    norm = cblas_dnrm2(n, v, 1);
    cblas_dscal(n, 1.0 / norm, v, 1);

    for (int iter = 0; iter < max_iter; iter++) {
        // Solve Ay = z using LU factorisation knowing that A = LU is stored in band format with k non-zero off-diagonals

        /*
        // forward substitution
        for (int i = 0; i < n; i++){
            double sum = v[i];
            for (int j = MAX(0, i - k); j < i; j++){
                sum -= LU[(j - i) + i * (2 * k + 1) + k] * v[j];
            }
            v[i] = sum;
        }

        // backward substitution
        for (int i = n-1; i >= 0; i--){
            double sum = v[i];
            for (int j = i+1; j < MIN(n, i + k + 1); j++){
                sum -= LU[(j - i) + i * (2 * k + 1) + k] * v[j];
            }
            v[i] = sum/LU[(i) * (2 * k + 1) + k];
        }
        */
        cblas_dtbsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit, n, k, LU, 2 * k + 1, v, 1);
        cblas_dtbsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, k, LU+k, 2 * k + 1, v, 1);

        // substract the shift Σ v_i v_i^T v to w
        for(int k = 0; k < idx; k++){
            double viTv = cblas_ddot(n, &eigvs[k * n], 1, v, 1);
            cblas_daxpy(n, -viTv, &eigvs[k * n], 1, v, 1);
        }

        // Compute the norm of v
        norm = cblas_dnrm2(n, v, 1);
        cblas_dscal(n, 1.0 / norm, v, 1);
        
        // Rayleigh quotient
        // Perform the matrix-vector multiplication A * v
        cblas_dgbmv(CblasRowMajor, CblasNoTrans, n, n, k, k, 1.0, A, 2 * k + 1, v, 1, 0.0, temp, 1);
        // Compute the dot product v^T * (A * v)
        lambda_new = cblas_ddot(n, v, 1, temp, 1);
        lambda_new = 1/lambda_new;

        // substract the shift Σ lambda_i (v_i^T v)²
        for(int k = 0; k < idx; k++){
            double viTv = cblas_ddot(n, &eigvs[k * n], 1, v, 1);
            lambda_new -= eigws[k] * SQUARE(viTv);
        }

        // Check for convergence
        if ((fabs(lambda_new - lambda_old) < tol) && (iter > 0)) {
            *mu = lambda_new;
            free(temp);
            return 0;
        }
        lambda_old = lambda_new;
    }
    free(temp);
    return 1;
}

/**
 * @brief Computes the nb smallest eigenvalues of a symmetric band matrix
 * using the deflation technique
 * @param L Laplacian matrix of desired storage
 * @param n Size of the square matrix A
 * @param k Number of non-zero off-diagonals of A
 * @param nb Number of eigenvalues to compute
 * @param eigws Array to store the eigenvalues (size nb)
 * @param eigvs Array to store the eigenvectors (size nb x n, row-major)
 * @return 0 if the algorithm converged, 1 otherwise
 */
int deflation(double *L, int n, int k, int nb, double *eigws, double *eigvs) {
    double *v = (double *)malloc(n * sizeof(double));
    double *LU = (double *)malloc(n * (2*k+1) * sizeof(double));
    // LU FACTO
    memcpy(LU, L, n * (2*k+1) * sizeof(double));
    
    LU_f_band(LU, n, k);
    
    for (int i = 0; i < nb; i++) {
        double mu = 0.0;
        int result = inverse_power_steroid(L, LU, n, k, &mu, v, eigws, eigvs, i);
        if (result == 0) {
            eigws[i] = mu;
            memcpy(&eigvs[i * n], v, n * sizeof(double));
        } else {
            free(v);
            free(LU);
            return 1;
        }
    }
    // inverse eigenvalues
    for (int i = 0; i < nb; i++) {
        eigws[i] = 1/eigws[i];
    }
    free(v);
    free(LU);
    return 0;
}

/**
 * @brief Indicates the storage format used
 * @return 0 if full, 1 if band, 2 if sym_band
 */
int get_storage_format() { return 1; }
