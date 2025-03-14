#include "utils.h"
#include "devoir_1.h"
#include <time.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQUARE(a) ((a) * (a))

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



/**
 * @brief Computes eigenvalues and eigenvectors of a discrete Laplacian matrix
 * @param lx Length of the domain in the x direction
 * @param ly Length of the domain in the y direction
 * @param nx Number of inner nodes in the x direction
 * @param ny Number of inner nodes in the y direction
 */
#define MAX_SIZE_TEST 25

int main(int argc, char *argv[]) {
    double time_array[MAX_SIZE_TEST] ;
    int size_array[MAX_SIZE_TEST] ;
    int band_array[MAX_SIZE_TEST] ;

    for(int i = 0; i < MAX_SIZE_TEST; i++){
        printf("iteration %d \n", i);
        
        double lx = 6.0 + 3*i;
        double ly = 6.0 + 3*i;
        int nx = 25;
        int ny = 25;
        int n = nx * ny;
        double *A;
    
        A = create_matrix(nx, ny, lx, ly, 2);
        clock_t start = clock();
        qr_eigs_band(A, n, nx, 1e-15, 1000);
        clock_t end = clock();
        float seconds = (float)(end - start) / CLOCKS_PER_SEC;

        time_array[i] = seconds;
        band_array[i] = nx;
        size_array[i] = n;

        printf ("Your calculations took %.6lf seconds to run.\n", seconds );

        free(A);

    }
    FILE *fp = fopen("timing_results-BandVar.csv", "w");
    if (fp == NULL) {
        perror("Unable to open file");
        return 1;
    }
    fprintf(fp, "Iteration,Time (seconds),Band Size,Array Size\n");

    for (int i = 0; i < MAX_SIZE_TEST; i++) {
        fprintf(fp, "%d,%.6lf,%d,%d\n", i, time_array[i], band_array[i], size_array[i]);
    }


    fclose(fp);
}