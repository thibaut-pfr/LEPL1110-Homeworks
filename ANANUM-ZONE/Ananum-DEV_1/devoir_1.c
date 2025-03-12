#include "utils.h"

#define get(A, i, j) (A[((j)-(i)) + (i) * ((k) + 1) + (k)])
#define get_safe(A, i, j) (((i) < 0 || (i) >= n || (j) < 0 || (j) >= n || abs((i) - (j)) > k) ? 0 : ((j) < (i) ? A[((j)-(i)) + ((i)) * ((k) + 1) + (k)] : A[((i)-(j)) + (j) * ((k) + 1) + (k)]))

#define set(A, i, j, val) (A[((j)-(i))+ (i) * ((k) + 1) + (k)] = (val))
#define set_safe(A, i, j, val) \
    if (!((i) < 0 || (i) >= n || (j) < 0 || (j) >= n || abs((i) - (j)) > k || (j) > (i))) { \
        A[((j)-(i)) + (i) * ((k) + 1) + (k)] = (val); \
    }

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))


double givens(double* A, int n, int k, int i, int j){
    // Rotation coefficients
    double a = get(A,i-1,j);
    double b = get(A,i,j);
    double r = sqrt(a * a + b * b);
    double c = a / r;
    double s = b / r;

    // Apply rotation to rows i and i-1
    int start = MAX(0,i-k);
    int stop = i;
    double temp_i = 0;
    double temp_j = 0;
    double temp_diag =  get(A, stop, i-1); // idem que get_safe(A,i-1,stop) par symétrie;

    for (int l = start; l < stop; l++) { 
        temp_i = c * get(A,i-1,l) + s * get(A,i,l);
        temp_j = - s * get(A,i-1,l) + c * get(A,i,l);
        set(A,i-1,l, temp_i);
        set(A,i,l, temp_j);
    }

    // stop iteration
    // We also need a way to remember the value of the last element of the row i-1 at the end of the loop 
    // that falls in the upper band of the matrix
    // Because we need it for the start of the column rotation
    // (as the matrix is in sym-band format, the value of the last element of the row i-1 is not stored in the matrix)
    // We need to store it in a temporary variable

    temp_i = c * temp_diag + s * get(A,i,stop);
    temp_j = - s * temp_diag + c * get(A,i,stop);
    temp_diag = temp_i;
    set(A,i,stop, temp_j);

    // Apply rotation to columns i and i-1
    start = i-1;
    stop = MIN(i+k, n-1);
    double temp_garbage = 0;

    // As previously said, we use the temp value stored at the first iteration of the loop
    // the other iterations does not suffer from this as their values are stored in the matrix

    // start iteration
    temp_i =  c * get(A,start,i-1) + s * temp_diag;
    temp_j = -s * get(A,start,i-1) + c * temp_diag;
    set(A,start,i-1, temp_i);

    // main iterations
    for (int l = start+1; l < stop; l++) {
        temp_i =  c * get(A,l,i-1) + s * get(A,l,i);
        temp_j = - s * get(A,l,i-1) + c * get(A,l,i);   
        set(A,l,i-1, temp_i);
        set(A,l,i, temp_j);
    }

    // stop iteration
    temp_i = c * get_safe(A,stop,i-1) + s * get(A,stop,i);   // Todo : if  stop == i+k, then get(A,stop,i-1) = 0
    temp_j = - s * get_safe(A,stop,i-1) + c * get(A,stop,i); // Todo ??
    
    temp_garbage = temp_i;                                   // Todo : This is needed if stop != n-1
    set_safe(A,stop,i-1, temp_i);                            // Todo : this is needed if stop = n-1
    set(A,stop,i, temp_j);

    return temp_garbage;
}


/*
** d²u/dx² = λu ~~ (Ui-1 - 2Ui + Ui+1)/h² = λUi
** Une fonction qui tridiagonalise une matrice symétrique bande par transformations de similitudes.
** Vous pouvez choisir le format de stockage et le type de transformations (réflexions de Householder
** ou rotations de Givens).: 
• A contient les entrées de la matrice dans un tableau Row-Major
– de taille n × n si [format] est full, ou
– de taille n × (k + 1) si [format] est band.
• n indique la taille de la matrice
• k indique le nombre de sous-diagonales non-nulles
• d est un tableau de taille n qui contient en sortie la diagonale principale de la matrice
tridiagonale (Hessenberg symétrique)
• e est un tableau de taille n qui contient en sortie la sous-diagonale de la matrice tridiagonale
dans ses n − 1 derniers éléments
*/
int tridiagonalize_band(double *A, int n, int k, double *d, double *e){
    // If matrix is already tridiagonal, copy then return
    if (k == 0){
        for(int i = 0; i < n; i++){
            d[i] = A[i * n + i];
        }
        return 0;
    }
    else if(k == 1){
        for(int i = 0; i < n; i++){
            d[i] = A[k + (k+1) * i ];
        }
        for(int i = 1; i < n; i++){
            e[i] = A[(i) * (k + 1) + k - 1];
        }
        return 0;
    }

    // Apply Givens rotations to zero out the sub-diagonals
    for (int j = 0; j < n - 2; j++) {
        for (int i = MIN(j + k, n-1); i > j + 1; i--) {
            if(get(A,i,j) == 0) continue;
            double temp_i, temp_j, a, b, r, c, s, temp_diag;
            int start, stop;
            double temp_garbage = givens(A, n, k, i, j); // just one givens rotation to get one zero

            // temp_garbage contains the garbage value below the band ! 
            // BAND - GIVENS ROTATION : GET RID OF PARASITIC ELEMENT
            // Be careful that when we make one zero appears in the matrix at the position (i,j)
            // The entry at position (i+k,i-1) which is out of the band becomes non zero
            // We will do a rotation to put it to zero, wich will make another element non zero to appear
            // Farther in the matrix, and so on... Until we reach the end of the matrix 

            int i_chase = i+k; int j_chase = i-1;

            while((i_chase < n)){
                // Rotation coefficients
                a = get(A,i_chase-1,j_chase); 
                b = temp_garbage;
                r = sqrt(a * a + b * b);
                c = a / r;
                s = b / r;
                start = i_chase-k-1;
                stop = i_chase;

                temp_diag =  get(A, stop, i_chase-1);

                // start iteration
                temp_i = c * get(A,i_chase-1,start) + s * temp_garbage;
                temp_j = - s * get(A,i_chase-1,start) + c * temp_garbage;
                set(A,i_chase-1,start, temp_i);
                temp_garbage = 0;

                for (int l = start+1; l < stop; l++) { 
                    temp_i = c * get(A,i_chase-1,l) + s * get(A,i_chase,l);
                    temp_j = - s * get(A,i_chase-1,l) + c * get(A,i_chase,l);
                    set(A,i_chase-1,l, temp_i);
                    set(A,i_chase,l, temp_j);
                }

                // stop iteration
                temp_i = c * temp_diag + s * get(A,i_chase,stop);
                temp_j = - s * temp_diag + c * get(A,i_chase,stop);
                temp_diag = temp_i;
                set(A,i_chase,stop, temp_j);

                // Apply rotation to columns i and i-1
                start = i_chase-1;
                stop = MIN(i_chase+k, n-1);
                
                // start iteration
                temp_i =  c * get(A,start,i_chase-1) + s * temp_diag;
                temp_j = -s * get(A,start,i_chase-1) + c * temp_diag;
                set(A,start,i_chase-1, temp_i);
    
                for (int l = start+1; l < stop; l++) {
                    temp_i = c * get(A,l,i_chase-1) + s * get(A,l,i_chase);
                    temp_j = - s * get(A,l,i_chase-1) + c * get(A,l,i_chase);
                    set(A,l,i_chase-1, temp_i);
                    set(A,l,i_chase, temp_j);
                }
                // stop iteration
                temp_i = c * get_safe(A,stop,i_chase-1) + s * get(A,stop,i_chase);   // Todo : if  stop == i_chase + k, then get(A,stop,i-1) = 0
                temp_j = - s * get_safe(A,stop,i_chase-1) + c * get(A,stop,i_chase); // Todo : ...
                temp_garbage = temp_i;
                set(A,stop,i_chase, temp_j);
                set_safe(A,stop,i_chase-1, temp_i);                                 // Todo : this is needed if stop = n-1

                int t = i_chase;
                i_chase = t+k;
                j_chase = t-1;
            }
        }
    }
    
    // Extract the diagonal and sub-diagonal elements
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (k + 1) + k];
    }
    for (int i = 1; i < n; i++) {
        e[i] = A[(i) * (k + 1) + k - 1];
    }
    return 0;
    
}
/*
Une fonction qui effectue une étape de l’algorithme QR avec un shift de Wilkinson μ sur une
matrice tridiagonale symétrique.
int step_qr_tridiag(double *d, double *e, int m, double eps);
• d et e sont définis comme précédemment, et sont mis à jour en sortie de telle sorte que
Hk+1 ← Q∗
kHkQk
• m indique la taille de la matrice active (m ≤ n), c’est-à-dire sur laquelle on n’a pas encore
isolé les valeurs propres
• eps est la tolérance pour déclarer qu’une valeur propre a été isolée : Golub et Van Loan
proposent |hp+1,p| ≤ ϵ(|hpp| + |hp+1,p+1|).
• En sortie, la fonction retourne la nouvelle dimension active ≤ m
*/
int step_qr_tridiag(double *d, double *e, int m, double eps){
    // 1. compute the Wilikinson shift

    // 2. Substract the shift to the diagonal

    // 3. Apply the Givens LEFT rotations to the sub-diagonal
    // Be careful that we will need to store the rotations performed to apply it after all the LEFT rotations

    // 4. Apply the Givens RIGHT rotations

    // 5. Check if the sub-diagonal is small enough to be considered as zero
    // If yes, reduce m by 1

    // 6. Add the shift back to the diagonal

    
    // For tests : compute_eigenvalues(A, n, k);
    return m;
}

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
    tridiagonalize_band(A, 10, 6, d, e);

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


int main() {
    return 0;
}
