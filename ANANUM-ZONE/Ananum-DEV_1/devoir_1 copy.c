/// BEFORE GARBAGE ROTATION

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

/*
** d²u/dx² = λu ~~ (Ui-1 - 2Ui + Ui+1)/h² = λUi
** Une fonction qui tridiagonalise une matrice symétrique bande par transformations de similitudes.
** Vous pouvez choisir le format de stockage et le type de transformations (réflexions de Householder
** ou rotations de Givens).: 
**  
• A contient les entrées de la matrice dans un tableau Row-Major
– de taille n × n si [format] est full, ou
– de taille n × (k + 1) si [format] est band.
• n indique la taille de la matrice
• k indique le nombre de sous-diagonales non-nulles
• d est un tableau de taille n qui contient en sortie la diagonale principale de la matrice
tridiagonale (Hessenberg symétrique)
• e est un tableau de taille n qui contient en sortie la sous-diagonale de la matrice tridiagonale
dans ses n − 1 premiers éléments
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
            e[i] = A[2*k + (k+1) * i];
        }
        return 0;
    }

    // Apply Givens rotations to zero out the sub-diagonals
    for (int j = 0; j < n - 2; j++) {
        for (int i = MIN(j + k, n-1); i > j + 1; i--) {
            // Rotation coefficients
            double a = get(A,i-1,j); 
            double b = get(A,i,j);
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = b / r;

            // Apply rotation to rows i and i-1
            int start = MAX(0,i-k-1);
            int stop = i;

            double temp_i = 0;
            double temp_j = 0;
            // We need to store this value because it will be overwritten 
            // As the rotation is applied in-place and we only store the diagonal and sub-diagonal elements
            // (This value is used in the last iteration of the loop)
            double temp_diag = get_safe(A,i-1,stop);  // Todo : Alternativement : get(A, stop, i-1)
            
            // Notice that we make use of get_safe to avoid accessing elements outside the matrix :
            // If we try to access to an element outside the matrix, the macro returns 0
            // And if we try to access to an element in the upper band of the matrix, the macro returns the value of the element
            // at the symetric position in the lower band of the matrix

            // Similarly, set_safe is used to avoid writing outside the matrix
            // If we try to write outside the matrix or in the upper band, the macro does nothing

            for (int l = start; l < stop+1; l++) { 
                temp_i = l < stop ? c * get_safe(A,i-1,l) + s * get_safe(A,i,l) : c * temp_diag + s * get_safe(A,i,l);
                temp_j = l < stop ? - s * get_safe(A,i-1,l) + c * get_safe(A,i,l) : - s * temp_diag + c * get_safe(A,i,l);

                // We also need a way to remember the value of the last element of the row i-1 at the end of the loop 
                // that falls in the upper band of the matrix
                // Because we need it for the start of the column rotation
                // (as the matrix is in sym-band format, the value of the last element of the row i-1 is not stored in the matrix)
                // We need to store it in a temporary variable

                if(l < stop){
                    set_safe(A,i-1,l, temp_i);
                }
                else{
                    temp_diag = temp_i;
                }
                set_safe(A,i,l, temp_j);
            }
            

            // Apply rotation to columns i and i-1
            start = i-1;
            stop = MIN(i+k, n-1);

            // As previously said, we use the temp value stored at the first iteration of the loop
            // the other iterations does not suffer from this as their values are stored in the matrix
            for (int l = start; l < stop+1; l++) {
                temp_i =  l > start ? c * get_safe(A,l,i-1) + s * get_safe(A,l,i) : c * get_safe(A,l,i-1) + s * temp_diag;
                temp_j = l > start ? - s * get_safe(A,l,i-1) + c * get_safe(A,l,i) : -s * get_safe(A,l,i-1) + c * temp_diag;
                set_safe(A,l,i-1, temp_i);
                set_safe(A,l,i, temp_j);
            }

            // Be careful that when we make one zero appears in the matrix at the position (i,j)
            // The entry at position (i+(k-j) + 1,i) which is out of the band becomes non zero      
        }
    }
    
    // Extract the diagonal and sub-diagonal elements
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (k + 1) + k];
    }
    for (int i = 1; i < n; i++) {
        e[i] = A[(i + 1) * (k + 1) + k - 1];
    }

    return 0;
    
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
    double d[10], e[10];
    e[0] = 0;

    // Print the original matrix
    print_sym_band(A, 10, 6, "Original A");

    // Call the tridiagonalize_band function
    tridiagonalize_band(A, 10, 6, d, e);

    // Print the tridiagonalized matrix
    print_sym_band(A, 10, 6, "Tridiagonalized A");

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
    double d[4], e[4]; e[0] = 0;

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

int smain(int argc, char *argv[]){
    // Test the tridiagonalize function
    double A[12] ={0, 0, 8, 
                      0, 1, 4,  
                         4, 3, 1,
                            1, 2, 9};
    double d[4], e[4]; e[0] = 0;
    print_sym_band(A, 4, 2, "A");
    tridiagonalize_band(A, 4, 2, d, e);
    print_sym_band(A, 4, 2, "A");
    printf("d = [%f, %f, %f, %f]\n", d[0], d[1], d[2], d[3]);
    printf("e = [%f, %f, %f]\n", e[0], e[1], e[2]);
    return 0;
}


int main() {
    test_tridiagonalize_band_small();
    return 0;
}