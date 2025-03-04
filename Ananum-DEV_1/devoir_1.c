#include "utils.h"


#define get(A, i, j) A[(j-i) + i * (k + 1) + k]
#define set(A, i, j, val) A[(j-i)+ i * (k + 1) + k] = val
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
        for(int i = 0; i < n - 1; i++){
            e[i] = A[2*k + (k+1) * i];
        }
        return 0;
    }
    
    // Apply Givens rotations to zero out the sub-diagonals
    for (int j = 0; j < n - 1; j++) {
        for (int i = j + k; i > j + 1; i--) {
            printf("i = %d, j = %d\n", i, j);
            print_sym_band(A, 4, 2, "A");

            double a = get(A,i-1,j);
            double b = get(A,i,j);
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = b / r;

            printf("a = %f, b = %f, r = %f, c = %f, s = %f\n", a, b, r, c, s);

            // Apply rotation to rows i and i-1
            double start = MAX(0,i-1-k);
            double stop = i-1;
            for (int l = start; l < stop; l++) {
                double temp_i = l == start ? c * get(A,i-1,l) : c * get(A,i-1,l) - s * get(A,i,l);
                set(A,i-1,l, temp_i);
            }

            start += 1;
            stop += 1;
            for (int l = start; l < stop; l++) {
                double temp_j = l == start ? s * get(A,i-1,l) : s * get(A,i-1,l) + c * get(A,i,l);
                set(A,i,l, temp_j);
            }

            // Apply rotation to columns i and i-1
            
    }

    // Extract the diagonal and sub-diagonal elements
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (k + 1) + k];
    }
    for (int i = 0; i < n - 1; i++) {
        e[i] = A[(i + 1) * (k + 1) + k - 1];
    }

    return 0;
}

int main(int argc, char *argv[]){
    // Test the tridiagonalize function
    double A[12] ={0, 0, 2, 
                      0, 1, 3,  
                         1, 3, 1,
                            1, 2, 1};
    double d[4], e[3];
    print_sym_band(A, 4, 2, "A");
    tridiagonalize_band(A, 4, 2, d, e);
    print_sym_band(A, 4, 2, "A");
//    printf("d = [%f, %f, %f, %f]\n", d[0], d[1], d[2], d[3]);
//    printf("e = [%f, %f, %f]\n", e[0], e[1], e[2]);
    return 0;
}