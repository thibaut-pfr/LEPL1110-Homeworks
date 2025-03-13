#ifndef CODE4_H
#define CODE4_H

double *create_matrix(int nx, int ny, double lx, double ly, int storage);
int inverse_power(double *A, int n, int k, double *mu, double *v);
int deflation(double *L, int n, int k, int nb, double *eigws, double *eigvs);
int get_storage_format();

#endif
