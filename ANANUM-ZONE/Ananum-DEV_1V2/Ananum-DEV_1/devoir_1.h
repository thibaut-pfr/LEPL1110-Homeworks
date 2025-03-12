#ifndef DEVOIR_1_H
#define DEVOIR_1_H

int tridiagonalize_band(double *A, int n, int k, double *d, double *e);
int step_qr_tridiag(double *d, double *e, int m, double eps);
int qr_eigs_band(double *A, int n, int k, double eps, int max_iter);

#endif