#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include <lapacke.h>

#define TOL_PRINT 0e-15
#define MAX_SIZE_PRINT 25

void print_sym_band(double *A, int n, int b, char *name);
void print_band(double *A, int n, int k, char *name);
void print_full(double *A, int n, char *name);
void compute_eigenvalues(double *A, int n, int k);
#endif