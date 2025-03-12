#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include <lapacke.h>
#include "devoir_1.h"

#define TOL_PRINT 0e-15
#define MAX_SIZE_PRINT 25

void print_sym_band(double *A, int n, int b, char *name);
void print_band(double *A, int n, int k, char *name);
void print_full(double *A, int n, char *name);
void compute_eigenvalues(double *A, int n, int k);
void print_tridiag(double *d, double *e, int n, char *name);
void test_tridiagonalize_band_small();
void test_tridiagonalize_band_large();
void test_tridiagonalize_band();
void test_qr_band();
void test_qr_band_small();
void test_step_qr_tridiag();

#endif