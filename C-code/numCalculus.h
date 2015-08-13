#ifndef NUMCALC_H //Prevent functions from being defined twice
#define NUMCALC_H

void trapz_R(int *N, double *x, double *y, double *result);
double trapz(int N, double *x, double *y);
void d2dx2_R(double *f, int *N, double* x, int *BC, double* result);
double d2dx2(double *f, int j, int N, double* x, int *BC);
double BC_dir(double *f, int j, int N, double* x, int *BC);
double BC_neu(double *f, int j, int N, double* x, int *BC);
double BC_lin(double *f, int j, int N, double* x, int *BC);

#endif