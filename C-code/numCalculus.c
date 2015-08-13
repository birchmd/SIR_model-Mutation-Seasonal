/* This file contains a routine to compute
   derivatives and integrals of a function which has been
   descritized onto a mesh.
 */
#include <math.h>
#include "numCalculus.h"

/* R interface for trapz function. */
void trapz_R(int *N, double *x, double *y, double *result){
	result[0] = trapz(N[0], x, y);
}

/* The following routine computes the integral of a function
   using the trapezoidal rule.
   -N is the length of the arrays x and y
   -x is the array of points at which the function is evaluated
   -y is the array of function values, y=f(x)
*/
double trapz(int N, double *x, double *y){
	double sum;
	int j;
	
	sum = (double)0.0;
	for (j = 1; j < N; j++){
		sum += (x[j] - x[j-1])*(y[j] + y[j-1]);
	}
	sum = sum * (double)0.5;
	return sum;
}

/* R interface for d2dx2. Returns array of derivatives
*/
void d2dx2_R(double *f, int *N, double* x, int* BC, double* result){
	int j;
	
	for (j = 0; j < N[0]; j++){
		result[j] = d2dx2(f, j, N[0], x, BC);
	}
}

/* This function compute the second derivative of the function f
   (which is descritized into an array of length N) 
   at the position x_j using the centred difference approximation
   on a non-uniform mesh. BC gives the boundary conditions:
   BC = 0 - Dirichlet
        1 - Neumann
		2 - Linearly approximate endpoints
*/
double d2dx2(double *f, int j, int N, double* x, int *BC){
	double d21, d20, d10;

	if (j == N){ return sqrt(-1); } //NaN, array index out of range
	if (j==0){
		if( BC[0] == 0 ){
			return BC_dir(f, j, N, x, BC);
		}else if( BC[0] == 1 ){
			return BC_neu(f, j, N, x, BC);
		}else if( BC[0] == 2 ){
			return BC_lin(f, j, N, x, BC);
		} else {
			return 0;
		}
	} else if (j == N-1){
		if( BC[1] == 0 ){
			return BC_dir(f, j, N, x, BC);
		}else if( BC[1] == 1 ){
			return BC_neu(f, j, N, x, BC);
		}else if( BC[1] == 2 ){
			return BC_lin(f, j, N, x, BC);
		} else {
			return 0;
		}
	}
	d21 = x[j+1] - x[j];
	d20 = x[j+1] - x[j-1];
	d10 = x[j] - x[j-1];
	return ((double)2.0/(d21*d20))*f[j+1] - ((double)2.0/(d21*d10))*f[j] + ((double)2.0/(d10*d20))*f[j-1];
}
double BC_dir(double *f, int j, int N, double* x, int *BC){
	double d21, d20, d10;

	if (j == N-1){ // points off array are zero
		d10 = x[j] - x[j-1];
		d21 = d10;
		d20 = (double)2.0*d10;
		return -((double)2.0/(d21*d10))*f[j] + ((double)2.0/(d10*d20))*f[j-1];
	} else if (j == 0){
		d21 = x[j+1] - x[j];
		d20 = (double)2.0*d21;
		d10 = d21;
		return ((double)2.0/(d21*d20))*f[j+1] - ((double)2.0/(d21*d10))*f[j];
	} else {
		return 0;
	}
}
double BC_neu(double *f, int j, int N, double* x, int *BC){
	double d21, d20, d10;
	
	if (j == N-1){ // df/dx = 0 at end points
		d10 = x[j] - x[j-1];
		d21 = d10;
		d20 = (double)2.0*d10;
		return ((double)2.0/(d21*d20))*f[j] - ((double)2.0/(d21*d10))*f[j] + ((double)2.0/(d10*d20))*f[j-1];
	} else if (j == 0){
		d21 = x[j+1] - x[j];
		d20 = (double)2.0*d21;
		d10 = d21;
		return ((double)2.0/(d21*d20))*f[j+1] - ((double)2.0/(d21*d10))*f[j] + ((double)2.0/(d10*d20))*f[j];
	} else {
		return 0;
	}
}
double BC_lin(double *f, int j, int N, double* x, int *BC){
	double d1, d2;
	if (j == N-1){ // approximate the end points linearly
		d1 = d2dx2(f, N-3, N, x, BC);
		d2 = d2dx2(f, N-2, N, x, BC);
		return ((x[N-1] - x[N-2])/(x[N-3] - x[N-2])) * (d1 - d2) + d2;
	} else if (j == 0){
		d1 = d2dx2(f, 1, N, x, BC);
		d2 = d2dx2(f, 2, N, x, BC);
		return (x[0] - x[1])/(x[2] - x[1]) * (d2 - d1) + d1;
	} else {
		return 0;
	}
}