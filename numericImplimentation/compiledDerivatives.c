/* This file contains the routines to compute
   the derivatives and Jacobian for the SI-model
   with mutation (modelled by diffusion)
 */
   
#include <stdio.h>
#include "../C-code/numCalculus.c"

#define PI 3.14159265358979323846264338327950288419716939937510582

static double parms[9];
#define nu parms[0]
#define mu parms[1]
#define betaMult parms[2]
#define betaExp parms[3]
#define alphaMin parms[4]
#define alphaMax parms[5]
#define diffConst parms[6]
#define seasAmp parms[7]
#define seasPer parms[8]

/* initializer */
void compiledDerivatives(void (* odeparms)(int *, double *))
{
	int N=9;
	odeparms(&N, parms);
}

/* This function fills the array arr with N values, from arr[0]=minVal
   to arr[N-1]=maxVal.
*/
void seq(double minVal, double maxVal, int N, double *arr){
	double dVal;
	int i;
	
	dVal = (maxVal - minVal)/((double)(N-1));
	arr[0] = minVal;
	for (i = 1; i < N; i++){
		arr[i] = arr[i-1] + dVal;
	}
	return;
}

/* beta(alpha) */
double beta(double alpha){
	return betaMult * pow(alpha, (double)1.0/betaExp);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
	double S;
	double *i;
	double alpha[*neq-1];
	double beta_times_i[*neq-1];
	double g; //seasonal forcing
	int j;
	int BC[2];
	
	BC[0] = 1; //Neumann boundary condition at 0
	BC[1] = 0; //Dirichlet boundary condition at other end

	//floor at zero to prevent negative numbers
	for (j = 0; j < *neq; j++){ 
		if(y[j] < 0.0){ y[j] = (double)0.0; }
	}
	
	S = y[0]; //get current number of susceptibles
	i = &y[1]; //save memory by only pointing to the infectious density values
	g = (double)1.0 + seasAmp*sin((double)(2.0*PI)*(t[0]/seasPer));
	seq(alphaMin,alphaMax,*neq-1,alpha); //create array of alpha values
	for (j = 0; j < *neq-1; j++){ 
		beta_times_i[j] = beta(alpha[j])*i[j]; //calculate beta(alpha)*i(alpha,t)
		ydot[j+1] = S*beta_times_i[j]*g - alpha[j]*i[j] - mu*i[j] + diffConst*d2dx2(i,j,*neq-1,alpha,BC); // di/dt
	}
	ydot[0] = nu - S*trapz(*neq-1,alpha,beta_times_i)*g - mu*S; // dS/dt
}