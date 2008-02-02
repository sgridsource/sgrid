#include <math.h>
#define NRANSI
#include "nrutil.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4


/* From Fortran version:
	SUBROUTINE RKQS(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c Fifth-order Runge-Kutta step with monitoring of local truncation error
c to ensure accuracy and adjust step size
c
c     INPUTS:
c         Y(N)     dependent variable vector at start of step
c         DYDX(N)  derivative vector at start of step
c         X        independent variable at start of step
c         HTRY     stepsize to be attempted
c         EPS      accuracy
c         YSCAL(N) vector against which the error is scaled
c     
c     OUTPUTS
c         Y(N)     dependent variable vector at end of step
c         X        independent variable at end of step
c         HDID     step size used
c         HNEXT    first guess at next step size
*/
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	void rkck(double y[], double dydx[], int n, double x, double h,
		double yout[], double yerr[], void (*derivs)(double, double [], double []));
	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=vector(1,n);
	ytemp=vector(1,n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) nrerror("stepsize underflow in rkqs");
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	free_vector(ytemp,1,n);
	free_vector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
