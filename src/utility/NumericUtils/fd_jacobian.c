/* (c) Wolfgang Tichy 6.6.2003 */
/* compute Jacobin with finite differences */

#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-4    /* should be approx square root of machine precision */
#define HMIN 1.0e-10  /* min h we use in fin. diff. computation of derivs */

void fd_jacobian(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []))
{
	int i,j;
	double h,temp,*f;

	f=vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h < HMIN) h=HMIN; /* make sure h is never below HMIN */
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_vector(f,1,n);
}

/* same as above, but pass on void *par */
void fd_jacobianP(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double [], void *par), void *par)
{
	int i,j;
	double h,temp,*f;

	f=vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h < HMIN) h=HMIN; /* make sure h is never below HMIN */
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f, par);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_vector(f,1,n);
}


#undef HMIN
#undef EPS
#undef NRANSI
