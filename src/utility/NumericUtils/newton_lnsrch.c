/* (c) Wolfgang Tichy 6.6.2003 */
/* adapted from numrec newt.c + tolerances changed */

#include <math.h>
#define NRANSI
#include "nrutil.h"

/* MAXITS and TOLF are passed in */
/*
#define MAXITS 200
#define TOLF 1.0e-4
*/

#define TOLMIN 1.0e-13
#define TOLX 1.0e-14
#define STPMX 100.0

/* global vars */
int newton_lnsrch_nn;
double *newton_lnsrch_fvec;
void (*newton_lnsrch_nrfuncv)(int n, double v[], double f[]);

#define FREERETURN {free_vector(newton_lnsrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newton_lnsrch(double x[], int n, int *check,
			void (*vecfunc)(int, double [], double []), 
			int MAXITS, double TOLF)
{
	void fd_jacobian(int n, double x[], double newton_lnsrch_fvec[], double **df,
		void (*vecfunc)(int, double [], double []));
	double f_to_minimize(double x[]);
	void lnsrch_double(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,its,j,*indx;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	newton_lnsrch_fvec=vector(1,n);
	newton_lnsrch_nn=n;
	newton_lnsrch_nrfuncv=vecfunc;
	f=f_to_minimize(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(newton_lnsrch_fvec[i]) > test) test=fabs(newton_lnsrch_fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fd_jacobian(n,x,newton_lnsrch_fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*newton_lnsrch_fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -newton_lnsrch_fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		lnsrch_double(n,xold,fold,g,p,x,&f,stpmax,check,f_to_minimize);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(newton_lnsrch_fvec[i]) > test) test=fabs(newton_lnsrch_fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
	}
	nrerror("MAXITS exceeded in newt");
}

/* same as above but return error code:
   returns its>=0 if ok, 
   returns -its<0 if error!!!    its = number of iterations done */
#undef FREERETURN
#define FREERETURN {free_vector(newton_lnsrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return its;}
#define FREERETURNERROR {free_vector(newton_lnsrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return -its;}
int newton_lnsrch_its(double x[], int n, int *check,
			void (*vecfunc)(int, double [], double []), 
			int MAXITS, double TOLF)
{
	void fd_jacobian(int n, double x[], double newton_lnsrch_fvec[], double **df,
		void (*vecfunc)(int, double [], double []));
	double f_to_minimize(double x[]);
	void lnsrch_double(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,*indx;
	int its=0;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	newton_lnsrch_fvec=vector(1,n);
	newton_lnsrch_nn=n;
	newton_lnsrch_nrfuncv=vecfunc;
	f=f_to_minimize(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(newton_lnsrch_fvec[i]) > test) test=fabs(newton_lnsrch_fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fd_jacobian(n,x,newton_lnsrch_fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*newton_lnsrch_fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -newton_lnsrch_fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		lnsrch_double(n,xold,fold,g,p,x,&f,stpmax,check,f_to_minimize);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(newton_lnsrch_fvec[i]) > test) test=fabs(newton_lnsrch_fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
	}
	/* nrerror("MAXITS exceeded in newt"); */
	FREERETURNERROR
}
#undef FREERETURNERROR
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
