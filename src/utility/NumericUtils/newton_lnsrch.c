/* (c) Wolfgang Tichy 6.6.2003 */
/* adapted from numrec newt.c + tolerances changed */

#include <stdio.h>
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

/* same as above with error code, but do not use global vars:
   returns its>=0   if ok, 
   returns -its-1<0 if error!!!    its = number of iterations done
   returns -its-MAXITS*10<0  if f=NAN                               */
#undef FREERETURN
#define FREERETURN {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return its;}
#define FREERETURNERROR {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return -its-1;}
int newton_linesrch_its(double x[], int n, int *check,
			void (*vecfunc)(int, double [], double []), 
			int MAXITS, double TOLF)
{
	void fd_jacobian(int n, double x[], double newton_linesrch_fvec[], double **df,
		void (*vecfunc)(int, double [], double []));
	double f_to_minimize_in_linesrch(double x[], int n, 
	                         double *newton_linesrch_fvec,
	                         void (*vecfunc)(int, double [], double []));
	void linesrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, 
		 double (*func)(double [], int, double *newton_linesrch_fvec,
		                void (*vecfunc)(int, double [], double [])),
		 double *newton_linesrch_fvec, 
		 void (*vecfunc)(int, double [], double []));
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,*indx;
	int its=0;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
        double *newton_linesrch_fvec;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	newton_linesrch_fvec=vector(1,n);
	f=f_to_minimize_in_linesrch(x, n, newton_linesrch_fvec, vecfunc);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fd_jacobian(n,x,newton_linesrch_fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*newton_linesrch_fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -newton_linesrch_fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		linesrch(n,xold,fold,g,p,x,&f,stpmax,check,
		         f_to_minimize_in_linesrch, 
		         newton_linesrch_fvec, vecfunc);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
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

/* same as newton_linesrch_its but now we use a vecfuncP with one more 
   argument par for parameters (a pointer to void):
   returns its>=0   if ok, 
   returns -its-1<0 if error!!!    its = number of iterations done
   returns -its-MAXITS*10<0  if f=NAN.
   The funtion vecfuncP has to be of the form:
    void vecfuncP(int n, double *x, double *fvec, void *par)
                  here x[1...n], fvec[1...n], par can point to something
                  that contains parameters.
   Then newton_linesrch_itsP will find x s.t. vecfuncP = 0.               */
#undef FREERETURN
#define FREERETURN {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return its;}
#define FREERETURNERROR {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return -its-1;}
int newton_linesrch_itsP(double x[], int n, int *check,
			 void (*vecfuncP)(int, double [], double [], void *par),
			 void *par, int MAXITS, double TOLF)
{
	void fd_jacobianP(int n, double x[], double newton_linesrch_fvec[], double **df,
		void (*vecfuncP)(int, double [], double [], void *par),
		void *par);
	double f_to_minimize_in_linesrchP(double x[], int n, 
	                         double *newton_linesrch_fvec,
	                         void (*vecfuncP)(int, double [], double [], void *par),
	                         void *par);
	void linesrchP(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, 
		 double (*func)(double [], int, double *newton_linesrch_fvec,
		                void (*vecfuncP)(int, double [], double [], void *par), 
                                void *par),
		 double *newton_linesrch_fvec, 
		 void (*vecfuncP)(int, double [], double [], void *par),
		 void *par);
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,*indx;
	int its=0;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
        double *newton_linesrch_fvec;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	newton_linesrch_fvec=vector(1,n);
	f=f_to_minimize_in_linesrchP(x, n, newton_linesrch_fvec, vecfuncP, par);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fd_jacobianP(n,x,newton_linesrch_fvec,fjac,vecfuncP, par);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*newton_linesrch_fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -newton_linesrch_fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		linesrchP(n,xold,fold,g,p,x,&f,stpmax,check,
		         f_to_minimize_in_linesrchP, 
		         newton_linesrch_fvec, vecfuncP, par);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
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


/***********************************************************/
/* functions to determine if we are within a certain range */
/***********************************************************/

/* check if vec is in a certain range: lo[i] <= vec[i] <= hi[i] */
int newton_lnsrch_vec_within_range(int n, double vec[],
                                   double hi[], double lo[])
{
  int i;
  int ret=1;
  
  for(i=1; i<=n; i++)
  {
    if(vec[i]>hi[i]) ret=0;
    if(vec[i]<lo[i]) ret=0;
  }
  return ret;
}

/* use vector vec and range lo[i] <= vec[i] <= hi[i] to obtain vector 
   vb on boundary and vi slightly inside this range */
void newton_lnsrch_set_vecs_for_lininterp(int n, double vec[], 
          double hi[], double lo[], double *vb, double *vi) 
{
  int i;
  double *cent;
  double *dir;
  double s, sc;

  /* alloc */
  cent=vector(1,n);
  dir =vector(1,n);

  for(i=1; i<=n; i++) 
  {
    /* pos. vec. of center */
    cent[i] = 0.5*(hi[i]+lo[i]);

    /* set directional vec dir */
    dir[i] = vec[i] - cent[i];
    
    /* find min s such that vb[i] = cent[i] + dir[i]*s is on boundary */
    s=1;
    if(dir[i]!=0)
    {
      if(vec[i]>hi[i])  sc = (hi[i] - cent[i])/dir[i];
      if(vec[i]<lo[i])  sc = (lo[i] - cent[i])/dir[i];
      if(sc<s) s=sc;
    }
  }

  /* set vb on boundary and vi a little bit inside boundary */
  for(i=1; i<=n; i++)
  {
    vb[i] = cent[i] + dir[i]*s;
    vi[i] = cent[i] + dir[i]*s*0.999999;
    if(vb[i]>hi[i]) vb[i]=hi[i]; /* make sure we stay in range */
    if(vb[i]<lo[i]) vb[i]=lo[i];
  }
  printf("newton_lnsrch_set_vecs_for_lininterp:\n");
  printf("vec=( ");
  for(i=1; i<=n; i++) printf("%g ", vec[i]);
  printf(")\n");
  printf(" vb=( ");
  for(i=1; i<=n; i++) printf("%g ", vb[i]);
  printf(")\n");
  printf(" vi=( ");
  for(i=1; i<=n; i++) printf("%g ", vi[i]);
  printf(")\n");

  /* free */
  free_vector(cent,1,n);
  free_vector(dir,1,n);
}

/* use function values fvb and fvi at vb and vi to find function value
   fvec at vec by linear interolation */
void newton_lnsrch_get_fvec_by_lininterp(int n, double vec[], 
          double vb[], double vi[], double *fvec, double fvb[], double fvi[]) 
{
  int i;
  double h, s;

  /* get length of vb[i] - vi[i] and vec[i] - vb[i]*/
  h=s=0;
  for(i=1; i<=n; i++)
  {
    double dum;
    
    dum = vb[i] - vi[i];
    h += dum*dum;
    
    dum = vec[i] - vb[i];
    s += dum*dum; 
  }
  h=sqrt(h);
  s=sqrt(s);

  /* use lengths to find fvec from fvb and fvi by lin. interp. */
  for(i=1; i<=n; i++) fvec[i] = fvb[i] + ((fvb[i]-fvi[i])/h)*s;

  printf("newton_lnsrch_get_fvec_by_lininterp:\n");
  printf("fvec=( ");
  for(i=1; i<=n; i++) printf("%g ", fvec[i]);
  printf(")\n");
  printf(" fvb=( ");
  for(i=1; i<=n; i++) printf("%g ", fvb[i]);
  printf(")\n");
  printf(" fvi=( ");
  for(i=1; i<=n; i++) printf("%g ", fvi[i]);
  printf(")\n");
}
