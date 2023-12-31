#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXSTP 10000
#define TINY 1.0e-250

/* odeintegrate.c is like odeint.c, except that no global vars are needed,
 and that it returns 
  x2 if all seems ok and status=0
  x  if "Step size too small in odeintegrate" and status=-1
  x  if "Too many steps in routine odeintegrate" and status=-2. */

/* From Fortran version:
c Integrator driver with adaptive stepsize control.  Integrate the NVAR 
c starting values YSTART from X1 to X2 with accuracy EPS, storing 
c intermediate results in the common block /PATH/.  H1 should be set as a 
c guessed first step size, HMIN as the minimum allowed step size (zero okay).
c On output, NOK and NBAD are the number of good and bad (but retried and 
c fixed) steps taken, and YSTART is replaced by values at the end of the
c integration interval.  DERIVS is the user-supplied subroutine for 
c calculating the right hand side derivative, while RKQS is the name of
c the stepper routine to be used.  PATH contains its own information about
c how often an intermediate value is to be stored.
*/
/* Variablen fuer odeintegrate.c
 double x1,x2;            intial and final x-point
 double *y;               The functions y1, y2, ...
 double *dy;              The functions' derivs dy1, dy2, ...
 int nvar;                The number of functions
 double eps, h1, hmin;    error, first step, minimum step
 int nok,nbad;            # of ok steps, # of bad steps
 int kmax;                max # of points outputed by odeint
 int kcount;              # of points outputed by odeint
 double *xp,**yp;         points outputed by odeint
 double dxsav;            approx. x-distance between points
Allocate y, dy, xp, yp like this:
 y =vector(1,nvar);
 dy=vector(1,nvar);
 xp=vector(1,kmax);
 yp=matrix(1,nvar,1,kmax);
*/
double odeintegrate(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status)
{
	int kount=*kcount;
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			*kcount = kount;
			*status = 0;
			return x2;
		}
		if(fabs(hnext) <= hmin) /* nrerror("Step size too small in odeintegrate"); */
		{
                  free_vector(dydx,1,nvar);
                  free_vector(y,1,nvar);
                  free_vector(yscal,1,nvar);
		  *status = -1;
		  return x;
		}
		h=hnext;
	}
	/* nrerror("Too many steps in routine odeintegrate"); */
	free_vector(dydx,1,nvar);
	free_vector(y,1,nvar);
	free_vector(yscal,1,nvar);
	*status = -2;
	return x;
}

/* version with extra parameter arg for derivs */
double odeintegrateP(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	int (*derivsP)(double x, const double *y, double *dy, void *p),
	void *par,
	void (*rkqsP)(double [], double [], int, double *, double, double, double [],
	double *, double *,
	int (*derP)(double x, const double *y, double *dy, void *p), void *p),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status)
{
	int kount=*kcount;
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivsP)(x,y,dydx, par);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqsP)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivsP, par);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			*kcount = kount;
			*status = 0;
			return x2;
		}
		if(fabs(hnext) <= hmin) /* nrerror("Step size too small in odeintegrate"); */
		{
                  free_vector(dydx,1,nvar);
                  free_vector(y,1,nvar);
                  free_vector(yscal,1,nvar);
		  *status = -1;
		  return x;
		}
		h=hnext;
	}
	/* nrerror("Too many steps in routine odeintegrate"); */
	free_vector(dydx,1,nvar);
	free_vector(y,1,nvar);
	free_vector(yscal,1,nvar);
	*status = -2;
	return x;
}

/* odeintegrate version from Hannes Ruter that uses a different yscal */
double odeintegrate_HR(double ystart[], int nvar, double x1, double x2,
    double eps, double h1, double hmin, int *nok, int *nbad,
    void (*derivs)(const double, const double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
    double *, double *, void (*)(const double, const double [], double [])),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status)
{
    int kount;
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;

	double *max_absy=vector(1,nvar);
    double max_abs_scale = 1e-5;

	for (i=1;i<=nvar;i++){
	    y[i]=ystart[i];
	    max_absy[i] = fabs(y[i]);
	}
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++) {
		  if(fabs(y[i])>max_absy[i]){
		    max_absy[i] = fabs(y[i]);
		  }
		}
		for (i=1;i<=nvar;i++)
			yscal[i]=max_abs_scale*max_absy[i]+fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) {
			  yp[i][kount]=y[i];
			}
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
        (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			*kcount = kount;
			*status = 0;
			return x2;
		}
		if(fabs(hnext) <= hmin) /* nrerror("Step size too small in odeintegrate"); */
		{
                  free_vector(dydx,1,nvar);
                  free_vector(y,1,nvar);
                  free_vector(yscal,1,nvar);
                  *kcount = kount;
		  *status = -1;
		  return x;
		}
		h=hnext;
	}
	/* nrerror("Too many steps in routine odeintegrate"); */
	free_vector(dydx,1,nvar);
	free_vector(y,1,nvar);
	free_vector(yscal,1,nvar);
	*kcount = kount;
	*status = -2;
	return x;
}

#undef MAXSTP
#undef TINY
#undef NRANSI
