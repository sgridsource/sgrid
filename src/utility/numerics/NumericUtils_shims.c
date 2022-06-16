/* NumericUtils_shims.c */
/* Wolfgang Tichy 6/2022 */

#include "sgrid.h"
#include "numerics.h"


/* Newton linesearch with vecfuncP that has an argument par for
   parameters (a pointer to void):
   returns its>=0   if ok, 
   returns -its-1<0 if error!!!    its = number of iterations done
   returns -its-itmax*10<0  if f=NAN.
   The funtion vecfuncP has to be of the form:
    void vecfuncP(int n, double *x, double *fvec, void *par)
                  here x[1...n], fvec[1...n], par can point to something
                  that contains parameters.
   Then newton_linesearch will find x s.t. vecfuncP = 0.               */
int newton_linesrch_itsP(double x[], int n, int *check,
                     void (*vecfuncP)(int, double [], double [], void *par),
                     void *par, int MAXITS, double TOLF)
{
  /* the old newton_linesrch_itsP had first enrty in x at x[1] */
  return newton_linesearch(n, x, check, vecfuncP,par,1, MAXITS, TOLF);
}



/* not done! */
int WT_newton(double *x, int n, int *check,
        void (*F_x)(int, double *x, double *Fx, void *par),
        void (*J_x_dx)(int, double *dx, double *Jdx, double *x, void *par),
        void *par, int MAXITS, double TOLF,
        int (*linSol)(int n, double *b, double *dx,
            void (*Jdx)(int, double *, double *, double *, void *),
            int (*precon)(int n, double *b, double *dx, double *x, void *par),
            double *x, void *par, int itmax, double tol),
        int (*precon)(int n, double *b, double *dx, double *x, void *par),
        int linitmax, double lintolfac)
{
  errorexit("not implemented in numerics");
  return -1;
}

/* bracket for root */
int zbrac_P(double (*func)(double,void *par), double *x1, double *x2,
            void *par)
{
  rt_brak(func,par, x1,x2, 50);
}

/* old way of finding root in bracket [x1,x2] */
int zbrent_itsP(double *x0, double (*func)(double,void *par),
                double x1, double x2, void *par, int ITMAX, double tol)
{

  return rtbrent_brak(x0, func, x1,x2, par, ITMAX, tol, 1);
}

/* lu_decomp but for a matrix that starts with a[1][1],
   indx also starts at indx[1].
   returns -i, if i=row with all zeros
   returns +j, if j=column with all zeros,
   0 otherwise */
int ludcmpSing(double **a, int n, int *indx, double *d)
{
  int parity, ret;
  double (*A)[n];

  A = (double (*)[n]) &(a[1][1]);

  ret = lu_decomp(n, A, indx+1, &parity);
  *d = parity;

  return ret;
}

/* LU solve for matrix, vector and indx, that all start at index 1 */
void lubksb(double **a, int n, int *indx, double b[])
{
  double (*A)[n];

  A = (double (*)[n]) &(a[1][1]);

  lu_solve(n, A, indx+1, b+1);
}

/* Brent's method */
double brent_with_pointer_to_pars(double ax, double bx, double cx, 
             double (*f)(double, void *ppointer), double tol,
	     double *xmin, void *parpointer)
{
  double fmin;

  minbrent_brak(f,parpointer, ax,bx,cx, 100,tol, xmin,&fmin, 1);
  return fmin;
}

/* bracket a min */
void mnbrak_with_pointer_to_pars(double *ax, double *bx, double *cx,
	    double *fa, double *fb, double *fc,
	    double (*func)(double, void *ppointer), void *parpointer)
{
  min_brak(func,parpointer, ax,bx,cx, fa,fb,fc, 1073741824);
}
