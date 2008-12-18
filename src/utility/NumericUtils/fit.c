/* fit.c */
/* Wolfgang Tichy, April 2005 */

/* test fitting routines from Numerical Recipes 

cc fit.c nrutil.c dsvdfit.c dsvdcmp.c dsvbksb.c dpythag.c -lm

*/

#include <math.h>
#include "nrutil.h"

void svdfit(double x[], double y[], double sig[], int ndata, 
	    double a[], int ma,
	    double **u, double **v, double w[], double *chisq,
	    void (*funcs)(double, double [], int));




/* recommended hack for multi-dimensional fits */
double xofi[100];
double yofi[100];




/* evaluate basis functions */
void poly2(double idouble, double p[], int np)
{
  int i = idouble;

  p[1] = 1;
  p[2] = xofi[i];
  p[3] = yofi[i];
  p[4] = p[2]*p[2];
  p[5] = p[3]*p[3];
  p[6] = p[3]*p[2];

  if (0) printf("poly2(%2d): %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
		i, p[1], p[2], p[3], p[4], p[5], p[6]);
}




/* test functions 2d */
double f2d(int i)
{
  static int n = 0;
  double f;
  double x = xofi[i];
  double y = yofi[i];

  f = pow(x*x + y*y + 1, 1);
  f = - exp(-x*x - y*y);

  if (0) printf("%3d:  f(%19.12e,%19.12e) = %19.12e\n", n++, x, y, f);
  return f;
}




/* main */
int main(int argc, char **argv)
{
  double min, fmin;
  int i, j, k;
  double h = 0.1;
  double chisq;
  double tol = 1e-0;

  int n = 16;                    /* number of point/value pairs */
  double *x = dvector(1, n);     /* points */
  double *f = dvector(1, n);     /* function values */
  double *sig = dvector(1, n);   /* sigmas */

  int ma = 6;                    /* number of coefficients */
  double *a = dvector(1, ma);    /* coefficients of fit */

  double **u = dmatrix(1, n, 1, ma);   // workspace
  double **v = dmatrix(1, ma, 1, ma);
  double *w = dvector(1, ma);

  /* initialize global 2d coordinates for centered 4x4 square */
  for (i = 1; i <= n; i++) {
    x[i] = i;
    j = (i-1) % 4;
    k = (i-1) / 4;
    xofi[i] = j*h - 1.5*h;
    yofi[i] = k*h - 1.5*h;

    f[i] = f2d(i);
    sig[i] = tol;
  }

  /* check */
  for (i = 1; i <= n; i++) {
    printf("i=%2d  x=%4.1f y=%4.1f  f(x,y) = %9.6f\n", 
	   i, xofi[i], yofi[i], f[i]);
  }
  
  /* fit */
  dsvdfit(x, f, sig, n, a, ma, u, v, w, &chisq, poly2);

  /* result */
  printf("chisq = %e\n", chisq);
  for (i = 1; i <= ma; i++) {
    printf("a%d = %10.3e\n", i-1, a[i]);
  }
    

  return 0;
}
