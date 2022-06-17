/* GSL_odeint.c */
/* Wolfgang Tichy 6/2022 */

#include "sgrid.h"
#include "numerics.h"


#ifdef GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#define G_SUCCESS GSL_SUCCESS

/* integrate from *xs to xe where we pass in the the n-dim vector y(xs)
   into y, here the first entries are y[0] and dt[0] */
int GSL_odeint(double *xs, double xe, int n, double *y,
        int (*derivsP)(double x, const double *y, double *dy, void *p),
        void *pars,
        double hstart, double epsabs, double epsrel,
        const gsl_odeiv2_step_type *step_type)
{
  int stat;

  gsl_odeiv2_system sys = { derivsP, NULL, n, pars };
  gsl_odeiv2_driver *drv;

  /* init GSL driver */
  drv = gsl_odeiv2_driver_alloc_y_new(&sys, step_type,
                                      hstart, epsabs, epsrel);

  /* integrate up to xe */
  stat = gsl_odeiv2_driver_apply(drv, xs, xe, y);

  /* free GSL driver memory */
  gsl_odeiv2_driver_free(drv);

  return stat;
}

/* integrate from *xs to xe where we give the the n-dim vector
   y(xs) = y with the Runge-Kutta Prince-Dormand (8, 9) method
   gsl_odeiv2_step_rk8pd */
int GSL_odeint_rk8pd(double *xs, double xe, int n, double *ystart,
        int (*derivsP)(double x, const double *y, double *dy, void *p),
        void *pars,
        double hstart, double epsabs, double epsrel)
{
  return GSL_odeint(xs, xe, n, ystart, derivsP,pars,
                    hstart,epsabs,epsrel, gsl_odeiv2_step_rk8pd);
}
#else
#define G_SUCCESS 0
#define gsl_odeiv2_step_type void
int GSL_odeint(double *xs, double xe, int n, double *y,
        int (*derivsP)(double x, const double *y, double *dy, void *p),
        void *pars,
        double hstart, double epsabs, double epsrel,
        gsl_odeiv2_step_type *step_type)
{
  errorexit("in order to compile with the GSL use MyConfig with\n"
            "DFLAGS += -DGSL\n"
            "SPECIALLIBS += -lgsl -lgslcblas\n");
  return 0;
}
int GSL_odeint_rk8pd(double *xs, double xe, int n, double *ystart,
        int (*derivsP)(double x, const double *y, double *dy, void *p),
        void *pars,
        double hstart, double epsabs, double epsrel)
{
  errorexit("in order to compile with the GSL use MyConfig with\n"
            "DFLAGS += -DGSL\n"
            "SPECIALLIBS += -lgsl -lgslcblas\n");
  return 0;
}
#endif


/*************************************************************************/
/* this chould also go into utility/numerics/NumericUtils_shims.c */
/*************************************************************************/

/* struct and conversion of derivs to derivsP for odeintegrate */
struct Shim_derivs {
  void (*derivs)(double x, double y[], double dy[]);
};
int Shim_derivs_to_derivsP(double x, const double *y, double *dy, void *p)
{
  double *yy = (double *) y;
  struct Shim_derivs *par = p;
  par->derivs(x, yy-1, dy-1);
  return G_SUCCESS;
}

/* here the entries are y[1], xp[1], yp[1][1] */
double odeintegrate(double y[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status)
{
  double xs, xf, xe;
  int stat, j, k;
  struct Shim_derivs par[1];

  par->derivs = derivs;

  *status = 0;
  *kcount = 0;;

  if(kmax>0)
  {
    /* save initial x and y in xp, yp */
    xs = x1;
    xp[1] = xs;
    for(j=1; j<=nvar; j++) yp[j][1] = y[j];
    *kcount = *kcount + 1;
  }

  //printf("x1=%g x2=%g\n", x1,x2);

  /* do the integral over the full region */
  xs = x1;
  xe = x2;
  stat = GSL_odeint_rk8pd(&xs, xe, nvar, y+1, Shim_derivs_to_derivsP,par,
                          h1, eps, eps);
  xf = xs; /* save final x */

  //printf("xs=%g xf=%g xe=%g\n ", xs, xf, xe);

  if(kmax>1)
  {
    /* now integrate again to create the damn nr-data (what a waste) */
    xs = x1;
    for(j=1; j<=nvar; j++) y[j] = yp[j][1];

    for(k=2; k<=kmax; k++)
    {
      xe = x1 + (k-1) * (xf-x1) / (kmax-1);

      //printf("k=%d: xs=%g xe=%g\n ", k, xs,xe);

      GSL_odeint_rk8pd(&xs, xe, nvar, y+1, Shim_derivs_to_derivsP,par,
                       h1, eps, eps);
      xp[k] = xs;
      for(j=1; j<=nvar; j++) yp[j][k] = y[j];
      *kcount = *kcount + 1;

      //printf("     xs=%g xe=%g\n ", xs,xe);
      if(xs<xe) break;
    }
  }
  //printf("kmax=%d *kcount=%d\n", kmax, *kcount);
  return xs;
}


/* convert from starting at y[0] to starting at y[1] */
struct Shim_y_to_ym1 {
  int (*derivsP)(double x, const double y[], double dy[], void *p);
  void *orig_par;
};
int Shim_derivsP_y_to_ym1(double x, const double *y, double *dy, void *p)
{
  struct Shim_y_to_ym1 *shim = p;
  return shim->derivsP(x, y-1, dy-1, shim->orig_par);
}

/* here the entries are y[1], xp[1], yp[1][1] */
double odeintegrateP(double y[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	int (*derivsP)(double x, const double *y, double *dy, void *p),
	void *par,
	void (*rkqsP)(double [], double [], int, double *, double, double, double [],
	double *, double *,
	int (*derP)(double x, const double *y, double *dy, void *p), void *p),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status)
{
  double xs, xf, xe;
  int stat, j, k;
  struct Shim_y_to_ym1 shim[1];

  shim->derivsP  = derivsP;
  shim->orig_par = par;

  *status = 0;
  *kcount = 0;;

  if(kmax>0)
  {
    /* save initial x and y in xp, yp */
    xs = x1;
    xp[1] = xs;
    for(j=1; j<=nvar; j++) yp[j][1] = y[j];
    *kcount = *kcount + 1;
  }

  //printf("x1=%.17g x2=%.17g\n", x1,x2);
  //printf("y[2]=%g\n", y[2]);

  /* do the integral over the full region */
  xs = x1;
  xe = x2;
  stat = GSL_odeint_rk8pd(&xs, xe, nvar, y+1, Shim_derivsP_y_to_ym1,shim,
                          h1, eps, eps);
  xf = xs; /* save final x */

  //printf("xs=%.17g xf=%.17g xe=%.17g\n ", xs, xf, xe);
  //printf("   ==> y[2]=%g\n", y[2]);

  if(kmax>1)
  {
    /* now integrate again to create the damn nr-data (what a waste) */
    xs = x1;
    for(j=1; j<=nvar; j++) y[j] = yp[j][1];

    for(k=2; k<=kmax; k++)
    {
      xe = x1 + (k-1) * (xf-x1) / (kmax-1);

      //printf("k=%d: xs=%.17g xe=%.17g\n ", k, xs,xe);
      //printf("      y[2]=%g\n", y[2]);

      GSL_odeint_rk8pd(&xs, xe, nvar, y+1, Shim_derivsP_y_to_ym1,shim,
                       h1, eps, eps);
      xp[k] = xs;
      for(j=1; j<=nvar; j++) yp[j][k] = y[j];
      *kcount = *kcount + 1;

      //printf("     xs=%.17g xe=%.17g\n ", xs,xe);
      //printf("     y[2]=%g\n", y[2]);
      if(xs<xe) break;
    }
  }
  //printf("kmax=%d *kcount=%d\n", kmax, *kcount);

  return xs;
}
