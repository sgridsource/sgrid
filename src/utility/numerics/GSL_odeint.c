/* GSL_odeint.c */
/* Wolfgang Tichy 6/2022 */

#include "sgrid.h"
#include "numerics.h"


#ifdef GSL
#define SUCCESS GSL_SUCCESS
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

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
#define SUCCESS 0
#define gsl_odeiv2_step_type int
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
  return SUCCESS;
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
  double xs;
  int stat, j, k;
  struct Shim_derivs par[1];

  par->derivs = derivs;

  *status = 0;
  *kcount = kmax;

printf("x1=%g x2=%g\n", x1,x2);
  xs = x1;
  for(k=1; k<=kmax; k++)
  {
    double xe = x1 + k * (x2-x1) / kmax;
printf("k=%d: xs=%g xe=%g\n ", k, xs,xe);

    stat = GSL_odeint_rk8pd(&xs, xe, nvar, y+1,
                            Shim_derivs_to_derivsP,par,
                            h1, eps, eps);
    xp[k] = xs;
    for(j=1; j<=nvar; j++) yp[j][k] = y[j];
printf("     xs=%g xe=%g\n ", xs,xe);
    if(xs<xe) break;
  }
  return xs;
}
