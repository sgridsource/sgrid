/* sgrid_NumericUtils.h */
/* Wolfgang Tichy */


#include "nrutil.h"


/* Functions */

/* for Newton-Raphson with line searches */
void newton_lnsrch(double x[], int n, int *check,
                   void (*vecfunc)(int, double [], double []),
               	   int MAXITS, double TOLF);
void fd_jacobian(int n, double x[], double fvec[], double **df,
  	         void (*vecfunc)(int, double [], double []));
int newton_lnsrch_its(double x[], int n, int *check,
		      void (*vecfunc)(int, double [], double []), 
		      int MAXITS, double TOLF);
  	  	               	   
/* 1D and 3D integrals */
double integral(double (*func)(double), double a, double b, double s, int max);
double rombintegral(double (*func)(double), double a, double b,
                    double eps, int max);
double integral3D(double (*int_meth)(double (*f_int)(double), 
                                     double a, double b,
                                     double eps, int max),
                  double (*func)(double, double, double), 
                  double x_limit1, double x_limit2,
                  double (*y_limit1)(double x), double (*y_limit2)(double x),
                  double (*z_limit1)(double x,double y), 
                  double (*z_limit2)(double x,double y),
                  double sx, double sy, double sz, 
                  int maxx, int maxy, int maxz);

/* Attenuation functions */
double Attenuation01(double x, double s, double p);

/* for ODEs */
double odeintegrate(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status);
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
