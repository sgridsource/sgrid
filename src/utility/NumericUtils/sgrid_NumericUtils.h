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
int newton_linesrch_its(double x[], int n, int *check,
			void (*vecfunc)(int, double [], double []), 
			int MAXITS, double TOLF);
int newton_linesrch_itsP(double x[], int n, int *check,
			 void (*vecfuncP)(int, double [], double [], void *par),
			 void *par, int MAXITS, double TOLF);
int newton_lnsrch_vec_within_range(int n, double vec[],
                                   double hi[], double lo[]);
void newton_lnsrch_set_vecs_for_lininterp(int n, double vec[], 
          double hi[], double lo[], double *vb, double *vi);
void newton_lnsrch_get_fvec_by_lininterp(int n, double vec[], 
          double vb[], double vi[], double *fvec, double fvb[], double fvi[]);

/* 1D root finding with bracketing */
int zbrac_P(double (*func)(double,void *par), double *x1, double *x2, void *par);
int zbrent_itsP(double *x0, double (*func)(double,void *par),
                double x1, double x2, void *par, int ITMAX, double tol);

  	  	               	   
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

/* from integrals_1Dgrid.c */
double integrate_simpson_1Dgrid(double* f, double dx, int i1, int i2);
double integrate_trapez_1Dgrid(double* f, double dx, int i1, int i2);

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

/* 1D Minimization */
double brent_with_pointer_to_pars(double ax, double bx, double cx,
                                  double (*f)(double, void *ppointer),
                                  double tol, double *xmin, void *parpointer);
/* Multidim. Minimization */
void powell(double p[], double **xi, int n, double ftol, int *iter,
	    double *fret, double (*func)(double []));
