/* NumericUtils.h */
/* (c) Wolfgang Tichy 6.6.2003 */
/* header file for NumericUtils functions */


/* Functions */

/* for Newton-Raphson with line searches */
void newton_lnsrch(double x[], int n, int *check,
                   void (*vecfunc)(int, double [], double []),
               	   int MAXITS, double TOLF);
void fd_jacobian(int n, double x[], double fvec[], double **df,
  	void (*vecfunc)(int, double [], double []));
  	  	               	   
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
