/* Spectral.h */
/* (c) Wolfgang Tichy 1.4.2005 */
/* header file for Spectral functions */


/* Functions from explicit_Cheb_trafos.c */

/* get Cheb coeffs c[0...n] of function f(X) (X in [a,b]) */
void cheb_getcoeff(double a, double b, double c[], int n, 
                   double (*func)(double));
/* find value of function at X (in [a,b]) from from Cheb coeffs c[0...n]  */
double cheb_eval(double a, double b, double c[], int n, double X);

/* compute Cheb coeffs of deriv cder[0...n] from Cheb coeffs c[0...n] */
void cheb_deriv(double a, double b, double c[], double cder[], int n);

/* compute Cheb coeffs c[0...n] from function u at the zeros of T_N(x).
   Note N=n+1                                                           */
void cheb_coeffs_fromZeros(double c[], double u[], int n);

/* compute Cheb coeffs c[0...N] from function u at the extrema of T_N(x). */
void cheb_coeffs_fromExtrema(double c[], double u[], int N);

/* find function u on the zeros of T_N(x).   Note N=n+1 */
void cheb_eval_onZeros(double c[], double u[], int n);

/* find function u on the extrema of T_N(X) */
void cheb_eval_onExtrema(double c[], double u[], int N);

/* filter: zero all c[j] with k<=j<=n */
void cheb_filter(double c[], int k, int n);


/* Functions from ... */
