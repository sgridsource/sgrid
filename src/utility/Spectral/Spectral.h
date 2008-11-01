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

/* compute Cheb coeffs of integral cint[0...n] from Cheb coeffs c[0...n] */
void cheb_int(double a, double b, double c[], double cint[], int n);

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

/* find value of Cheb. basis function T_n at X (in [a,b]) */
double cheb_basisfunc(void *aux, double a, double b, int n, int n1, double X);

/* find value of Cheb. basis function T_n at X (in [a,b]) */
double cheb_basisfunc_FromSum(void *aux, double a, double b, int n, int n1, double X);


/* Functions from explicit_Four_trafos.c */

/* compute Four coeffs of deriv cder[0...N-1] from Four coeffs c[0...N-1] */
void four_deriv(double a, double b, double c[], double cder[], int N);

/* compute Cheb coeffs of integral cint[0...n] from Cheb coeffs c[0...n] */
void four_int(double a, double b, double c[], double cint[], int n);

/* compute Four coeffs c[0...N-1] from function u at x_k = k/N, k=0,...,N-1 */
void four_coeffs(double c[], double u[], int N);

/* find function u from Four coeffs c[0...N-1] */
void four_eval(double c[], double u[], int N);

/* filter: zero all c[j] with k<=j<=n */
void four_filter(double c[], int k, int N);

/* find value of Fourier basis function B_n at X (in [a,b]) */
double four_basisfunc(void *aux, double a, double b, int n, int N, double X);


/* Functions from finite_differences.c */

/* compute coeffs of deriv cder[0...n] from coeffs c for a periodic grid:
 x_j = j (b-a)/(n+1) + a ,  j=0, ..., n */
void fd2_deriv_periodic(double a, double b, double c[], double cder[], int n);

/* compute coeffs of deriv cder[0...n] from coeffs c for a non-periodic grid:
 x_j = j (b-a)/n + a ,  j=0, ..., n */
void fd2_deriv_onesidedBC(double a, double b, double c[], double cder[], int n);

/* compute coeffs of 2nd deriv cder[0...n] from coeffs c for a periodic grid:
 x_j = j (b-a)/(n+1) + a ,  j=0, ..., n */
void fd2_2ndderiv_periodic(double a, double b, double c[], double cder[], int n);

/* compute coeffs of 2nd deriv cder[0...n] from coeffs c for a non-periodic 
   grid: x_j = j (b-a)/n + a ,  j=0, ..., n */
void fd2_2ndderiv_onesidedBC(double a, double b, 
                             double c[], double cder[], int n);

/* compute Four coeffs c[0...n] from function u */
void fd2_coeffs(double c[], double u[], int n);

/* find function u from Four coeffs c[0...n] */
void fd2_eval(double c[], double u[], int n);

/* filter: zero all c[j] with k<=j<=n */
void fd2_filter(double c[], int k, int n);

/* compute coeffs cder[0...n] of centered deriv from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdcentered_deriv_onesidedBC(double x[], double c[], double cder[], int n);

/* compute coeffs cder[0...n] of deriv D^+ from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdDp_deriv_onesidedBC(double x[], double c[], double cder[], int n);

/* compute coeffs cder[0...n] of deriv D^- from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdDm_deriv_onesidedBC(double x[], double c[], double cder[], int n);

/* basis functions for Kronecker delta basis */
double fd_basis1(void *aux, double a, double b, int k, int N, double X);
double fd_basis2(void *aux, double a, double b, int k, int N, double Y);
double fd_basis3(void *aux, double a, double b, int k, int N, double Z);


/* Functions from FFTs_for_sgrid.c */
void four_coeffs_numrecFFT(double c[], double u[], int n);
void four_eval_numrecFFT(double c[], double u[], int n);
void four_int_numrecFFT(double a, double b, double c[], double cint[], int n);
void cheb_coeffs_fromZeros_numrecFFT(double c[], double u[], int n);
void cheb_coeffs_fromExtrema_numrecFFT(double c[], double u[], int n);
void cheb_eval_onZeros_numrecFFT(double c[], double u[], int n);
void cheb_eval_onExtrema_numrecFFT(double c[], double u[], int N);
void four_coeffs_FFTW3(double *c, double *u, int n);
void four_eval_FFTW3(double *c, double *u, int n);
