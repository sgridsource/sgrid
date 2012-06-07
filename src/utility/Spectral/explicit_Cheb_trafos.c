/* explicit_Cheb_trafos.c */
/* does explicit slow Chebycheff trafos: */
/* Wolfgang Tichy 3/2004 */

#include "sgrid.h"
#include "Spectral.h"

/* define PI */
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876


/* Note here we have X = [a,b] and x=[1,-1]
   X = (a-b)*x/2 + (a+b)/2   so that: x=-1 => X=b,  x=1 => X=a
   x = (2X-b-a)/(a-b)
   This does not imply that the result of all Cheb. sums 
   (e.g. u[j] = sum_k c[k] T_k[x_j]
         c[j] ~ sum_k u[x_k] T_j[x_k] )
   gets a minus sign compared to what is done in numerical recipes, 
   even though T_k(-x) = (-1)^k T_k(x).
*/
/* NOT NEEDED: */
/* implementation of sign flip */
#define FlipSignOfOdd(c,n) { int j;\
        for(j=1; j<=(n); j+=2) (c)[j] = -(c)[j]; }


/* get Cheb coeffs c[0...n] of function f(X) (X in [a,b]) */
void cheb_getcoeff(double a, double b, double c[], int n, 
                   double (*func)(double))
{
  int k,j;
  int np1=n+1;
  double fac, bpa, amb, *f;

  f = (double*) calloc(np1,sizeof(double));
  
  amb=0.5*(a-b);
  bpa=0.5*(b+a);
  
  /* use function at zeros: */
  for(k=0;k<=n;k++)
  {
    double y=cos(PI*(k+0.5)/np1);
    f[k]=(*func)(y*amb+bpa);
  }
  fac=2.0/np1;
  for(j=0;j<=n;j++)
  {
    double sum=0.0;
    for (k=0;k<=n;k++)
      sum += f[k]*cos(PI*j*(k+0.5)/np1);
    c[j]=fac*sum;
  }
  free(f);
}


/* find value of function at X (in [a,b]) from Cheb coeffs c[0...n] */
double cheb_eval(double a, double b, double c[], int n, double X)
{
  double d=0.0, dd=0.0, sv, y, y2;
  int j;

  if ((X-a)*(X-b) > 0.0) printf("X not in range in routine cheb_eval\n");
  y=(2.0*X-a-b)/(a-b);
  y2=2.0*y;
  for(j=n;j>=1;j--)
  {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}


/* compute Cheb coeffs of deriv cder[0...n] from Cheb coeffs c[0...n] */
void cheb_deriv(void *aux,
                double a, double b, double c[], double cder[], int n)
{
  int j;
  double con;

  cder[n]=0.0;
  if(n>=1) cder[n-1]=2*(n)*c[n];
  for (j=n-2;j>=0;j--)
    cder[j]=cder[j+2]+2*(j+1)*c[j+1];

  /* convert to interval [a,b] */
  con=2.0/(a-b);
  for (j=0;j<n;j++)
    cder[j] *= con;
}


/* compute Cheb coeffs of integral cint[0...n] from Cheb coeffs c[0...n] */
void cheb_int(double a, double b, double c[], double cint[], int n)
{
  int j;
  double sum=0.0,fac=1.0,con;

  con=0.25*(a-b);
  for (j=1;j<=n-1;j++)
  {
    cint[j]=con*(c[j-1]-c[j+1])/j;
    sum += fac*cint[j];
    fac = -fac;
  }
  cint[n]=con*c[n-1]/(n);
  sum += fac*cint[n];
  cint[0]=2.0*sum;  /* <--arbitrary const picked as in numrec */
}


/* compute Cheb coeffs c[0...n] from function u at the zeros of T_N(x).
   Note N=n+1                                                           */
void cheb_coeffs_fromZeros(double c[], double u[], int n)
{
  int k, j;
  int N=n+1;
  double fac, sum, PIoN;

  PIoN=PI/N;
  fac=2.0/N;
  
  for(j=0;j<=n;j++)
  {
    sum=0.0;
    for(k=0;k<=n;k++)
      sum += cos(j*PIoN*(k+0.5))*u[k];
    c[j]=fac*sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* compute Cheb coeffs c[0...N] from function u at the extrema of T_N(x). */
void cheb_coeffs_fromExtrema(double c[], double u[], int N)
{
  int k, j;
  double fac, sum, PIoN;

  PIoN=PI/N;
  fac=2.0/N;
  
  for(j=0;j<=N;j++)
  {
    sum = 0.5 * (u[0] + cos(j*PI)*u[N]);
    for(k=1;k<N;k++)
      sum += cos(j*PIoN*k)*u[k];
    c[j]=fac*sum;
  }
  c[N]*=0.5;
  /* we should use a FFT for eveything above this line */
}


/* find function u on the zeros of T_N(x).   Note N=n+1 */
void cheb_eval_onZeros(double c[], double u[], int n)
{
  int k, j;
  int N=n+1;
  double fac, sum, PIoN;

  PIoN=PI/N;
  fac=2.0/N;
  
  for(j=0;j<=n;j++)
  {
    sum=0.5*c[0];
    for(k=1;k<=n;k++)
      sum += c[k]*cos(k*PIoN*(j+0.5));
    u[j]=sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* find function u on the extrema of T_N(X) */
void cheb_eval_onExtrema(double c[], double u[], int N)
{
  int k, j;
  double fac, sum, PIoN;

  PIoN=PI/N;
  fac=2.0/N;
  
  for(j=0;j<=N;j++)
  {
    sum = 0.5*c[0];
    for(k=1;k<=N;k++)
      sum += c[k]*cos(k*PIoN*j);
    u[j]=sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* filter: zero all c[j] with k<=j<=n */
void cheb_filter(double c[], int k, int n)
{
  int j;
  
  for(j=k;j<=n;j++) c[j]=0;
}


/* find value of Cheb. basis function T_n at X (in [a,b]) */
double cheb_basisfunc_FromSum(void *aux, double a, double b, int n, int n1, double X)
{
  double d=0.0, dd=0.0, sv, y, y2;
  int j;

  if(n==0) return 0.5;  /* in numrec T_0 / 2 is used as basisfunc # 0 */

  if ((X-a)*(X-b) > 0.0) printf("X not in range in routine cheb_basis\n");
  y=(2.0*X-a-b)/(a-b);
  y2=2.0*y;

  sv=d;
  d=y2*d-dd+ 1.0;
  dd=sv;
        
  for(j=n-1;j>=1;j--)
  {
    sv=d;
    d=y2*d-dd;
    dd=sv;
  }
  return y*d-dd;
}


/* find value of Cheb. basis function T_n at X (in [a,b]) */
double cheb_basisfunc(void *aux, double a, double b, int n, int n1, double X)
{
  double y;
  if(n==0) return 0.5;  /* in numrec T_0 / 2 is used as basisfunc # 0 */

  /* y=(2.0*X-a-b)/(a-b); */
  y=((X-a)+(X-b))/(a-b);
  return cos(n*acos(y));
}
