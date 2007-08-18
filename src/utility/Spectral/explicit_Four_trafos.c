/* explicit_Four_trafos.c */
/* does explicit slow Fourier trafos: */
/* Wolfgang Tichy 4/2004 */

#include "sgrid.h"
#include "Spectral.h"

/* define PI */
#define PI  3.1415926535897932
#define PIh 1.5707963267948966
#define PIq 0.78539816339744831


/* Note here we have X = [a,b] and x=[0,1]
   X = (b-a)*x + a   so that: x=0 => X=a,  x=1 => X=b
   x = (X-a)/(b-a)
*/



/* compute Four coeffs of deriv cder[0...n] from Four coeffs c[0...n] */
void four_deriv(double a, double b, double c[], double cder[], int n)
{
  int j;
  double PI2_con;
  int N=n+1;

  PI2_con = 2.0*PI/(b-a);

  for (j=1;j<=N/2;j++)
  {
    cder[2*j-1] =  c[2*j] * PI2_con*j;
    if(2*j<N)
      cder[2*j]   = -c[2*j-1] * PI2_con*j;
  }
  cder[0] = 0.0;
  if( N%2 == 0 ) cder[N-1] = 0.0;
  
  //else  /* this is what I had in the old explicit_Four_trafos.c version */
  //{ cder[N-2] = 0.0;  cder[N-1] = 0.0; }
}


/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 */
void four_coeffs(double c[], double u[], int n)
{
  int k, j;
  double Re_c_j, Im_c_j, PI2oN, toN;
  int N=n+1;

  toN=2.0/N;
  PI2oN=PI*toN;

  /* first j=0 , i.e. c_0 */
  c[0] = 0.0;
  for(k=0;k<N;k++)  c[0] += u[k];
  c[0] = c[0]/N;

  for(j=1; j<=N/2; j++)
  {
    Re_c_j = Im_c_j = 0.0;
    for(k=0;k<N;k++)
    {
      Re_c_j += cos(j*PI2oN*k)*u[k];
      Im_c_j += sin(j*PI2oN*k)*u[k];
    }
    c[2*j-1] = Re_c_j * toN; /* real part of c_j */
    if(2*j<N)
      c[2*j] = Im_c_j * toN; /* imaginary part of c_j */
  }
  if( N%2 == 0 )
    c[N-1] = 0.5*c[N-1];
  /* we should use a FFT for everything above this line */
}


/* find function u from Four coeffs c[0...n], computed with four_coeffs */
void four_eval(double c[], double u[], int n)
{
  int k, j;
  double sum, Re_c_k, Im_c_k, PI2oN;
  int N=n+1;

  PI2oN=2.0*PI/N;
  
  for(j=0; j<N; j++)
  {
    sum = 0.0;
    for(k=1;k<=N/2;k++)
    {
      Re_c_k   = c[2*k-1];
      if(2*k<N)
        Im_c_k = c[2*k];
      else
        Im_c_k = 0.0;
      sum += cos(k*PI2oN*j)*Re_c_k + sin(k*PI2oN*j)*Im_c_k;
    }
    if( N%2 == 0 )
      u[j] = ( c[0] + sum - c[N-1]*cos(PI*j) );
    else
      u[j] = ( c[0] + sum );
  }
  /* we should use a FFT for eveything above this line */
}


/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 
NOTE: four_coeffsN returns c[] that are N times of those of four_coeffs */
void four_coeffsN(double c[], double u[], int n)
{
  int k, j;
  double Re_c_j, Im_c_j, PI2oN;
  int N=n+1;

  PI2oN=2.0*PI/N;

  /* first j=0 , i.e. c_0 */
  c[0] = 0.0;
  for(k=0;k<N;k++)  c[0] += u[k];

  for(j=1; j<=N/2; j++)
  {
    Re_c_j = Im_c_j = 0.0;
    for(k=0;k<N;k++)
    {
      Re_c_j += cos(j*PI2oN*k)*u[k];
      Im_c_j += sin(j*PI2oN*k)*u[k];
    }
    c[2*j-1] = Re_c_j; /* real part of c_j */
    if(2*j<N)
      c[2*j] = Im_c_j; /* imaginary part of c_j */
  }
  /* we should use a FFT for everything above this line */
}


/* find function u from Four coeffs c[0...n], computed with four_coeffsN */
void four_evalN(double c[], double u[], int n)
{
  int k, j;
  double sum, Re_c_k, Im_c_k, PI2oN;
  int N=n+1;

  PI2oN=2.0*PI/N;
  
  for(j=0; j<N; j++)
  {
    sum = 0.0;
    for(k=1;k<=N/2;k++)
    {
      Re_c_k   = c[2*k-1];
      if(2*k<N)
        Im_c_k = c[2*k];
      else
        Im_c_k = 0.0;
      sum += cos(k*PI2oN*j)*Re_c_k + sin(k*PI2oN*j)*Im_c_k;
    }
    if( N%2 == 0 )
      u[j] = ( c[0] + 2.0*sum - c[N-1]*cos(PI*j) )/N;
    else
      u[j] = ( c[0] + 2.0*sum )/N;
  }
  /* we should use a FFT for eveything above this line */
}


/* filter: zero all c[j] with k<=j<=n */
void four_filter(double c[], int k, int n)
{
  int j;

  for(j=k;j<=n;j++) c[j]=0;
}


/* find value of Fourier basis function B_n at X (in [a,b]) */
double four_basisfunc(double a, double b, int n, double X)
{
  double K = 2.0*PI/(b-a);
  int j = n/2 + n%2;

  if(n%2!=0) return cos(j*K*X);
  if(n==0)   return 1.0;
  return sin(j*K*X);
}
