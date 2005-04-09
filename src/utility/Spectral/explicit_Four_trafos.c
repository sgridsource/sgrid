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

  for (j=1;j<N/2;j++)
  {
    cder[2*j]   = c[2*j+1] * PI2_con*j;
    cder[2*j+1] = -c[2*j] * PI2_con*j;
  }
  cder[0] = cder[1] = 0.0;
}



/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 */
void four_coeffs(double c[], double u[], int n)
{
  int k, j;
  double Re_c_j, Im_c_j, PI2oN;
  int N=n+1;

  PI2oN=2.0*PI/N;
  
  for(j=0; j<N/2; j++)
  {
    Re_c_j = Im_c_j = 0.0;
    for(k=0;k<N;k++)
    {
      Re_c_j += cos(j*PI2oN*k)*u[k];
      Im_c_j += sin(j*PI2oN*k)*u[k];
    }
    c[2*j]  = Re_c_j; /* real part of c_j */
    c[2*j+1]= Im_c_j; /* imaginary prat of c_j */
  }
  /* Now j=N/2 */
  j=N/2;
  Re_c_j = Im_c_j = 0.0;
  for(k=0;k<N;k++)
  {
    /*
    Re_c_j += cos(j*PI2oN*k)*u[k];
    Im_c_j += sin(j*PI2oN*k)*u[k]; // Note Im_c_j will be automatically zero
    */
    Re_c_j += cos(PI*k)*u[k];
  }
  c[1] = Re_c_j; /* Note we save real part of c_N/2 in c[1] */
  /* we should use a FFT for eveything above this line */
}



/* find function u from Four coeffs c[0...n] */
void four_eval(double c[], double u[], int n)
{
  int k, j;
  double sum, Re_c_k, Im_c_k, PI2oN;
  int N=n+1;

  PI2oN=2.0*PI/N;
  
  for(j=0; j<N; j++)
  {
    sum = 0.0;
    for(k=1;k<=N/2-1;k++)
    {
      Re_c_k = c[2*k];
      Im_c_k = c[2*k+1];
      sum += cos(k*PI2oN*j)*Re_c_k + sin(k*PI2oN*j)*Im_c_k;
    }
    u[j] = 2.0*sum + c[0] + c[1]*cos(PI*j);
  }
  /* we should use a FFT for eveything above this line */
}


/* filter: zero all c[j] with k<=j<=n */
void four_filter(double c[], int k, int n)
{
  int j;
  int N=n+1;

  c[1]=0;  
  for(j=k/2;j<N/2;j++)
  {
    c[2*j] = c[2*j+1] = 0.0;
  }
}
