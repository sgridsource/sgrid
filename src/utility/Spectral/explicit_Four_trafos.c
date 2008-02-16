/* explicit_Four_trafos.c */
/* does explicit slow Fourier trafos: */
/* Wolfgang Tichy 4/2004 */

#include "sgrid.h"
#include "Spectral.h"

/* define PI */
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876


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

  for(j=1; 2*j<N; j++)
  {
    cder[2*j-1] =  c[2*j] * PI2_con*j;
    cder[2*j]   = -c[2*j-1] * PI2_con*j;
  }
  cder[0] = 0.0;
  if( N%2 == 0 ) cder[N-1] = 0.0;
  
  //else  /* this is what I had in the old explicit_Four_trafos.c version */
  //{ cder[N-2] = 0.0;  cder[N-1] = 0.0; }
}


/* compute Cheb coeffs of integral cint[0...n] from Cheb coeffs c[0...n] */
void four_int(double a, double b, double c[], double cint[], int n)
{
  int j;
  double PI2_con, L;
  int N=n+1;

  L = b-a;
  PI2_con = 2.0*PI/L;

  /* get terms coming from integrating c[0] */
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] = -0.5*L*c[0]/((double) N);
    cint[2*j]   = -0.5*L*c[0]/((double) N); // WRONG!!!!!
    // integrate the func 1/N and find its coeffs instead!!!
    // multiply this by c[0]
    if(N!=4) errorexit("four_int is wrong");
  }
  if( N%2 == 0 ) cint[N-1] = -0.5*L*c[0]/((double) N);

  cint[0] = 0.5*n*L*c[0]/((double) N);  

  /* add terms coming from integrating everything but the c[0] term */
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] += -c[2*j] / (PI2_con*j);
    cint[2*j]   +=  c[2*j-1] / (PI2_con*j);
  }
  if( N%2 == 0 ) cint[N-1] += 0.0;
}


/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 */
void four_coeffs_alt(double c[], double u[], int n)
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


/* find function u from Four coeffs c[0...n], computed with four_coeffs_alt */
void four_eval_alt(double c[], double u[], int n)
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
NOTE: four_coeffs returns c[] that are N times of those of four_coeffs_alt */
void four_coeffs(double c[], double u[], int n)
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
double four_basisfunc_alt(void *aux, double a, double b, int n, int N, double X)
{
  double K = 2.0*PI/(b-a);
  int j = n/2 + n%2;

  if(n%2!=0) return cos(j*K*X);
  if(n==0)   return 1.0;
  return sin(j*K*X);
}


/* find value of Fourier basis function B_k at X (in [a,b]) */
double four_basisfunc(void *aux, double a, double b, int k, int N, double X)
{
  double K = 2.0*PI/(b-a);
  int j = k/2 + k%2;

  if(k==0)              return 1.0/N;
  if(k==N-1 && N%2==0)  return cos(j*K*X)*1.0/N;;
  if(k%2!=0)            return cos(j*K*X)*2.0/N;
  return sin(j*K*X)*2.0/N;
}
