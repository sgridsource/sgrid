/* finite_differences.c */
/* does finite differencing using a Kronecker delta basis:
   u(x_j) = sum_k c[k] f_k(x_j)    c_k = sum_j u(x_j) f_j(x_k)
   with f_k(x_j) = delta_kj
   Then \nabla f_k(x_j) = ( f_{k-1}(x_j) - f_{k+1}(x_j) )/(2h)
   and  c'[k] = (c[k+1] - c[k-1])/(2h)    */
/* Wolfgang Tichy 5/2005 */

#include "sgrid.h"
#include "Spectral.h"


/* compute coeffs of deriv cder[0...n] from coeffs c for a periodic grid:
 x_j = j (b-a)/(n+1) + a ,  j=0, ..., n */
void fd2_deriv_periodic(double a, double b, double c[], double cder[], int n)
{
  int j;
  double con = 0.5*(n+1)/(b-a);

  for (j=1;j<n;j++)
    cder[j] =  (c[j+1] - c[j-1]) * con;

  cder[0] = (c[1] - c[n]) * con;
  cder[n] = (c[0] - c[n-1]) * con;
}


/* compute coeffs of deriv cder[0...n] from coeffs c for a non-periodic grid:
 x_j = j (b-a)/n + a ,  j=0, ..., n */
void fd2_deriv_onesidedBC(double a, double b, double c[], double cder[], int n)
{
  int j;
  double con = 0.5*n/(b-a);

  for (j=1;j<n;j++)
    cder[j] =  (c[j+1] - c[j-1]) * con;

  cder[0] = (-3.0*c[0] + 4.0*c[1]   - c[2]) * con;
  cder[n] = ( 3.0*c[n] - 4.0*c[n-1] + c[n-2]) * con;
}


/* compute coeffs of 2nd deriv cder[0...n] from coeffs c for a periodic grid:
 x_j = j (b-a)/(n+1) + a ,  j=0, ..., n */
void fd2_2ndderiv_periodic(double a, double b, double c[], double cder[], int n)
{
  int j;
  double ooh = (n+1)/(b-a);
  double con = ooh*ooh;

  for (j=1;j<n;j++)
    cder[j] =  (c[j+1] - 2.0*c[j] + c[j-1]) * con;

  cder[0] = (c[1] - 2.0*c[0] + c[n]) * con;
  cder[n] = (c[0] - 2.0*c[n] + c[n-1]) * con;
}


/* compute coeffs of 2nd deriv cder[0...n] from coeffs c for a non-periodic 
   grid: x_j = j (b-a)/n + a ,  j=0, ..., n */
void fd2_2ndderiv_onesidedBC(double a, double b, 
                             double c[], double cder[], int n)
{
  int j;
  double ooh = n/(b-a);
  double con = ooh*ooh;

  for (j=1;j<n;j++)
    cder[j] =  (c[j+1] - 2.0*c[j] + c[j-1]) * con;

  cder[0] = ( 2.0*c[0] - 5.0*c[1]   + 4.0*c[2]   - c[3]) * con;
  cder[n] = ( 2.0*c[n] - 5.0*c[n-1] + 4.0*c[n-2] - c[n-3]) * con;
}


/* compute Four coeffs c[0...n] from function u */
void fd2_coeffs(double c[], double u[], int n)
{
  int k;

  for(k=0;k<=n;k++)  c[k] = u[k];
}


/* find function u from Four coeffs c[0...n] */
void fd2_eval(double c[], double u[], int n)
{
  int k;
  
  for(k=0;k<=n;k++)  u[k] = c[k];
}


/* filter: zero all c[j] with k<=j<=n */
void fd2_filter(double c[], int k, int n)
{
  /* not yet done !*/
  return ;
}


/* compute coeffs cder[0...n] of centered deriv from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdcentered_deriv_onesidedBC(double x[], double c[], double cder[], int n)
{
  int j;

  for (j=1;j<n;j++)
    cder[j] =  0.5 * (c[j+1] - c[j-1]) / (x[j+1] - x[j-1]);

  cder[0] = (c[1] - c[0])  /(x[1] - x[0]);
  cder[n] = (c[n] - c[n-1])/(x[n] - x[n-1]);
}


/* compute coeffs cder[0...n] of deriv D^+ from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdDp_deriv_onesidedBC(double x[], double c[], double cder[], int n)
{
  int j;

  for (j=0;j<n;j++)
    cder[j] =  (c[j+1] - c[j]) / (x[j+1] - x[j]);

  cder[n] = (c[n] - c[n-1])/(x[n] - x[n-1]);
}


/* compute coeffs cder[0...n] of deriv D^- from coeffs c[0...n] for
   a non-periodic grid with non-uniform grid spacing: x_j,  j=0, ..., n */
void fdDm_deriv_onesidedBC(double x[], double c[], double cder[], int n)
{
  int j;

  for (j=1;j<=n;j++)
    cder[j] =  (c[j] - c[j-1]) / (x[j] - x[j-1]);

  cder[0] = (c[1] - c[0])/(x[1] - x[0]);
}
