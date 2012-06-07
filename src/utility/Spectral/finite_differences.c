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
void fd2_deriv_periodic(void *aux,
                        double a, double b, double c[], double cder[], int n)
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
void fd2_deriv_onesidedBC(void *aux,
                          double a, double b, double c[], double cder[], int n)
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
void fd2_2ndderiv_periodic(void *aux, double a, double b, 
                           double c[], double cder[], int n)
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
void fd2_2ndderiv_onesidedBC(void *aux, double a, double b, 
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
    cder[j] = (c[j+1] - c[j-1]) / (x[j+1] - x[j-1]);

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


/* basis functions for Kronecker delta basis */ 
/* basis function in direction1 for Kronecker delta basis */
double fd_basis1(void *aux, double a, double b, int k, int N, double X)
{
  tBox *box = (tBox *) aux;

  if(box==NULL)
  {
    if(X<=a+dequaleps)      return (k == 0);
    else if(X>=b-dequaleps) return (k == N-1);
    else                    return (k == (int) (N*(X-a)/(b-a)));
  }
  else
  {
    int n1 = box->n1;
    int n2 = box->n2;
    /* int n3 = box->n3; */
    int iX;
    int indX=Ind("X");
    double XvarL, XvarR, Xd;

    for(iX=0; iX<N-1; iX++)
    {
      XvarL = box->v[indX][Index(iX,0,0)];
      XvarR = box->v[indX][Index(iX+1,0,0)];
      Xd = (XvarL + XvarR)*0.5;
      if(X<Xd) return k==iX;
    }
    return k==N-1;
  }
}
/* basis function in direction2 for Kronecker delta basis */
double fd_basis2(void *aux, double a, double b, int k, int N, double Y)
{
  tBox *box = (tBox *) aux;

  if(box==NULL)
  {
    if(Y<=a+dequaleps)      return (k == 0);
    else if(Y>=b-dequaleps) return (k == N-1);
    else                    return (k == (int) (N*(Y-a)/(b-a)));
  }
  else
  {
    int n1 = box->n1;
    int n2 = box->n2;
    /* int n3 = box->n3; */
    int iY;
    int indY=Ind("Y");
    double YvarL, YvarR, Yd;

    for(iY=0; iY<N-1; iY++)
    {
      YvarL = box->v[indY][Index(0,iY,0)];
      YvarR = box->v[indY][Index(0,iY+1,0)];
      Yd = (YvarL + YvarR)*0.5;
      if(Y<Yd) return k==iY;
    }
    return k==N-1;
  }
}
/* basis function in direction3 for Kronecker delta basis */
double fd_basis3(void *aux, double a, double b, int k, int N, double Z)
{
  tBox *box = (tBox *) aux;

  if(box==NULL)
  {
    if(Z<=a+dequaleps)      return (k == 0);
    else if(Z>=b-dequaleps) return (k == N-1);
    else                    return (k == (int) (N*(Z-a)/(b-a)));
  }
  else
  {
    int n1 = box->n1;
    int n2 = box->n2;
    /* int n3 = box->n3; */
    int iZ;
    int indZ=Ind("Z");
    double ZvarL, ZvarR, Zd;

    for(iZ=0; iZ<N-1; iZ++)
    {
      ZvarL = box->v[indZ][Index(0,0,iZ)];
      ZvarR = box->v[indZ][Index(0,0,iZ+1)];
      Zd = (ZvarL + ZvarR)*0.5;
      if(Z<Zd) return k==iZ;
    }
    return k==N-1;
  }
}
