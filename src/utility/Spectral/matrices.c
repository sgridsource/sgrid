/* matrices.c */
/* Wolfgang Tichy 7/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* this is just a matrix multiplication Mu = M u */
void matrix_times_vector(double *M, double *u, double *Mu, int n)
{
  int i,j;
  double sum;
  
  for(i=0; i<n; i++)
  {
    sum=0.0;
    for(j=0; j<n; j++)  sum += M[n*i + j] * u[j];
    Mu[i] = sum;
  }
}

/* multiply two matricies M and D: MD = M D ,   MD_ik = M_ij D_jk */
void matrix_times_matrix(double *M, double *D, double *MD, int n)
{
  int i,j,k;
  double sum;
  
  for(i=0; i<n; i++)
    for(k=0; k<n; k++)
    {
      sum=0.0;
      for(j=0; j<n; j++)  sum += M[n*i + j] * D[n*j + k];
      MD[n*i + k] = sum;
    }
}

/* multiply vector transpose with matrix uM = u^T M ,  uM_j = u_i M_ij */
void vector_times_matrix(double *u, double *M, double *uM, int n)
{
  int i,j;
  double sum;
  
  for(j=0; j<n; j++)
  {
    sum=0.0;
    for(i=0; i<n; i++)  sum += u[i] * M[n*i + j];
    uM[j] = sum;
  }
}

/* scalarproduct of 2 vectors */
double scalarproduct_vectors(double *v, double *w, int n)
{
  int i;
  double sum=0.0;
  for(i=0; i<n; i++)  sum += v[i]*w[i];
  return sum;
}
