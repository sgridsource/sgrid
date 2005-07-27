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

