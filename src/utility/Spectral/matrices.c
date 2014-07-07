/* matrices.c */
/* Wolfgang Tichy 7/2005 */

#include "sgrid.h"
#include "Spectral.h"
#ifdef CBLAS
#include <cblas.h>
#endif


/* do we use BLAS??? */
#ifdef CBLAS
/* this is just a matrix multiplication Mu = M u */
void matrix_times_vector(double *M, double *u, double *Mu, int n)
{
  /* call dgemv (matrix times vector) from CBLAS: */
  /* Note: we need something like this in MyConfig if we use CBLAS from ATLAS:
     DFLAGS += -DCBLAS 
     CBLASDIR = /home/wolf/Packages/ATLAS3.8.3/my_build_dir/
     #SPECIALINCS += -I$(CBLASDIR)/include
     SPECIALLIBS += -L$(CBLASDIR)/lib -lcblas -latlas   */  
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, M, n,
              u, 1, 0.0, Mu, 1);
/* this is what we need if we call the Fortran DGEMV from original
   Fortan BLAS: */
/*
  char trans[] = "T";
  double done = 1.0;
  double dzero= 0.0;
  int one = 1;
  dgemv_(trans, &n, &n, &done, M, &n, u, &one, &dzero, Mu, &one);
*/
}
#else
/* this is just a matrix multiplication Mu = M u */
void matrix_times_vector(double *M, double *u, double *Mu, int n)
{
  int i;
  double *M_i = M; /* we'll increment M_i so that matrix M_i[j] = M_{ij} */

  for(i=0; i<n; i++, M_i+=n)
  {
    int j;
    double sum=0.0;
    /* for(j=0; j<n; j++)  sum += M[n*i + j] * u[j]; */
    for(j=0; j<n; j++)  sum += M_i[j] * u[j];
    Mu[i] = sum;
  }
}
#endif

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
