/* lu.c */
/* Crout's algorithm for LU decomposition and solve
   Wolfgang Tichy 6/2022 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* WT: do LU-decomp for matrix A[i][j] in place, here i,j \in [0,n-1]
   lu_decomp writes into A:
     U is in upper triangular part plus diagonal of A[i][j]
     L is in lower triangular part of A[i][j], with L_{ii}=1
   writes the permutations into p_idx
   writes the parity of the permutions in parity
   returns -(i+1), if i=row with all zeros
   returns +(j+1), if j=column with all zeros,
   returns 0 otherwise */
#define TINY 1.0e-300;
int lu_decomp(int n, double A[][n], int *p_idx, int *parity)
{
  int i,j,k;
  double big, tmp;
  double *vscal = malloc(sizeof(double)*n);

  if(!vscal)
  {
    printf("allocation failed!");
    abort();
    exit(1);
  }

  /* check if column with all zeros in a_ij */
  for(j=0; j<n; j++)
  {
    big = 0.;
    for(i=0; i<n; i++) /* loop over rows of col j */
      if((tmp=fabs(A[i][j])) > big) big = tmp;
    if(big == 0.) /* Singular matrix!!! */
    {
      /* A[i][j] is still unmodified here */
      free(vscal);
      return j+1; /* j=col with all zeros */
    }
  }

  /* check if row with all zeros in a_ij, and save scaling in vscal */
  for(i=0; i<n; i++)
  {
    big = 0.;
    for(j=0; j<n; j++)
      if((tmp=fabs(A[i][j])) > big) big = tmp;
    if(big == 0.) /* Singular matrix!!! */
    {
      /* A[i][j] is still unmodified here */
      free(vscal);
      return -(i+1); /* i=row with all zeros */
    }
    vscal[i] = 1./big;
  }

  /* set initial parity of the permutions */
  *parity=1;

  /* Crout's method */
  for(j=0; j<n; j++)
  {
    int imax;
    double sum, tmp2;

    for(i=0; i<j; i++)
    {
      sum = A[i][j];
      for(k=0; k<i; k++) sum -= A[i][k]*A[k][j];
      A[i][j] = sum;
    }

    big=0.;
    for(i=j; i<n; i++)
    {
      sum=A[i][j];
      for(k=0; k<j; k++) sum -= A[i][k]*A[k][j];
      A[i][j]=sum;
      if( (tmp2=vscal[i]*fabs(sum)) >= big )
      {
        big=tmp2;
        imax=i;
      }
    }
    if(j != imax)
    {
      for(k=0; k<n; k++)
      {
        tmp2 = A[imax][k];
        A[imax][k] = A[j][k];
        A[j][k] = tmp2;
      }
      *parity = -(*parity); /* record parity of the permutions */
      vscal[imax] = vscal[j];
    }
    p_idx[j] = imax; /* record permutions */
    if(A[j][j] == 0.) A[j][j] = TINY;
    if(j != n)
    {
      tmp2 = 1./(A[j][j]);
      for(i=j+1; i<n; i++) A[i][j] *= tmp2;
    }
  }
  free(vscal);
  return 0; /* return 0 if all is ok */
}
#undef TINY


/* Overwrite b with solution x = A^{-1} b
   here A and p_idx have to be set by a prior call to lu_decomp:
   I.e. A = LU plus A is permuted */
void lu_solve(int n, double A[][n], int *p_idx, double b[])
{
  int i,j, p_i;
  double sum;
  int i_nz;

  /* find y from L y = b, we write y into b */
  i_nz = -1; /* init to negative val. */
  for(i=0; i<n; i++)
  {
    p_i    = p_idx[i];  /* deal         */
    sum    = b[p_i];    /* with         */
    b[p_i] = b[i];      /* permutations */
    if(i_nz>=0)      for(j=i_nz; j<i; j++) sum -= A[i][j]*b[j];
    else if(sum!=0.) i_nz = i; /* set to i of first b[i]!=0 */
    b[i] = sum; /* Note: L_{ii}=1  */
  }

  /* find x from U x = y, (y is in b), we write x into b */
  for(i=n-1; i>=0; i--)
  {
    sum = b[i];
    for(j=i+1; j<n; j++) sum -= A[i][j]*b[j];
    b[i] = sum/A[i][i];
  }
}


/* calc det of LU decomposed A, i.e. A and parity have to be set
   by lu_decomp */
double lu_det(int n, double A[][n], int parity)
{
  double det=parity;
  int i;
  for(i=0; i<n; i++) det *= A[i][i];
  return det;
}


/* primitive print */
void pr_C_matrix(int n, double A[][n])
{
  int i,j;
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++) printf("%.15g  ", A[i][j]);
    printf("\n");
  }
}
void pr_C_vector(int n, double x[])
{
  int i;
  for(i=0; i<n; i++) printf("%.15g\n", i, x[i]);
}


/*
#define N 2
int main()
{
  double a_data[][N] = { {0.77949826, 0.56263848},
                         {0.91241331, 0.65857599} };
  double b[] = { 0.21685978, 0.25383732 };
  double (*A)[N];
  int parity, p_idx[N];
  
  A = a_data;

  pr_C_matrix(N,A);
  pr_C_vector(N,b);

  lu_decomp(N,A, p_idx, &parity);
  lu_solve(N,A, p_idx, b);

  pr_C_matrix(N,A);
  pr_C_vector(N,b);

  printf("det=%.15g\n", lu_det(N,A, parity));
}
#undef N
*/

/*
#define N 3
int main()
{
  double a_data[][N] = { {0.77949826, 0.56263848, 1.},
                         {0.91241331, 0.65857599, 2.},
                         {3.,         4.,         5.} };
  double b[] = { 0.21685978, 0.25383732, 3. };
  double (*A)[N];
  int parity, p_idx[N];

  A = a_data;

  pr_C_matrix(N,A);
  pr_C_vector(N,b);

  lu_decomp(N,A, p_idx, &parity);
  lu_solve(N,A, p_idx, b);

  pr_C_matrix(N,A);
  pr_C_vector(N,b);

  printf("det=%.19g\n", lu_det(N,A, parity));
}
*/
