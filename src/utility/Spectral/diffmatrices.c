/* diffmatrices.c */
/* Wolfgang Tichy 3/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* init a n1*n1 diff. matrix */
void initdiffmatrix(double a, double b, double *D, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) )
{
  int i,j;
  double *u;
  double *c;
  double *d;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));
  d = (double *) calloc(n1, sizeof(double));

  if( !(u && c && d) ) errorexit("initdiffmatrix: out of memory for u, c, d");


  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    u[j]=1.0;
    //cheb_coeffs_fromExtrema(c, u, n1-1);
    get_coeffs(c, u, n1-1);
    
    cheb_deriv(a, b,  c, d, n1-1);
    
    //cheb_eval_onExtrema(d, c, n1-1);
    eval_onPoints(d, c, n1-1);
    
    /* set D */
    for(i=0; i<n1; i++) D[n1*i + j] = c[i];
    
    cheb_deriv(a, b,  d, c, n1-1);

    //cheb_eval_onExtrema(c, d, n1-1);
    eval_onPoints(c, d, n1-1);

    /* set DD */
    for(i=0; i<n1; i++) DD[n1*i + j] = d[i];

    u[j]=0.0;
  }

  free(u);
  free(c);
  free(d);
}


/* take 1d deriv of u along a line */
/* this is just a matrix multiplication */
void diffmat_deriv(double *D, double *u, double *du, int n)
{
  int i,j;
  double sum;
  
  for(i=0; i<n; i++)
  {
    sum=0.0;
    for(j=0; j<n; j++)  sum += D[n*i + j] * u[j];
    du[i] = sum;
  }
}

