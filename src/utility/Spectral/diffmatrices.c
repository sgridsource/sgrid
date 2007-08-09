/* diffmatrices.c */
/* Wolfgang Tichy 3/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* init the n1*n1 diff. matrices D and DD */
void initdiffmatrix(double a, double b, double *D, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_deriv)(double, double, double *,double *, int),
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
    
    //cheb_deriv(a, b,  c, d, n1-1);
    coeffs_of_deriv(a, b,  c, d, n1-1);
    
    //cheb_eval_onExtrema(d, c, n1-1);
    eval_onPoints(d, c, n1-1);
    
    /* set D */
    for(i=0; i<n1; i++) D[n1*i + j] = c[i];
    
    //cheb_deriv(a, b,  d, c, n1-1);
    coeffs_of_deriv(a, b,  d, c, n1-1);

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


/* init the diff. matrix DD for second derivs only */
void initdiffmatrix2(double a, double b, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_2ndderiv)(double, double, double *,double *, int),
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
    get_coeffs(c, u, n1-1);
    coeffs_of_2ndderiv(a, b,  c, d, n1-1);
    eval_onPoints(d, c, n1-1);
    
    /* set DD */
    for(i=0; i<n1; i++) DD[n1*i + j] = c[i];
    
    u[j]=0.0;
  }

  free(u);
  free(c);
  free(d);
}


/* take 1d deriv of u along a line */
/* this is just a matrix multiplication */
/*
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
*/


/* init the special n1*n1 diff. matrices D and DD using finite differencing
   for a grid with non-uniform grid spacing */
void init_fdcentered_diffmatrix(double *x, double *D, int n1,
              void (*fd_deriv)(double *, double *,double *, int) )
{ /* fd_deriv=fdcentered_deriv_onesidedBC,fdDp_deriv_onesidedBC,fdDm_deriv_onesidedBC */
  int i,j;
  double *u;
  double *c;
  double *d;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));
  d = (double *) calloc(n1, sizeof(double));

  if( !(u && c && d) )
    errorexit("init_fdcentered_diffmatrix: out of memory for u, c, d");

  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    u[j]=1.0;

    /* use centered first derivs */
    fd2_coeffs(c, u, n1-1);
    // fdcentered_deriv_onesidedBC(x, c, d, n1-1);
    // or: fdDp_deriv_onesidedBC(x, c, d, n1-1);
    // or: fdDm_deriv_onesidedBC(x, c, d, n1-1);
    fd_deriv(x, c, d, n1-1);
    fd2_eval(d, c, n1-1);

    /* set D */
    for(i=0; i<n1; i++) D[n1*i + j] = c[i];

/*
    fd2_coeffs(c, u, n1-1);
    fdDp_deriv_onesidedBC(x, c, d, n1-1);
    fdDm_deriv_onesidedBC(x, d, c, n1-1);
    // ^Note: If I switch fdDp_deriv_onesidedBC and fdDm_deriv_onesidedBC 
    //        c will be different at this point and thus DD wil change 
    fd2_eval(c, d, n1-1);

    // set DD 
    for(i=0; i<n1; i++) DD[n1*i + j] = d[i];
*/
    u[j]=0.0;
  }

  free(u);
  free(c);
  free(d);
}
