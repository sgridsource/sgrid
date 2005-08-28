/* spec_coeffs.c */
/* Wolfgang Tichy 7/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* init a n1*n1 matrix M used to compute coeffs */
void initMatrix_ForCoeffs(double *M, int n1,
                          void (*get_coeffs)(double *,double *, int))
{
  int i,j;
  double *u;
  double *c;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));

  if( !(u && c) ) errorexit("initMatrix_ForCoeffs: out of memory for u, c");


  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    u[j]=1.0;

    get_coeffs(c, u, n1-1);
    
    /* set M */
    for(i=0; i<n1; i++) M[n1*i + j] = c[i];
    
    u[j]=0.0;
  }
  free(u);
  free(c);
}

/* init a n1*n1 matrix M used to compute function u from coeffs */
void initMatrix_ToEvaluate(double *M, int n1,
                           void (*eval_onPoints)(double *,double *, int))
{
  int i,j;
  double *u;
  double *c;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));

  if( !(u && c) ) errorexit("initMatrix_ForCoeffs: out of memory for u, c");


  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    c[j]=1.0;

    eval_onPoints(c, u, n1-1);
    
    /* set M */
    for(i=0; i<n1; i++) M[n1*i + j] = u[i];
    
    c[j]=0.0;
  }
  free(u);
  free(c);
}



/* get the coeffs c of 3d var u in direction direc on the whole box */
void spec_analysis1(tBox *box, int direc, double *M, double *u, double *c)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *cline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    cline = (double*) realloc(cline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n1);
        put_memline(c, cline, 1, j, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n2);
        put_memline(c, cline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n3);
        put_memline(c, cline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_analysis1: possible values for direction direc are 1,2,3.");
}




/* compute the 3d var u from the coeffs c in direction direc  
   on the whole box                                         */
void spec_synthesis1(tBox *box, int direc, double *M, double *u, double *c)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *cline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    cline = (double*) realloc(cline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        get_memline(c, cline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n1);
        put_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(c, cline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n2);
        put_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(c, cline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n3);
        put_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_synthesis1: possible values for direction direc are 1,2,3.");
}



/* get all relevant function points in one box in direction direc */
void get_spec_functionpointers(tBox *box, int direc,
     void (**get_coeffs)(double *,double *, int),
     void (**coeffs_of_deriv)(double, double, double *,double *, int),
     void (**coeffs_of_2ndderiv)(double, double, double *,double *, int),
     void (**eval_onPoints)(double *,double *, int),
     void (**filter_coeffs)(double *, int, int) )
{       
  char str[1000];

  *get_coeffs=NULL;
  *coeffs_of_deriv=NULL;
  *coeffs_of_2ndderiv=NULL;
  *eval_onPoints=NULL;
  *filter_coeffs=NULL;

  snprintf(str, 999, "box%d_basis%d", box->b, direc);
  //printf("%s=%s\n", str, Gets(str));
  if( Getv(str, "ChebExtrema") )
  {
    *get_coeffs = cheb_coeffs_fromExtrema;
    *coeffs_of_deriv = cheb_deriv;
    *eval_onPoints = cheb_eval_onExtrema;
    *filter_coeffs = cheb_filter;
  }
  else if( Getv(str, "Fourier") )
  {
    *get_coeffs = four_coeffs;
    *coeffs_of_deriv = four_deriv;
    *eval_onPoints = four_eval;
    *filter_coeffs = four_filter;
  }
  else if( Getv(str, "fd2_onesided") )
  {
    *get_coeffs = fd2_coeffs;
    *coeffs_of_deriv = fd2_deriv_onesidedBC;
    *eval_onPoints = fd2_eval;
    *filter_coeffs = fd2_filter;
    *coeffs_of_2ndderiv = fd2_2ndderiv_onesidedBC;
  }
  else if( Getv(str, "fd2_periodic") )
  {
    *get_coeffs = fd2_coeffs;
    *coeffs_of_deriv = fd2_deriv_periodic;
    *eval_onPoints = fd2_eval;
    *filter_coeffs = fd2_filter;
    *coeffs_of_2ndderiv = fd2_2ndderiv_periodic;
  }
}
