/* filters.c */
/* Wolfgang Tichy 3/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* init a n1*n1 filter. matrix */
void initfiltermatrix(double *F, int k, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*filter_coeffs)(double *, int, int),
                    void (*eval_onPoints)(double *,double *, int) )
{
  int i,j;
  double *u;
  double *c;
  double *uf;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));
  uf = (double *) calloc(n1, sizeof(double));

  if( !(u && c && uf) ) errorexit("initfiltermatrix: out of memory for u, c, uf");


  /* read matrix from functions */
  if(k<=n1) 
  for(j=0; j<n1; j++)
  {
    u[j]=1.0;
    //cheb_coeffs_fromExtrema(c, u, n1-1);
    get_coeffs(c, u, n1-1);

    /* filter: zero all c[j] with k<=j<=n1-1 */
    //cheb_filter(c, k-1, n1-1);
    filter_coeffs(c, k-1, n1-1);    
    
    //cheb_eval_onExtrema(uf, c, n1-1);
    eval_onPoints(c, uf, n1-1);
    
    /* set F */
    for(i=0; i<n1; i++) F[n1*i + j] = uf[i];
    
    u[j]=0.0;
  }
  else
    for(i=0; i<n1; i++) F[n1*i + i] = 1.0;

  free(u);
  free(c);
  free(uf);
}


/* this is just a matrix multiplication */
/*
void filtermatrix(double *F, double *u, double *uf, int n)
{
  int i,j;
  double sum;
  
  for(i=0; i<n; i++)
  {
    sum=0.0;
    for(j=0; j<n; j++)  sum += F[n*i + j] * u[j];
    uf[i] = sum;
  }
}
*/


/* Filter with filter matrices */

/* filter 3d var u in dirction direc on a box */
void spec_filter1(tBox *box, int direc, double *u)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *ufline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    ufline = (double*) realloc(ufline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        get_memline(u, uline,  1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->F1, uline, ufline, box->n1);
        put_memline(u, ufline, 1, j, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline,  2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->F2, uline, ufline, box->n2);
        put_memline(u, ufline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline,  3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->F3 , uline, ufline, box->n3);
        put_memline(u, ufline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("cheb_filter1: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(ufline);
    uline = ufline = NULL;
  */
}

