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
     void (**filter_coeffs)(double *, int, int),
     double (**basisfunc)(double a, double b, int k, double X) )
{       
  char str[1000];

  *get_coeffs=NULL;
  *coeffs_of_deriv=NULL;
  *coeffs_of_2ndderiv=NULL;
  *eval_onPoints=NULL;
  *filter_coeffs=NULL;
  *basisfunc=NULL;

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
  else if( Getv(str, "ChebZeros") )
  {
    *get_coeffs = cheb_coeffs_fromZeros;
    *coeffs_of_deriv = cheb_deriv;
    *eval_onPoints = cheb_eval_onZeros;
    *filter_coeffs = cheb_filter;
  }
  else
  {
    printf("get_spec_functionpointers: %s=%s is unknown!\n", str, Gets(str));
    errorexits("get_spec_functionpointers: don't know what to do "
               "with %s" , Gets(str));
  }
}


/* compute B_k(x) M_ki     <-- B_k is basis function k
   Cu_k = M_ki u_i         <-- M_ki is coeff matrix
   u(x) = Cu_k B_k(x) = B_k(x) M_ki u_i = BM_i u_i      */
double spec_Basis_times_CoeffMatrix(tBox *box, int direc, 
                                    double *BM, double X)
{
  int linelen;
  double *M=NULL;
  double *B=NULL;
  int i;
  double a,b;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, double X)=NULL;
  /* basisfunc is something like cheb_basisfunc */

  get_spec_functionpointers(box, direc, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
 
  if(direc==1)      { linelen = box->n1; a=box->bbox[0]; b=box->bbox[1]; }
  else if(direc==2) { linelen = box->n2; a=box->bbox[2]; b=box->bbox[3]; }
  else if(direc==3) { linelen = box->n3; a=box->bbox[4]; b=box->bbox[5]; }
  else
    errorexit("spec_Basis_times_CoeffMatrix: possible values for direction direc are 1,2,3.");
    

  /* initialize the matrix M used to compute coeffs */
  M = (double *) calloc(linelen*linelen, sizeof(double));
  initMatrix_ForCoeffs(M, linelen, get_coeffs);

  /* initialize basis functions at point X */
  B = (double *) calloc(linelen, sizeof(double));
  for(i = 0; i < linelen; i++)  B[i] = basisfunc(a,b, i, X); // B[i]=B_i(X)

  /* Cu_k = M_ki u_i   <-- M is coeff matrix
     u(x) = Cu_k B_k(x) = B_k(x) M_ki u_i = BM_i u_i */

  /* get BM_i = B_k(x) M_ki */
  vector_times_matrix(B, M, BM, linelen);

  /* free memory for  matrix M and basis funcs B */
  free(M);
  free(B);
}
