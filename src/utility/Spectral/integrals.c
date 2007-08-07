/* integrals.c */
/* Wolfgang Tichy 8/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* compute U = integral of 3d var u over full interval 
   in direction direc on a box                         */
void spec_Integral1(tBox *box, int direc, double *u, double *U)
{
  int linelen;
  double *M=NULL;
  int i,j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, double X)=NULL;

  get_spec_functionpointers(box, direc, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  if(direc==1)
  {
    linelen = n1;

    /* initialize the matrix M used to compute coeffs */
    M = (double *) calloc( linelen*linelen, sizeof(double));
    initMatrix_ForCoeffs(M, linelen, get_coeffs);

    /* write spec coeffs into U */
    spec_analysis1(box, direc, M, u, U);

    /* write cheb-integral from a to b into U */
    if(get_coeffs==cheb_coeffs_fromExtrema)
      for(k = 0; k < n3; k++)
        for(j = 0; j < n2; j++)
        {
          double sum;
          double L = box->bbox[1] - box->bbox[0];

          sum = 0.5*U[Index(0,j,k)]*L;
          for(i = 2; i < n1; i+=2)
            sum += U[Index(i,j,k)]* L/(1.0-i*i);

          for(i = 0; i < n1; i++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(get_coeffs==four_coeffs)
      for(k = 0; k < n3; k++)
        for(j = 0; j < n2; j++)
        {
          double U0LoN = U[Index(0,j,k)]*(box->bbox[1]-box->bbox[0])/n1;

          for(i = 0; i < n1; i++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else if(direc==2)
  {
    linelen = n2;

    /* initialize the matrix M used to compute coeffs */
    M = (double *) calloc( linelen*linelen, sizeof(double));
    initMatrix_ForCoeffs(M, linelen, get_coeffs);

    /* write spec coeffs into U */
    spec_analysis1(box, direc, M, u, U);

    /* write cheb-integral from a to b into U */
    if(get_coeffs==cheb_coeffs_fromExtrema)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[3] - box->bbox[2];

          sum = 0.5*U[Index(i,0,k)]*L;
          for(j = 2; j < n2; j+=2)
            sum += U[Index(i,j,k)]* L/(1.0-j*j);

          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(get_coeffs==four_coeffs)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double U0LoN = U[Index(i,0,k)]*(box->bbox[3]-box->bbox[2])/n2;

          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else if(direc==3)
  {
    linelen = n3;

    /* initialize the matrix M used to compute coeffs */
    M = (double *) calloc( linelen*linelen, sizeof(double));
    initMatrix_ForCoeffs(M, linelen, get_coeffs);

    /* write spec coeffs into U */
    spec_analysis1(box, direc, M, u, U);

    /* write cheb-integral from a to b into U */
    if(get_coeffs==cheb_coeffs_fromExtrema)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[5] - box->bbox[4];

          sum = 0.5*U[Index(i,j,0)]*L;
          for(k = 2; k < n3; k+=2)
            sum += U[Index(i,j,k)]* L/(1.0-k*k);

          for(k = 0; k < n3; k++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(get_coeffs==four_coeffs)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++)
        {
          double U0LoN = U[Index(i,j,0)]*(box->bbox[5]-box->bbox[4])/n3;

          for(k = 0; k < n3; k++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else
    errorexit("spec_Integral1: possible values for direction direc are 1,2,3.");

  /* free memory for  matrix M */
  free(M);
}


/* compute surface integral over surfaces normal to norm */
void spec_2dIntegral(tBox *box, int norm, double *u, double *U)
{
  if(norm==1)
  {
    spec_Integral1(box, 2, u, U);
    spec_Integral1(box, 3, U, U);
  }
  else if(norm==2)
  {
    spec_Integral1(box, 1, u, U);
    spec_Integral1(box, 3, U, U);
  }
  else if(norm==3)
  {
    spec_Integral1(box, 1, u, U);
    spec_Integral1(box, 2, U, U);
  }
  else
    errorexit("spec_2dIntegral: possible values for normal norm are 1,2,3.");
}


/* compute volume integral over 3d var u */
double spec_3dIntegral(tBox *box, double *u, double *U)
{
  spec_Integral1(box, 1, u, U);
  spec_Integral1(box, 2, U, U);
  spec_Integral1(box, 3, U, U);
  
  return U[0];
}


/* compute U = 2d integral of var u over theta and phi */
void spec_sphericalDF2dIntegral(tBox *box, double *u, double *U)
{
  int linelen;
  double *M=NULL;
  int i,j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, double X)=NULL;
  double *pX = box->v[Ind("X")];

  /* do phi integral */
  spec_Integral1(box, 3, u, U);

  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  {
    linelen = n2;

    /* initialize the matrix M used to compute coeffs */
    M = (double *) calloc( linelen*linelen, sizeof(double));
    initMatrix_ForCoeffs(M, linelen, get_coeffs);

    /* write spec coeffs into U */
    spec_analysis1(box, 2, M, U, U);

    /* write four-integral from a to b into U */
    if(get_coeffs==four_coeffs)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[3] - box->bbox[2];
          int n;
          int N = box->n2;
          /* double theta = thm + PI/((1+N%2)*N); */
          double d = 1.0/(2.0*(1+N%2)*N);
          double PI2 = 2.0*PI;
          double Re_c_n, Im_c_n;

          /* sum over all theta-integrated terms */
          sum = (1.0/PI) * U[Index(i,0,k)];
          sum += 0.5*( sin(PI2*d) * U[Index(i,1,k)]
                      +cos(PI2*d) * U[Index(i,2,k)] );
          for(n=2;n<N/2;n+=2)
          {
            Re_c_n = U[Index(i, 2*n-1, k)]; /* c[2*n-1]; */
            Im_c_n = U[Index(i, 2*n, k)];   /* c[2*n];   */
            sum += 2.0*( cos(PI2*n*d)/((1-n*n)*PI) * Re_c_n 
                        +sin(PI2*n*d)/((n*n-1)*PI) * Im_c_n );
          }
          if( N%4 == 0 )
            sum += cos(PI2*(N/2)*d)/((1-N*N/4)*PI) * U[Index(i,N-1,k)];

          /* adjust sum for L and N to obtain integral over theta */
          sum *= L/N;

          /* write integral into U along direction 2, and multiply by r^2  */
          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = sum*pX[Index(i,j,k)]*pX[Index(i,j,k)];
        }
    else errorexit("spec_sphericalDF2dIntegral: you need Fourier in direction 2");
  }

  /* free memory for  matrix M */
  free(M);
}


/* compute volume integral over var u if we use sphericalDF */
double spec_sphericalDF3dIntegral(tBox *box, double *u, double *U)
{
  spec_sphericalDF2dIntegral(box, u, U);
  spec_Integral1(box, 1, U, U);
  
  return U[0];
}
