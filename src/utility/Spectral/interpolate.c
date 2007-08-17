/* interpolate.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Spectral.h"



/* get c=c_ijk of var u */
void spec_Coeffs(tBox *box, double *u, double *c)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double *M;
  double *d;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, double X)=NULL;
  int m3 = max3(n1,n2,n3);

  /* memory for matrix and temp var d */
  M = (double*) calloc( m3*m3, sizeof(double) );
  d = (double*) calloc( n1*n2*n3, sizeof(double) );
  if( !(M && d) ) errorexit("spec_Coeffs: out of memory for M, d");
         
  /* get c_ijk by calling spec_analysis1 in all 3 directions */
  /* direction 1 */
  get_spec_functionpointers(box, 1, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n1, get_coeffs);
  spec_analysis1(box, 1, M, u, c);

  /* direction 2 */
  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n2, get_coeffs);
  spec_analysis1(box, 2, M, c, d);

  /* direction 3 */
  get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n3, get_coeffs);
  spec_analysis1(box, 3, M, d, c);

  /* free memory for matrix M and temp var d */
  free(M);
  free(d);
}


/* use coeffs c to interpolate to X,Y,Z in box 
   (set coeffs c of u by calling spec_Coeffs(box, u, c); */
double spec_interpolate(tBox *box, double *c, double X, double Y, double Z)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;
  double sum=0.0;
  double a1=box->bbox[0];
  double b1=box->bbox[1];
  double a2=box->bbox[2];
  double b2=box->bbox[3];
  double a3=box->bbox[4];
  double b3=box->bbox[5];

  if(box->basis1==NULL || box->basis2==NULL || box->basis3==NULL)
    errorexiti("spec_interpolate: box%d: one box->basis1/2/3 is NULL",box->b);

  forallijk(i,j,k)
  {
    sum += c[Index(i,j,k)] * box->basis1(a1,b1, i, X) * 
           box->basis2(a2,b2, j, Y) * box->basis3(a3,b3, k, Z);
  }
  return sum;
}
