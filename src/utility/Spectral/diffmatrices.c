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


/* reset all matrices and basis funcs on grid to do finite differencing */
void convert_grid_to_fd_onesidedBC(tGrid *grid)
{
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int vind;
    int i;
    double *temp;

    /* direction 1 */
    temp=box->Mcoeffs1; /* use box->Mcoeffs1 as temp storage, will be overwritten soon after */
    vind=Ind("X");
    for(i=0; i<n1; i++)  temp[i] = box->v[vind][Index(i,0,0)];
    init_fdcentered_diffmatrix(temp, box->D1, n1, fdcentered_deriv_onesidedBC);
    box->basis1=fd_basis1;
    initMatrix_ForCoeffs(box->Mcoeffs1, n1, fd2_coeffs);
    initMatrix_ToEvaluate(box->Meval1,  n1, fd2_eval);

    /* direction 2 */
    temp=box->Mcoeffs2; /* use box->Mcoeffs2 as temp storage, will be overwritten soon after */
    vind=Ind("Y");
    for(i=0; i<n2; i++)  temp[i] = box->v[vind][Index(0,i,0)];
    init_fdcentered_diffmatrix(temp, box->D2, n2, fdcentered_deriv_onesidedBC);
    box->basis2=fd_basis2;
    initMatrix_ForCoeffs(box->Mcoeffs2, n2, fd2_coeffs);
    initMatrix_ToEvaluate(box->Meval2,  n2, fd2_eval);

    /* direction 3 */
    temp=box->Mcoeffs3; /* use box->Mcoeffs3 as temp storage, will be overwritten soon after */
    vind=Ind("Z");
    for(i=0; i<n3; i++)  temp[i] = box->v[vind][Index(0,0,i)];
    init_fdcentered_diffmatrix(temp, box->D3, n3, fdcentered_deriv_onesidedBC);
    box->basis3=fd_basis3;
    initMatrix_ForCoeffs(box->Mcoeffs3, n3, fd2_coeffs);
    initMatrix_ToEvaluate(box->Meval3,  n3, fd2_eval);
  }      
}
