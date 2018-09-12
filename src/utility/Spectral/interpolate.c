/* interpolate.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Spectral.h"



/* get c=c_ijk of var u */
void spec_Coeffs(tBox *box, double *u, double *c)
{
  /* get c_ijk by calling spec_analysis1 in all 3 directions */
  spec_analysis1(box, 1, u, c);
  spec_analysis1(box, 2, c, c);
  spec_analysis1(box, 3, c, c);
}


/* get var u from coeffs c=c_ijk */
void spec_Eval(tBox *box, double *u, double *c)
{
  /* get u by calling spec_synthesis1 in all 3 directions */
  spec_synthesis1(box, 1, u, c);
  spec_synthesis1(box, 2, u, u);
  spec_synthesis1(box, 3, u, u);
}


/* old version of spec_Coeffs */
void spec_Coeffs_old(tBox *box, double *u, double *c)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double *M;
  double *d;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int n1, double X)=NULL;
  int m3 = max3(n1,n2,n3);

  /* memory for matrix and temp var d */
  M = (double*) calloc( m3*m3, sizeof(double) );
  d = (double*) calloc( n1*n2*n3, sizeof(double) );
  if( !(M && d) ) errorexit("spec_Coeffs: out of memory for M, d");
         
  /* get c_ijk by calling spec_analysis1 in all 3 directions */
  /* direction 1 */
  get_spec_functionpointers(box, 1, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n1, get_coeffs);
  spec_analysis1_diffmatrix(box, 1, M, u, c);

  /* direction 2 */
  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n2, get_coeffs);
  spec_analysis1_diffmatrix(box, 2, M, c, d);

  /* direction 3 */
  get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  initMatrix_ForCoeffs(M, n3, get_coeffs);
  spec_analysis1_diffmatrix(box, 3, M, d, c);

  /* free memory for matrix M and temp var d */
  free(M);
  free(d);
}


/* use coeffs c to interpolate to X,Y,Z in box 
   (set coeffs c of u by calling spec_Coeffs(box, u, c); */
double spec_interpolate(tBox *box, double *c, double X, double Y, double Z)
{
  double *B1;
  double *B2;
  double *B3;
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

  /* allocate memory for B1/2/3 */
  B1 = (double*) malloc(n1 * sizeof(double));
  B2 = (double*) malloc(n2 * sizeof(double));
  B3 = (double*) malloc(n3 * sizeof(double));

  /* set basis func values at X,Y,Z */
  for(i = n1-1; i >=0; i--)  B1[i]=box->basis1((void *) box, a1,b1, i,n1, X);
  for(j = n2-1; j >=0; j--)  B2[j]=box->basis2((void *) box, a2,b2, j,n2, Y);
  for(k = n3-1; k >=0; k--)  B3[k]=box->basis3((void *) box, a3,b3, k,n3, Z);

  /* interpolate to X,Y,Z */
  /*  for(k = n3-1; k >=0; k--)
      for(j = n2-1; j >=0; j--)
      for(i = n1-1; i >=0; i--)
        sum += c[Index(i,j,k)] * B1[i] * B2[j] * B3[k]; */
  SGRID_LEVEL3_Pragma(omp parallel for reduction(+:sum)) 
  for(k = n3-1; k >=0; k--)
  { int j,i;
    for(j = n2-1; j >=0; j--)
    for(i = n1-1; i >=0; i--)
      sum += c[Index(i,j,k)] * B1[i] * B2[j] * B3[k];
  }

  /* free memory for B1/2/3 */
  free(B1);
  free(B2);
  free(B3);

  return sum;
}


/* get coeffs c of var u in plane N with index p */
void spec_Coeffs_inplaneN(tBox *box, int N, int p, double *u, double *c)
{
  /* get c calling spec_analysis1_inplaneN in both directions */
  if(N==3)
  {
    spec_analysis1_inplaneN(box, 1, N,p, u, c);
    spec_analysis1_inplaneN(box, 2, N,p, c, c);
  }
  else if(N==2)
  {
    spec_analysis1_inplaneN(box, 1, N,p, u, c);
    spec_analysis1_inplaneN(box, 3, N,p, c, c);
  }
  else if(N==1)
  {
    spec_analysis1_inplaneN(box, 2, N,p, u, c);
    spec_analysis1_inplaneN(box, 3, N,p, c, c);
  }
  else
    errorexit("spec_Coeffs_inplaneN: N must be 1, 2, or 3.");
}

/* get var u from coeffs c in plane N with index p */
void spec_Eval_inplaneN(tBox *box, int N, int p, double *u, double *c)
{
  /* get c calling spec_synthesis1_inplaneN in both directions */
  if(N==3)
  {
    spec_synthesis1_inplaneN(box, 1, N,p, u, c);
    spec_synthesis1_inplaneN(box, 2, N,p, u, u);
  }
  else if(N==2)
  {
    spec_synthesis1_inplaneN(box, 1, N,p, u, c);
    spec_synthesis1_inplaneN(box, 3, N,p, u, u);
  }
  else if(N==1)
  {
    spec_synthesis1_inplaneN(box, 2, N,p, u, c);
    spec_synthesis1_inplaneN(box, 3, N,p, u, u);
  }
  else
    errorexit("spec_Eval_inplaneN: N must be 1, 2, or 3.");
}

/* use coeffs c in plane N with index p to interpolate to (X1,X2) = (X,Y), (X,Z) or
   (Y,Z) in this plane. This is 2d interpolation.
   [set coeffs c of u by calling spec_Coeffs_inplaneN(box, N,p, u, c); ]   */
double spec_interpolate_inplaneN(tBox *box, int N, int p, double *c,
                                 double X1, double X2)
{
  double *B1;
  double *B2;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int nmax = max3(n1, n2, n3);
  int i,j,k;
  double sum=0.0;
  double a1=box->bbox[0];
  double b1=box->bbox[1];
  double a2=box->bbox[2];
  double b2=box->bbox[3];
  double a3=box->bbox[4];
  double b3=box->bbox[5];

  if(box->basis1==NULL || box->basis2==NULL || box->basis3==NULL)
    errorexiti("spec_interpolate_inplaneN: box%d: one box->basis1/2/3 is NULL",box->b);

  /* allocate memory for B1/2 */
  B1 = (double*) malloc(nmax * sizeof(double));
  B2 = (double*) malloc(nmax * sizeof(double));

  /* set basis func values at X1,X2 */
  if(N==3)
  {
    for(i = n1-1; i >=0; i--)  B1[i]=box->basis1((void *) box, a1,b1, i,n1, X1);
    for(j = n2-1; j >=0; j--)  B2[j]=box->basis2((void *) box, a2,b2, j,n2, X2);

    /* interpolate to X,Y */
    forplane3(i,j,k, n1,n2,n3, p)
      sum += c[Index(i,j,k)] * B1[i] * B2[j];
  }
  else if(N==2)
  {
    for(i = n1-1; i >=0; i--)  B1[i]=box->basis1((void *) box, a1,b1, i,n1, X1);
    for(k = n3-1; k >=0; k--)  B2[k]=box->basis3((void *) box, a3,b3, k,n3, X2);

    /* interpolate to X,Z */
    forplane2(i,j,k, n1,n2,n3, p)
      sum += c[Index(i,j,k)] * B1[i] * B2[k];
  }
  else if(N==1)
  {
    for(j = n2-1; j >=0; j--)  B1[j]=box->basis2((void *) box, a2,b2, j,n2, X1);
    for(k = n3-1; k >=0; k--)  B2[k]=box->basis3((void *) box, a3,b3, k,n3, X2);

    /* interpolate to Y,Z */
    forplane1(i,j,k, n1,n2,n3, p)
      sum += c[Index(i,j,k)] * B1[j] * B2[k];
  }
  else
    errorexit("spec_interpolate_inplaneN: N must be 1, 2, or 3.");

  /* free memory for B1/2 */
  free(B2);
  free(B1);

  return sum;
}

/* use coeffs c along a line in direction dir at index
   (i1,i2) = (j,k), (i,k), or (i,j)
   to interpolate to X3 = X, Y, or Z along this same line.
   This is 1d interpolation.
   [set coeffs c of u by calling spec_analysis1 or spec_analysis1_inplaneN ] */
double spec_interpolate_in_dir_at_i1_i2(tBox *box, int dir, int i1, int i2,
                                        double *c, double X3)
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
    errorexiti("spec_interpolate_in_dir_at_i1_i2: box%d: one box->basis1/2/3 is NULL",box->b);

  if(dir==1)
  {
    /* interpolate to X */
    for(i = n1-1; i >=0; i--)
      sum += c[Index(i,i1,i2)] * box->basis1((void *) box, a1,b1, i,n1, X3);
  }
  else if(dir==2)
  {
    /* interpolate to Y */
    for(j = n2-1; j >=0; j--)
      sum += c[Index(i1,j,i2)] * box->basis2((void *) box, a2,b2, j,n2, X3);
  }
  else if(dir==3)
  {
    /* interpolate to Z */
    for(k = n3-1; k >=0; k--)
      sum += c[Index(i1,i2,k)] * box->basis3((void *) box, a3,b3, k,n3, X3);
  }
  else
    errorexit("spec_interpolate_in_dir_at_i1_i2: dir must be 1, 2, or 3.");

  return sum;
}


/* get vlc=c_ijk of varlist vlu in box */
void spec_Coeffs_varlist(tBox *box, tVarList *vlu, tVarList *vlc)
{
  int j;

  for(j=0; j<vlu->n; j++)
  {
    double *u = box->v[vlu->index[j]]; 
    double *c = box->v[vlc->index[j]]; 
    spec_Coeffs(box, u, c);
  }
}

/* get varlist vlu from coeff varlist vlc=c_ijk */
void spec_Eval_varlist(tBox *box, tVarList *vlu, tVarList *vlc)
{
  int j;

  for(j=0; j<vlc->n; j++)
  {
    double *u = box->v[vlu->index[j]]; 
    double *c = box->v[vlc->index[j]]; 
    spec_Eval(box, u, c);
  }
}

/* interpolate one Var with index vind from grid2 to grid1, use var with 
   index tempind on grid2 as temporary storage for coefficients */
void spec_interpolate_Var_from_grid2_to_grid1(tGrid *grid1, tGrid *grid2,
                                              int vind, int tempind)
{
  int cind = tempind;
  int Xind = Ind("X");
  int Yind = Ind("Y");
  int Zind = Ind("Z");
  int b, i;

  /* save coeffs of vind on grid2 in cind */
  forallboxes(grid2, b)
  {
    tBox *box = grid2->box[b];
    if(box->v[vind]!=NULL)
      spec_Coeffs(box, box->v[vind], box->v[cind]);
  }

  /* loop over grid1 and interpolate vind from grid2 */
  forallboxes(grid1, b)
  {
    tBox *box = grid1->box[b];
    int b2 = b; /* use same box on grid2 as on grid1 */
    double *pX = box->v[Xind];
    double *pY = box->v[Yind];
    double *pZ = box->v[Zind];
    double *pv = box->v[vind];

    /* loop over points on grid1 */
    if(pv!=NULL) forallpoints(box, i)
    {
      double X = pX[i];
      double Y = pY[i];
      double Z = pZ[i];

      /* get var vind at point X,Y,Z by interpolation */
      pv[i] = spec_interpolate(grid2->box[b2], grid2->box[b2]->v[cind], X,Y,Z);
    }
  }
}
