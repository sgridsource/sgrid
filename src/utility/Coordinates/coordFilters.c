/* coordFilters.c */
/* Wolfgang Tichy 7/2005 */
/* Filter in a way that depends on the coordinate location */

#include "sgrid.h"
#include "Coordinates.h"



/* filter a var list unew */
void coordinateDependentFilter(tVarList *unew)
{
  tGrid *grid = unew->grid;
  int bi;

  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", bi);
    if( Getv(str, "SphericalDF") )
      coordinateDependentFilter_SphericalDF(box, unew);
  }
}

/* Filter unew such that the the upper 1-sin(theta) portion of 
   the Fourier coeffs in the phi direction is zero.            */
void coordinateDependentFilter_SphericalDF(tBox *box, tVarList *unew)
{
  int vi;
  int i,j,k, n1,n2,n3;
  double *u;
  double *F;
  double *B;
  double *c = box->v[Ind("temp1")]; /* we store the coeffs in the variable ADMVars temp1 */
  double *thm = box->v[Ind("Y")];

  n1 = box->n1;
  n2 = box->n2;
  n3 = box->n3;

  /* initialize the matrix used to compute Fourier coeffs */
  F = (double *) calloc(n3, sizeof(double));
  initMatrix_ForCoeffs(F, n3, four_coeffs);

  /* initialize the matrix used to get var u from Fourier coeffs */
  B = (double *) calloc(n3, sizeof(double));
  initMatrix_ForCoeffs(B, n3, four_eval);

  /* for all variables */
  for(vi = 0; vi < unew->n; vi++)
  {
    u = box->v[unew->index[vi]];
    
    /* get spectral coeffs in c */
    spec_analysis1(box, 3, F, u, c);
    
    /* set the upper 1-sin(theta) portion of the coeffs to zero */
    /*
    for(k = 0; k < n3; k++)
      for(j = 0; j < n2; j++)
        for(i = 0; i < n1; i++)
        {
          double theta = thm[Index(i,j,k)] + PI/((1+n2%2)*n2);
          
          if( k > sin(theta)*n3 )  c[Index(i,j,k)] = 0.0;
        }
    */

    /* set the upper 1-sin(theta) portion of the coeffs c to zero */
    for(j = 0; j < n2; j++)
    {
      double theta = thm[Index(0,j,0)] + PI/((1+n2%2)*n2);
      int ks = n3*sin(theta);
      
      if( ks%2 == 0 ) ks++;
      for(k = ks; k < n3; k++)
        for(i = 0; i < n1; i++)
          c[Index(i,j,k)] = 0.0;
    }

    /* get new u from new c */
    spec_synthesis1(box, 3, B, u, c);
  }
  free(F);
  free(B);
}
