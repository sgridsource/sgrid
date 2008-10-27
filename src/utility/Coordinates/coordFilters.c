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
    else if(Getv(str, "Spherical")  ||
            Getv(str, "Spherical2") ||
            Getv(str, "Spherical3") )
      coordinateDependentFilter_Spherical(box, unew);
  }
}


/* Filter unew such that the the upper 1-sin(theta) portion of 
   the Fourier coeffs in the phi direction is zero.            */
void coordinateDependentFilter_SphericalDF(tBox *box, tVarList *unew)
{
  int vi;
  int i,j,k, n1,n2,n3;
  double *u;
  double *c = box->v[Ind("temp1")]; /* we store the coeffs in the variable ADMVars temp1 */
  double *thm = box->v[Ind("Y")];

  n1 = box->n1;
  n2 = box->n2;
  n3 = box->n3;

  /* for all variables */
  for(vi = 0; vi < unew->n; vi++)
  {
    u = box->v[unew->index[vi]];
    
    /* get spectral coeffs in c */
    spec_analysis1(box, 3, u, c);
      
    /* set the upper 1-sin(theta) portion of the coeffs c to zero */
    for(j = 0; j < n2; j++)
    {
      double theta = thm[Index(0,j,0)] + PI/((1+n2%2)*n2);
      int ms, ks;
      
      /* first m which we eliminate is m = ms = [N fabs(sin(theta)] + 1, 
         where N=n3/2                                                    */
      ms = n3/2;
      ms *= fabs(sin(theta));
      ms++;
      ks = 2*ms-1; // before we had: ks = n3*fabs(sin(theta)); if( ks%2 == 0 ) ks++;
      if( ks >= n3) continue;
      for(k = ks; k < n3; k++)
        for(i = 0; i < n1; i++)
          c[Index(i,j,k)] = 0.0;
    }

    /* get new u from new c */
    spec_synthesis1(box, 3, u, c);
  }
}


/* Filter unew such that the the upper 1-sin(theta) portion of 
   the Fourier coeffs in the phi direction is zero.            */
void coordinateDependentFilter_Spherical(tBox *box, tVarList *unew)
{
  int vi;
  int i,j,k, n1,n2,n3;
  double *u;
  double *c = box->v[Ind("temp1")]; /* we store the coeffs in the variable ADMVars temp1 */
  double *xp = box->v[Ind("x")];
  double *yp = box->v[Ind("y")];
  double *zp = box->v[Ind("z")];

  n1 = box->n1;
  n2 = box->n2;
  n3 = box->n3;

  /* for all variables */
  for(vi = 0; vi < unew->n; vi++)
  {
    u = box->v[unew->index[vi]];
    
    /* get spectral coeffs in c */
    spec_analysis1(box, 3, u, c);
      
    /* set the upper 1-sin(theta) portion of the coeffs c to zero */
    for(j = 0; j < n2; j++)
    {
      int I0j0 = Index(0,j,0);
      double x=xp[I0j0];
      double y=yp[I0j0];
      double z=zp[I0j0];
      double r = sqrt(x*x + y*y + z*z);
      double sintheta = sin(acos(z/r));
      int ms, ks;

      /* first m which we eliminate is m = ms = [N fabs(sin(theta)] + 1, 
         where N=n3/2                                                    */
      ms = n3/2;
      ms *= fabs(sintheta);
      ms++;
      ks = 2*ms-1; // before we had: ks = n3*fabs(sin(theta)); if( ks%2 == 0 ) ks++;
      if( ks >= n3) continue;
      for(k = ks; k < n3; k++)
        for(i = 0; i < n1; i++)
          c[Index(i,j,k)] = 0.0;
    }

    /* get new u from new c */
    spec_synthesis1(box, 3, u, c);
  }
}
