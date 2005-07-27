/* doubleCovering.c */
/* Wolfgang Tichy 7/2005 */

/* deal with points which are covered twice */

#include "sgrid.h"
#include "Coordinates.h"


/* reset double covered points. So far onlt done for SphericalDF */
void reset_doubleCoveredPoints(tVarList *unew)
{
  tGrid *grid = unew->grid;
  int bi;

  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", bi);
    if( Getv(str, "SphericalDF") )
      reset_doubleCoveredPoints_SphericalDF(box, unew);
  }
}

/* Set the double covered points in SphericalDF to value found in normal
   region of 0<theta<PI, 0<phi<2PI:
   theta = 2PI j/n2 + PI/n2
    u(theta,phi) = u(2PI-theta,phi+PI) = u(2PI-theta,phi-PI)
=>  u[i,j,k]     = u[i,n2-j-1,k+n3/2]  = u[i,n2-j-1,k-n3/2] 
if both n2 and n3 are even!!!
*/
void reset_doubleCoveredPoints_SphericalDF(tBox *box, tVarList *unew)
{
  int vi;
  int i,j,k, n1,n2,n3;
  double *u;
  double av;

  /* for all variables */
  for(vi = 0; vi < unew->n; vi++)
  {
    u = box->v[unew->index[vi]];
    n1 = box->n1;
    n2 = box->n2;
    n3 = box->n3;

    if( n2%2 || n3%2 ) 
      errorexit("reset_doubleCoveredPoints_SphericalDF: "
                "n2 and n3 must be even!");

    /* copy points into double covered regions */
    /*
    for(k = 0; k < n3/2; k++)
      for(j = n2/2; j < n2; j++)
        for(i = 0; i < n1; i++)
          u[Index(i,j,k)] = u[Index(i,n2-j-1,k+n3/2)];

    for(k = n3/2; k < n3; k++)
      for(j = n2/2; j < n2; j++)
        for(i = 0; i < n1; i++)
          u[Index(i,j,k)] = u[Index(i,n2-j-1,k-n3/2)];
    */
    /* use average of u at the double covered points */
    for(k = 0; k < n3/2; k++)
      for(j = n2/2; j < n2; j++)
        for(i = 0; i < n1; i++)
        {
          av = ( u[Index(i,j,k)] + u[Index(i,n2-j-1,k+n3/2)] )*0.5;
          u[Index(i,j,k)] = av;
          u[Index(i,n2-j-1,k+n3/2)] = av;;
        }

    for(k = n3/2; k < n3; k++)
      for(j = n2/2; j < n2; j++)
        for(i = 0; i < n1; i++)
        {
          av = ( u[Index(i,j,k)] + u[Index(i,n2-j-1,k-n3/2)] )*0.5;
          u[Index(i,j,k)] = av;
          u[Index(i,n2-j-1,k-n3/2)] = av;
        }
  }
}
