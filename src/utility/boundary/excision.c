/* excision.c */
/* Wolfgang Tichy 5/2005 */
/* simple excsion boundary condition */

#include "sgrid.h"
#include "boundary.h"


/* set excision boundary for one variable */
void set_boundary_simpleExcision(tGrid *grid, int unew, int upre)
{
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *varnew = box->v[unew];
    double *varpre = box->v[upre];

//if(v==0) printf("v=0 ");

    /* NOTE: This works only if excsion bound is at i=0: */
    forPointList_inbox(simpleExcisionBoundaryPointList, box, pi , ijk)
      varnew[ijk] =  varpre[ijk] + varnew[ijk+1] - varpre[ijk+1];
  }
}
