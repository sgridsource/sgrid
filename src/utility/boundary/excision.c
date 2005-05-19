/* excision.c */
/* Wolfgang Tichy 5/2005 */
/* simple excsion boundary condition */

#include "sgrid.h"
#include "boundary.h"


/* set simple excision boundary for one variable */
void set_boundary_simpleExcision(tGrid *grid, int unew, int upre)
{
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *varnew = box->v[unew];
    double *varpre = box->v[upre];

    /* NOTE: This works only if excsion boundary is at i=0: */
    forPointList_inbox(ExcisionBoundaryPointList, box, pi , ijk)
      varnew[ijk] =  varpre[ijk] + varnew[ijk+1] - varpre[ijk+1];
  }
}


/* set VonNeumann excision boundary for one variable */
void set_boundary_VonNeumannExcision(tPointList *PL, int unew)
{
  /* NOTE: This works only if excsion boundary is at i=0: */
  set_boundary_normalderiv_leftBound(PL, 1, unew, 0.0);
}
