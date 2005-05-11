/* boundary.c */
/* Wolfgang Tichy 5/2005 */

#include "sgrid.h"
#include "boundary.h"

/* here we use:
   tPointList *radiativeBoundaryPointList;
   tPointList *simpleExcisionBoundaryPointList; */
/* initialize our Boundary PointLists */
int initialize_BoundaryPointLists(tGrid *grid)
{
  char str[100];
  int bi;
  int i,j,k;
  int n1,n2,n3;
  int nboxes = grid->nboxes;
  
  radiativeBoundaryPointList = AllocatePointList(grid);
  simpleExcisionBoundaryPointList = AllocatePointList(grid);
  
  /* radiative boundary condition */
  if( Getv("boundary", "radiative") )
  {
    bi=nboxes-1;
    snprintf(str, 99, "box%d_Coordinates", bi);
    if(Getv(str, "SphericalDF"))
    {
       n1=grid->box[bi]->n1;
       n2=grid->box[bi]->n2;
       n3=grid->box[bi]->n3;
        
       forplane1(i,j,k, n1,n2,n3, n1-1)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));
    }
    else
    {
       n1=grid->box[bi]->n1;
       n2=grid->box[bi]->n2;
       n3=grid->box[bi]->n3;
        
       forplane1(i,j,k, n1,n2,n3, 0)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));

       forplane1(i,j,k, n1,n2,n3, n1-1)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));

       forplane2(i,j,k, n1,n2,n3, 0)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));

       forplane2(i,j,k, n1,n2,n3, n2-1)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));

       forplane3(i,j,k, n1,n2,n3, 0)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));

       forplane3(i,j,k, n1,n2,n3, n3-1)
         AddToPointList(radiativeBoundaryPointList, bi, Index(i,j,k));
    }
    printf("radiativeBoundaryPointList:\n");
    prPointList(radiativeBoundaryPointList);
  }
  return 0;
}


/* apply all active boundary conditions after a standard evolution step */
void set_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  tGrid *grid = unew->grid;
  int i,j;
  //double harmonic_f = 1;
  double v0, var0, varfalloff;
  
  /* radiative boundary condition */
  if( Getv("boundary", "radiative") )
  {
    /* monopole fall off rate: not implemented yet */
    varfalloff = 0;

    /* for all variables */
    for (j = 0; j < unew->n; j++)
    {
      i    = unew->index[j];
      v0   = VarPropSpeed(i);
      var0 = VarFarLimit(i);

      /* keep this variable constant, just copy it */
      if (Getv("boundary_radconstant", VarName(i))) v0 = 0;

      /* do it */
      set_boundary_radiative(grid, 
	unew->index[j], upre->index[j], c, ucur->index[j], var0, v0);
    }
  }

  /* symmetry boundary 
     has to be after radiative and excision
     can be before or after synchronization, but before is more efficient
  */
  //set_boundary_symmetry(level, unew); 
}

