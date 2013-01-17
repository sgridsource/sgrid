/* wrappers_for_SuiteSparseQR.c */
/* Wolfgang Tichy 1/2013 */


#include "sgrid.h"
#include "GridIterators.h"



/* call all relevant routines to set up a SuiteSparseQR call */
int SuiteSparseQR_solve_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  double drop = Getd("GridIterators_setABStozero_below");
  int col, ncols, nvars, INFO;
  tSparseVector **Acol;

  if(pr) printf("SuiteSparseQR_solve_wrapper: ");

  /* allocate Acol */
  nvars=x->n; 
  ncols = 0;
  forallboxes(grid,bi)  ncols+=(grid->box[bi]->nnodes)*nvars;
  Acol = calloc(ncols, sizeof(*Acol));
  if(Acol) { if(pr) printf("allocated %d matrix columns\n", ncols); }
  else       errorexit("no memory for Acol");
  for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();

  /* set Acol */                
  SetMatrixColumns_slowly(Acol, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

  /* solve A x = b with umfpack */
  INFO=SuiteSparseQR_solve_fromAcolumns(Acol, x, b, drop, pr);
  if(pr) 
    printf("SuiteSparseQR_solve_wrapper: "
           "SuiteSparseQR_solve_fromAcolumns returned INFO=%d\n",INFO);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}
