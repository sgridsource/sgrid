/* LinSolve_withLAPACK.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "Poisson.h"


int LinSolve_withLAPACK(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int bi;
  int pr = 1;
  int line, nlines, nvars, ok;
  tSparseVector **Aline;

  if(pr) printf("LinSolve_withLAPACK: ");

  /* allocate Aline */
  nvars=x->n; 
  nlines = 0;
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;
  Aline = calloc(nlines, sizeof(*Aline));
  if(Aline)  if(pr) printf("allocated %d matrix lines\n", nlines);
  else       errorexit("no memory for Aline");
  for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();

  /* set Aline */                
  SetMatrixLines_slowly(Aline, lop, r, x, c1, c2, pr);

  /* solve A x = b with lapack's dgesv */
  ok=lapack_dgesv(Aline, x, b, pr);
  if(pr) printf("LinSolve_withLAPACK: lapack_dgesv returned ok=%d\n", ok);

  return ok;
}
