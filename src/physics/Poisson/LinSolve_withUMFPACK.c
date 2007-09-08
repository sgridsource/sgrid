/* LinSolve_withUMFPACK.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "Poisson.h"


int LinSolve_withUMFPACK(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int bi;
  int pr = 1;
  int line, nlines, nvars, INFO;
  tSparseVector **Aline;

  if(pr) printf("LinSolve_withUMFPACK: ");

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
  if(pr&&0) 
    for(line=0; line<nlines; line++) prSparseVector(Aline[line]);

  /* solve A x = b with umfpack */
  INFO=umfpack_solve(Aline, x, b, pr);
  if(pr) printf("LinSolve_withUMFPACK: umfpack_solve returned INFO=%d\n",INFO);

  /* free matrix Aline */
  for(line=0; line<nlines; line++)  FreeSparseVector(Aline[line]);
  free(Aline);

  return INFO;
}
