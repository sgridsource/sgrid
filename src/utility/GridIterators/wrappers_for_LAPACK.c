/* wrappers_for_LAPACK.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "GridIterators.h"


int LAPACK_dgesv_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  int line, nlines, nvars, INFO;
  tSparseVector **Aline;

  if(pr) printf("LAPACK_dgesv_wrapper: ");

  /* allocate Aline */
  nvars=x->n; 
  nlines = 0;
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;
  Aline = calloc(nlines, sizeof(*Aline));
  if(Aline) { if(pr) printf("allocated %d matrix lines\n", nlines); }
  else        errorexit("no memory for Aline");
  for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();

  /* set Aline */                
  SetMatrixLines_slowly(Aline, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(line=0; line<nlines; line++) prSparseVector(Aline[line]);

  /* solve A x = b with lapack's dgesv */
  INFO=lapack_dgesv(Aline, x, b, pr);
  if(pr) printf("LAPACK_dgesv_wrapper: lapack_dgesv returned INFO=%d\n", INFO);

  /* free matrix Aline */
  for(line=0; line<nlines; line++)  FreeSparseVector(Aline[line]);
  free(Aline);

  return INFO;
}
