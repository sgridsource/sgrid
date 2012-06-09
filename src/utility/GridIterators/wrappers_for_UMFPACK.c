/* UMFPACK_solve_wrapper.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "GridIterators.h"


/* global vars */
/* UMFPACK matrix */
int *Ap;
int *Ai;
double *Ax;


/* switch to fin. diff. and use UMFPACK as lin solver as precon */
void precon_fd_UMFPACK(tVarList *vlx, tVarList *vlb,
                       tVarList *vlc1, tVarList *vlc2)
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;

  /* use umfpack for solve */
  INFO = umfpack_solve_from_Ap_Ai_Ax(Ap,Ai,Ax, vlx, vlb, pr);
  if(pr) 
    printf("precon_fd_UMFPACK: umfpack_solve_from_Ap_Ai_Ax returned INFO=%d\n",INFO);
}

/* call all relevant routines to set up an UMFPACK call */
int bicgstab_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  tGrid *grid_bak=make_empty_grid(grid->nvariables, 0);
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  double drop = Getd("GridIterators_setABStozero_below");
  int col, ncols, nvars, INFO;
  tSparseVector **Acol;
  int nz=0;

  if(pr) printf("bicgstab_with_fd_UMFPACK_precon: ");

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

  /* count number of entries in sparse matrix */
  for(col = 0; col < ncols; col++) nz+=Acol[col]->entries;

  /* allocate memory for matrix in UMFPACK format */
  allocate_umfpack_matrix(&Ap, &Ai, &Ax, ncols, nz);

  /* save current grid in grid_bak and then convert grid to fin. diff. */
  copy_grid_withoutvars(grid, grid_bak, 0);
  convert_grid_to_fd(grid);

  /* now set fin. diff. matrix for UMFPACK */
  set_umfpack_matrix_from_columns(Ap, Ai, Ax, Acol, ncols, drop, pr);
  if(pr)
    printf("bicgstab_with_fd_UMFPACK_precon:\n"
           "  %d entries of magnitude <= %g were dropped\n",
           nz-Ap[ncols], drop);

  /* restore grid to spectral */
  copy_grid_withoutvars(grid_bak, grid, 0);
  free_grid(grid_bak);

  /* solve A x = b with bicgstab and the Precon precon_fd_UMFPACK */
  INFO = bicgstab(x, b, r,c1,c2, itmax,tol,normres, lop, precon_fd_UMFPACK);

  /* free matrix A in UMFPACK format */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  free_umfpack_matrix(Ap, Ai, Ax);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}


/* call all relevant routines to set up an umfpack call */
int UMFPACK_solve_wrapper(tVarList *x, tVarList *b, 
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

  if(pr) printf("UMFPACK_solve_wrapper: ");

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
  INFO=umfpack_solve_fromAcolumns(Acol, x, b, drop, pr);
  if(pr) 
    printf("UMFPACK_solve_wrapper: umfpack_solve_fromAcolumns returned INFO=%d\n",INFO);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}


/* call all relevant routines to set up an umfpack call using
   umfpack_solve_forSortedVars_fromAcolumns */
int UMFPACK_solve_forSortedVars_wrapper(tVarList *x, tVarList *b, 
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

  if(pr) printf("UMFPACK_solve_forSortedVars_wrapper: ");

  /* allocate Acol */
  nvars=x->n; 
  ncols = 0;
  forallboxes(grid,bi)  ncols+=(grid->box[bi]->nnodes)*nvars;
  Acol = calloc(ncols, sizeof(*Acol));
  if(Acol)  if(pr) printf("allocated %d matrix columns\n", ncols);
  else       errorexit("no memory for Acol");
  for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();

  /* set Acol */                
  SetMatrixColumns_forSortedVars_slowly(Acol, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

  /* solve A x = b with umfpack */
  INFO=umfpack_solve_forSortedVars_fromAcolumns(Acol, x, b, drop, pr);
  if(pr) 
    printf("UMFPACK_solve_forSortedVars_wrapper: umfpack_solve_fromAcolumns returned INFO=%d\n",INFO);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}
