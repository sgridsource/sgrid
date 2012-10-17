/* UMFPACK_solve_wrapper.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "GridIterators.h"


/* global vars */
/* UMFPACK matrix */
int *umfAp;
int *umfAi;
LONGINT *LumfAp;
LONGINT *LumfAi;
double *umfAx;


/* use UMFPACK (with matrix set in global umfAp,umfAi,umfAx) as lin solver as precon */
void precon_Ap_Ai_Ax_UMFPACK(tVarList *vlx, tVarList *vlb,
                             tVarList *vlc1, tVarList *vlc2)
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;

  /* use umfpack for solve */
  if(Getv("GridIterators_UMFPACK_version", "di"))
    INFO = umfpack_di_solve_from_Ap_Ai_Ax(umfAp,umfAi,umfAx, vlx, vlb, pr);
  else
    INFO = umfpack_dl_solve_from_Ap_Ai_Ax(LumfAp,LumfAi,umfAx, vlx, vlb, pr);
    
  if(pr) 
    printf("precon_Ap_Ai_Ax_UMFPACK: umfpack_*_solve_from_Ap_Ai_Ax returned INFO=%d\n",INFO);
}

/* do a linear solve with precon_Ap_Ai_Ax_UMFPACK, where global
   umfAp,umfAi,umfAx is set by temporarily switching to finite differencing */
int linSolve_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
            int (*lsolver)(tVarList *, tVarList *, tVarList *, tVarList *,tVarList *,
	                   int imax, double tl, double *res,
	                   void (*Lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	                   void (*Prec)(tVarList *, tVarList *, tVarList *, tVarList *)),
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *)) 
{
  tGrid *grid = b->grid;
  tGrid *grid_bak=make_empty_grid(grid->nvariables, 0);
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  double drop = Getd("GridIterators_setABStozero_below");
  int col, ncols, nvars, INFO;
  tSparseVector **Acol;
  LONGINT nz=0;

  if(pr) printf("lsolver_with_fd_UMFPACK_precon:\n");

  /* allocate Acol */
  nvars=x->n; 
  ncols = 0;
  forallboxes(grid,bi)  ncols+=(grid->box[bi]->nnodes)*nvars;
  Acol = calloc(ncols, sizeof(*Acol));
  if(Acol) { if(pr) printf("allocated %d matrix columns\n", ncols); }
  else       errorexit("no memory for Acol");
  for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();

  /* save current grid in grid_bak and then convert grid to fin. diff. */
  copy_grid_withoutvars(grid, grid_bak, 0);
  convert_grid_to_fd(grid);

  /* set Acol */                
  SetMatrixColumns_slowly(Acol, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

  /* restore grid to spectral */
  copy_grid_withoutvars(grid_bak, grid, 0);
  free_grid(grid_bak);

  /* count number of entries in sparse matrix */
  for(col = 0; col < ncols; col++) nz+=Acol[col]->entries;

  if(Getv("GridIterators_UMFPACK_version", "di"))
  {
    /* allocate memory for matrix in UMFPACK format */
    allocate_umfpack_di_matrix(&umfAp, &umfAi, &umfAx, ncols, nz);
    /* now set fin. diff. matrix for UMFPACK */
    set_umfpack_di_matrix_from_columns(umfAp, umfAi, umfAx, Acol, ncols, drop, pr);
    if(pr)
      printf("lsolver_with_fd_UMFPACK_precon:\n"
             "  %d entries of magnitude <= %g were dropped\n",
             (int) (nz-umfAp[ncols]), drop);
  }
  else
  {
    /* allocate memory for matrix in UMFPACK format */
    allocate_umfpack_dl_matrix(&LumfAp, &LumfAi, &umfAx, ncols, nz);
    /* now set fin. diff. matrix for UMFPACK */
    set_umfpack_dl_matrix_from_columns(LumfAp, LumfAi, umfAx, Acol, ncols, drop, pr);
    if(pr)
      printf("lsolver_with_fd_UMFPACK_precon:\n"
             "  %ld entries of magnitude <= %g were dropped\n",
             (LONGINT) (nz-umfAp[ncols]), drop);
  }

  /* solve A x = b with lsolver and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = lsolver(x, b, r,c1,c2, itmax,tol,normres, lop, precon_Ap_Ai_Ax_UMFPACK);

  /* free matrix A in UMFPACK format */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  if(Getv("GridIterators_UMFPACK_version", "di"))
    free_umfpack_di_matrix(umfAp, umfAi, umfAx);
  else
    free_umfpack_dl_matrix(LumfAp, LumfAi, umfAx);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}

/* do a linear solve with bicgstab and precon_Ap_Ai_Ax_UMFPACK */
int bicgstab_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("bicgstab_with_fd_UMFPACK_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_fd_UMFPACK_precon(x, b, r,c1,c2, bicgstab,
                                         itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_gmres_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_gmres_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_gmres_wrapper_with_fd_UMFPACK_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_fd_UMFPACK_precon(x, b, r,c1,c2, templates_gmres_wrapper,
                                         itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicgstab_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_bicgstab_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicgstab_wrapper_with_fd_UMFPACK_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_fd_UMFPACK_precon(x, b, r,c1,c2, templates_bicgstab_wrapper,
                                         itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_cgs_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_cgs_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_cgs_wrapper_with_fd_UMFPACK_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_fd_UMFPACK_precon(x, b, r,c1,c2, templates_cgs_wrapper,
                                         itmax,tol,normres, lop);
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
  if(Getv("GridIterators_UMFPACK_version", "di"))
    INFO=umfpack_di_solve_fromAcolumns(Acol, x, b, drop, pr);
  else
    INFO=umfpack_dl_solve_fromAcolumns(Acol, x, b, drop, pr);
  if(pr) 
    printf("UMFPACK_solve_wrapper: "
           "umfpack_*_solve_fromAcolumns returned INFO=%d\n",INFO);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}


/* call all relevant routines to set up an umfpack call using
   umfpack_di_solve_forSortedVars_fromAcolumns */
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
  if(Getv("GridIterators_UMFPACK_version", "di"))
    INFO=umfpack_di_solve_forSortedVars_fromAcolumns(Acol, x, b, drop, pr);
  else /* INFO=umfpack_dl_solve_forSortedVars_fromAcolumns(Acol, x, b, drop, pr); */
    errorexit("implement umfpack_dl_solve_forSortedVars_fromAcolumns");
    
  if(pr) 
    printf("UMFPACK_solve_forSortedVars_wrapper: umfpack_di_solve_fromAcolumns returned INFO=%d\n",INFO);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}
