/* wrappers_for_JacobiPrecon.c */
/* Wolfgang Tichy 1/2013 */


#include "sgrid.h"
#include "GridIterators.h"


/* global vars */
/* Diagonal Matrix elements*/
double *DiagM;


/* Diag contains the diagonal of a matrix set with SetMatrixColumns_slowly */
void Jacobi_Preconditioner_from_DiagM(tVarList *vlx, tVarList *vlb,
                                      tVarList *vlc1, tVarList *vlc2)
{
  tGrid *grid = vlx->grid;
  int bi, line;
  int pr = Getv("GridIterators_verbose", "yes");

  /* set vlx */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *xp = box->v[vlx->index[j]];
        double *bp = box->v[vlb->index[j]];
        xp[i] = bp[i]/DiagM[line];
        line++;
      }
  }
  /*
  if(pr)
  { 
    printf("Jacobi_Preconditioner_from_DiagM: vector vlx=%p is now set!\n", vlx);
    fflush(stdout);
  }
  */
}


/* do a linear solve with Jacobi_Preconditioner_from_DiagM, where global
   DiagM is set by temporarily switching to finite differencing */
int linSolve_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
            int (*lsolver)(tVarList *, tVarList *, tVarList *, tVarList *,tVarList *,
	                   int imax, double tl, double *res,
	                   void (*Lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	                   void (*Prec)(tVarList *, tVarList *, tVarList *, tVarList *)),
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *)) 
{
  tGrid *grid = b->grid;
  tGrid *grid_bak;
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  int col, ncols, nvars, ent, INFO;
  tSparseVector **Acol;
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");

  if(pr) printf("linSolve_with_Jacobi_precon\n");

  /* allocate Acol */
  nvars=x->n; 
  ncols = 0;
  forallboxes(grid,bi)  ncols+=(grid->box[bi]->nnodes)*nvars;
  Acol = calloc(ncols, sizeof(*Acol));
  if(Acol) { if(pr) printf("allocated %d matrix columns\n", ncols); }
  else       errorexit("no memory for Acol");
  for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();

  if(use_fd)
  {
    /* save current grid in grid_bak and then convert grid to fin. diff. */
    grid_bak = make_empty_grid(grid->nvariables, 0);
    copy_grid_withoutvars(grid, grid_bak, 0);
    convert_grid_to_fd(grid);
    if(pr) printf("Using finite differencing to set matrix Acol...\n");
  }

  /* set Acol */                
  SetMatrixColumns_slowly(Acol, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

  if(use_fd)
  {
    /* restore grid to spectral */
    copy_grid_withoutvars(grid_bak, grid, 0);
    free_grid(grid_bak);
  }

  /* allocate memory for diagonal of matrix in DiagM */
  DiagM = calloc(ncols, sizeof(*DiagM));

  /* now set DiagM to diagonal of matrix Acol */
  for(col=0; col<ncols; col++)
    for(ent = 0; ent < Acol[col]->entries; ent++)
      if(Acol[col]->pos[ent] == col) 
      {
        DiagM[col] = Acol[col]->val[ent];
        if(DiagM[col]==0.0) errorexit("DiagM is singular!!!");
        break;
      }
  
  /* solve A x = b with lsolver and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = lsolver(x, b, r,c1,c2, itmax,tol,normres, lop, Jacobi_Preconditioner_from_DiagM);

  /* free matrix DiagM */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  free(DiagM);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}

/* do a linear solve with bicgstab and precon_Ap_Ai_Ax_UMFPACK */
int bicgstab_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("bicgstab_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, bicgstab,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_gmres_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_gmres_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_gmres_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_gmres_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicgstab_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_bicgstab_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicgstab_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_bicgstab_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_cgs_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_cgs_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_cgs_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_cgs_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}
