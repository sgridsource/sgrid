/* SOR.c */
/* Wolfgang Tichy 1/2013 */

#include "sgrid.h"
#include "GridIterators.h"

/* global vars */
/* Matrix for SOR */
tSparseVector **SOR_Aline;
int SOR_nlines;


/* basic SOR step */
void SOR_step(tSparseVector **Aline, int nlines,
              double *x, double *b, double omega)
{
  int i;

  /* loop over x[i] */
  for(i=0; i<nlines; i++)
  {
    int ent;
    double Aii, sig;

    /* compute Aij * x[j] */
    sig = Aii = 0.0;
    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       int j;
       double Aij;

       j   = Aline[i]->pos[ent];
       Aij = Aline[i]->val[ent];
       if(j!=i)  sig += Aij * x[j];
       else      Aii = Aij;
    }
    if(Aii==0.0) errorexit("SOR_step failed, because Aii=0.");
    
    /* set x[i] = (1-omega) x[i] + (omega/Aii)(b[i] - sig) */
    x[i] = (1.0 - omega)*x[i] + (omega/Aii)*(b[i] - sig);
  }
}

/* SOR for a matrix set with SetMatrixLines */
int matrix_SOR(tSparseVector **Aline, int nlines,
               double *x, double *b, double omega,
               int itmax, double tol, double *normres)
{
  int i, it;
  int pr = Getv("GridIterators_verbose", "yes");
  double *f;

  f  = calloc(nlines, sizeof(double));
  if(f==NULL) errorexit("matrix_SOR: no memory for f");

  /* compute the residual f_i */
  /* f_i = sum_j A_ij x_j - b_i */
  SparseMatrixLines_times_vector(Aline, nlines, x, f);
  for(i=0; i<nlines; i++)  f[i] -= b[i];
  *normres = sqrt(scalarproduct_vectors(f,f, nlines));
  if(pr) 
  { 
    printf("matrix_SOR: nlines=%d  omega=%g  itmax=%d  tol=%g\n", 
           nlines, omega, itmax, tol);
    printf("matrix_SOR: step %d  residual=%g\n", 0, *normres); fflush(stdout);
  }
  if(*normres <= tol) return 0;

  /* start iterations */
  for(it=1; it<=itmax; it++)
  { 
    /* make one SOR step */
    SOR_step(Aline, nlines, x, b, omega);
    /* compute the residual f_i */
    /* f_i = sum_j A_ij x_j - b_i */
    SparseMatrixLines_times_vector(Aline, nlines, x, f);
    for(i=0; i<nlines; i++)  f[i] -= b[i];
    *normres = sqrt(scalarproduct_vectors(f,f, nlines));
    if(pr)
    { 
      printf("matrix_SOR: step %d  residual=%g\n", it, *normres);
      fflush(stdout);
    }
    if(*normres <= tol)  break;
  } /* end iterations loop */

  free(f);
  return it;
}

/* prepare matrix and vectors and then call my own matrix SOR implementation */
int SOR_Iterator(tVarList *vlx, tVarList *vlb, 
                 tVarList *r, tVarList *c1,tVarList *c2,
    	         int itmax, double tol, double *normres,
	         void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	         void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = vlb->grid;
  int bi;
  int line, nlines, nvars, INFO;
  int pr = Getv("GridIterators_verbose", "yes");
  double omega = Getd("GridIterators_SOR_omega");
  double *x;
  double *b;
  tSparseVector **Aline;
  void (*SetMatrixLines)(tSparseVector **Aline,
       void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
       tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
  void (*copy_vl_into_array)(tVarList *vlx, double *x);
  void (*copy_array_into_vl)(double *x, tVarList *vlx);

  /* decide how we order the matrix lines */
  if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_slowly"))
  {
    SetMatrixLines = SetMatrixLines_slowly;
    copy_vl_into_array = copy_varlist_into_array;
    copy_array_into_vl = copy_array_into_varlist;
  }
  else if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_forSortedVars_slowly"))
  {
    SetMatrixLines = SetMatrixLines_forSortedVars_slowly;
    copy_vl_into_array = copy_varlist_into_array_forSortedVars;
    copy_array_into_vl = copy_array_into_varlist_forSortedVars;
  }
  else
    errorexit("SOR_Iterator: unknown GridIterators_SOR_matrix");

  if(pr) printf("sgrid_SOR:\n");

  /* allocate Aline */
  nvars=vlx->n; 
  nlines = 0;
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;
  Aline = calloc(nlines, sizeof(*Aline));
  if(Aline)  { if(pr) printf("allocated %d matrix lines\n", nlines); }
  else       errorexit("no memory for Aline");
  for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();

  /* set Aline */
  SetMatrixLines(Aline, lop, r, vlx, c1, c2, pr);
  if(pr&&0) prSparseVectorArray(Aline,nlines);
  if(Getv("GridIterators_verbose", "very"))
  {
    static unsigned long int id=0;
    char name[1000];
    snprintf(name, 999, "%s/lop_matrix%lu_%.0fs.mtx",
             Gets("outdir"), id, getTimeIn_s());
    write_SparseVectorArray_inMatrixMarketFormat(name, Aline,nlines, 0);
    id++;
  }

  /* allocate x,b */
  x = calloc(nlines, sizeof(double));
  b = calloc(nlines, sizeof(double));
  if(x==NULL || b==NULL) errorexit("sgrid_SOR: no memory for x,b");
     
  /* set x = vlx , b = vlb */
  copy_vl_into_array(vlx, x);
  copy_vl_into_array(vlb, b);

  /* solve A x = b with WTiterator */
  INFO=matrix_SOR(Aline, nlines, x, b, omega, itmax, tol, normres);
                 
  /* set vlx = x */
  copy_array_into_vl(x, vlx);

  /* free matrix Aline and x,b */
  for(line=0; line<nlines; line++)  FreeSparseVector(Aline[line]);
  free(Aline);
  free(x);
  free(b);

  return INFO;
}


/* SOR_Aline contains the matrix set with SetMatrixLines */
void SOR_Preconditioner_from_SOR_Aline(tVarList *vlx, tVarList *vlb,
                                       tVarList *vlc1, tVarList *vlc2)
{
  /* tGrid *grid = vlb->grid; */
  int INFO;
  int pr = Getv("GridIterators_verbose", "yes");
  int itmax = Getd("GridIterators_Preconditioner_itmax");
  double rtol = Getd("GridIterators_Preconditioner_reltol");
  double omega = Getd("GridIterators_SOR_omega");
  double normb, res;
  double *x;
  double *b;
  void (*copy_vl_into_array)(tVarList *vlx, double *x);
  void (*copy_array_into_vl)(double *x, tVarList *vlx);

  /* decide how we order the matrix lines */
  if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_slowly"))
  {
    copy_vl_into_array = copy_varlist_into_array;
    copy_array_into_vl = copy_array_into_varlist;
  }
  else if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_forSortedVars_slowly"))
  {
    copy_vl_into_array = copy_varlist_into_array_forSortedVars;
    copy_array_into_vl = copy_array_into_varlist_forSortedVars;
  }
  else
    errorexit("SOR_Iterator: unknown GridIterators_SOR_matrix");

  if(pr) printf("SOR_Preconditioner_from_SOR_Aline:\n");

  /* allocate x,b */
  x = calloc(SOR_nlines, sizeof(double));
  b = calloc(SOR_nlines, sizeof(double));
  if(x==NULL || b==NULL) errorexit("SOR_Preconditioner_from_SOR_Aline: no memory for x,b");
     
  /* set x = vlx , b = vlb */
  copy_vl_into_array(vlx, x);
  copy_vl_into_array(vlb, b);

  normb = sqrt(scalarproduct_vectors(b,b, SOR_nlines));

  /* solve only if b!=0 (if b=0 => x=0) */
  if(normb>0.0)
  {
    double tol = normb*rtol;

    /* solve A x = b with SOR */
    INFO=matrix_SOR(SOR_Aline, SOR_nlines, x, b, omega, itmax, tol, &res);
  }

  /* set vlx = x */
  copy_array_into_vl(x, vlx);

  /* free x,b */
  free(x);
  free(b);
  if(pr)
  { 
    printf("SOR_Preconditioner_from_SOR_Aline: vector vlx=%p is now set!\n", vlx);
    fflush(stdout);
  }
}


/* do a linear solve with SOR_Preconditioner_from_SOR_Aline, where global
   SOR_Aline is set by temporarily switching to finite differencing */
int linSolve_with_SOR_precon(tVarList *x, tVarList *b, 
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
  int line, nvars, INFO;
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");
  void (*SetMatrixLines)(tSparseVector **Aline,
       void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
       tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);

  /* decide how we order the matrix lines */
  if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_slowly"))
    SetMatrixLines = SetMatrixLines_slowly;
  else if(Getv("GridIterators_SOR_matrix", "SetMatrixLines_forSortedVars_slowly"))
    SetMatrixLines = SetMatrixLines_forSortedVars_slowly;
  else
    errorexit("SOR_Iterator: unknown GridIterators_SOR_matrix");

  if(pr) printf("linSolve_with_SOR_precon\n");

  /* allocate SOR_Aline */
  nvars=x->n; 
  SOR_nlines = 0;
  forallboxes(grid,bi)  SOR_nlines+=(grid->box[bi]->nnodes)*nvars;
  SOR_Aline = calloc(SOR_nlines, sizeof(*SOR_Aline));
  if(SOR_Aline) { if(pr) printf("allocated %d matrix lines\n", SOR_nlines); }
  else      errorexit("no memory for SOR_Aline");
  for(line=0; line<SOR_nlines; line++)  SOR_Aline[line]=AllocateSparseVector();

  if(use_fd)
  {
    /* save current grid in grid_bak and then convert grid to fin. diff. */
    grid_bak = make_empty_grid(grid->nvariables, 0);
    copy_grid_withoutvars(grid, grid_bak, 0);
    convert_grid_to_fd(grid);
    if(pr) printf("Using finite differencing to set matrix SOR_Aline...\n");
  }

  /* set SOR_Aline */
  SetMatrixLines(SOR_Aline, lop, r, x, c1, c2, pr);                
  if(pr&&0) 
    for(line=0; line<SOR_nlines; line++) prSparseVector(SOR_Aline[line]);

  if(use_fd)
  {
    /* restore grid to spectral */
    copy_grid_withoutvars(grid_bak, grid, 0);
    free_grid(grid_bak);
  }

  /* solve A x = b with lsolver and the Precon SOR_Preconditioner_from_SOR_Aline */
  INFO = lsolver(x, b, r,c1,c2, itmax,tol,normres, lop, SOR_Preconditioner_from_SOR_Aline);

  /* free matrix SOR_Aline */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  for(line=0; line<SOR_nlines; line++)  FreeSparseVector(SOR_Aline[line]);
  free(SOR_Aline);
  SOR_nlines=0;

  return INFO;
}

/* do a linear solve with bicgstab and precon_Ap_Ai_Ax_UMFPACK */
int bicgstab_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("bicgstab_with_SOR_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_SOR_precon(x, b, r,c1,c2, bicgstab,
                                  itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_gmres_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_gmres_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_gmres_wrapper_with_SOR_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_SOR_precon(x, b, r,c1,c2, templates_gmres_wrapper,
                                  itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicgstab_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_bicgstab_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicgstab_wrapper_with_SOR_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_SOR_precon(x, b, r,c1,c2, templates_bicgstab_wrapper,
                                  itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_cgs_wrapper and precon_Ap_Ai_Ax_UMFPACK */
int templates_cgs_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_cgs_wrapper_with_SOR_precon: using ");

  /* solve A x = b with bicgstab and the Precon precon_Ap_Ai_Ax_UMFPACK */
  INFO = linSolve_with_SOR_precon(x, b, r,c1,c2, templates_cgs_wrapper,
                                  itmax,tol,normres, lop);
  return INFO;
}
