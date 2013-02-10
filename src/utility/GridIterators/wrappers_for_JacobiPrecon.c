/* wrappers_for_JacobiPrecon.c */
/* Wolfgang Tichy 1/2013 */


#include "sgrid.h"
#include "GridIterators.h"


/* global vars */
/* inverse Diagonal Matrix elements */
double *DiagMinv_JacobiPrecon;

/* blocks of a block diagonal matrix */
struct
{
  int nblocks; /* # of blocks */
  int type;    /* 0 means arrays of tSparseVector **, 
                  1 means UMFPACK Ap,Ai,Ax arrays */
  int *blockdims;           /* array of dims of blocks 0 to nblocks-1 */
  tSparseVector ***Mblock;  /* array of matrices for blocks 0 to nblocks-1 */
  LONGINT **Ap; /* array of Ap's */ // currently not used
  LONGINT **Ai; /* array of Ai's */ // currently not used
  double  **Ax; /* array of Ax's */ // currently not used
} Blocks_JacobiPrecon;


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
        xp[i] = bp[i]*DiagMinv_JacobiPrecon[line];
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
   DiagMinv_JacobiPrecon is set by temporarily switching to finite differencing */
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

  /* allocate memory for diagonal of matrix in DiagMinv_JacobiPrecon */
  DiagMinv_JacobiPrecon = calloc(ncols, sizeof(*DiagMinv_JacobiPrecon));

  /* now set DiagMinv_JacobiPrecon to diagonal of matrix Acol */
  for(col=0; col<ncols; col++)
    for(ent = 0; ent < Acol[col]->entries; ent++)
      if(Acol[col]->pos[ent] == col) 
      {
        double DiagM = Acol[col]->val[ent];
        if(DiagM==0.0) errorexit("DiagM is singular!!!");
        DiagMinv_JacobiPrecon[col] = 1.0/DiagM;
        break;
      }
  
  /* solve A x = b with lsolver and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = lsolver(x, b, r,c1,c2, itmax,tol,normres, lop, Jacobi_Preconditioner_from_DiagM);

  /* free matrix DiagMinv_JacobiPrecon */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  free(DiagMinv_JacobiPrecon);

  /* free matrix Acol */
  for(col=0; col<ncols; col++)  FreeSparseVector(Acol[col]);
  free(Acol);

  return INFO;
}

/* do a linear solve with bicgstab and Jacobi_Preconditioner_from_DiagM */
int bicgstab_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("bicgstab_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, bicgstab,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_gmres_wrapper and Jacobi_Preconditioner_from_DiagM */
int templates_gmres_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_gmres_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with GMRES and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_gmres_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicgstab_wrapper and Jacobi_Preconditioner_from_DiagM */
int templates_bicgstab_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicgstab_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_bicgstab_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_cgs_wrapper and Jacobi_Preconditioner_from_DiagM */
int templates_cgs_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_cgs_wrapper_with_Jacobi_precon: using ");

  /* solve A x = b with CGS and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = linSolve_with_Jacobi_precon(x, b, r,c1,c2, templates_cgs_wrapper,
                                     itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_qmr_wrapper and Jacobi_Preconditioner_from_DiagM */
int templates_qmr_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_qmr_wrapper_with_Jacobi_precon:\n");

  /* solve A x = b with qmr and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = templates_qmr_wrapper(x, b, r,c1,c2, itmax,tol,normres,
                               lop, Jacobi_Preconditioner_from_DiagM);
  return INFO;
}

/* do a linear solve with templates_bicg_wrapper and Jacobi_Preconditioner_from_DiagM */
int templates_bicg_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicg_wrapper_with_Jacobi_precon:\n");

  /* solve A x = b with bicg and the Precon Jacobi_Preconditioner_from_DiagM */
  INFO = templates_bicg_wrapper(x, b, r,c1,c2, itmax,tol,normres,
                                lop, Jacobi_Preconditioner_from_DiagM);
  return INFO;
}


/* ********************************************************************* */
/* Block diagonal Preconditioners                                        */
/* ********************************************************************* */

/* example of a solver for a block */
void Matrix_BlockJacobi_Solve(LONGINT *Ap, LONGINT *Ai, double *Ax,
                              double *x, double *b, int ncols)
{
  int pr = 0; // Getv("GridIterators_verbose", "yes");

  umfpack_dl_solve_from_Ap_Ai_Ax_x_b(Ap,Ai,Ax, x, b, ncols, pr);
}

/* Blocks_JacobiPrecon contains a block diagonal matrix set with 
   SetMatrixColumns_ForOneVarInOneBox_slowly */
void BlockJacobi_Preconditioner_from_Blocks(tVarList *vlx, tVarList *vlb,
                                            tVarList *vlc1, tVarList *vlc2)
{
  tGrid *grid = vlx->grid;
  int bi, vi;
  int sbi,sbj,sbk;
  int nsb1 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb1");
  int nsb2 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb2");
  int nsb3 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb3");

  /* loop over vars, boxes and subboxes */
  /* we want: #pragma omp parallel for collapse(5) 
     but icc 10.1 says: syntax error in omp clause */
  SGRID_TOPLEVEL_Pragma(omp parallel for collapse(5))
  forallVarsBoxesAndSubboxes(vlx, vi,bi, sbi,sbj,sbk, nsb1,nsb2,nsb3)
  {
    tBox *box = grid->box[bi];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i1,i2, j1,j2, k1,k2; 
    double *x;
    double *b;
    //tSparseVector **Acol;
    int ncols;
    LONGINT *Ap;
    LONGINT *Ai;
    double *Ax;
    /* blocki = sbi + nsb1 * sbj + nsb1*nsb2 * sbk
                + nsb1*nsb2*nsb3 * bi + nsb1*nsb2*nsb3*(grid->nboxes) * vi */
    int blocki = sbi + nsb1 * (sbj + nsb2 * (sbk + 
                 nsb3 * (bi + (grid->nboxes) * vi)));

    /* allocate mem for x,b */
    IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
    x = dmalloc((i2-i1)*(j2-j1)*(k2-k1));
    b = dmalloc((i2-i1)*(j2-j1)*(k2-k1));

    /* set x,b */
    copy_varlistCompInSubbox_into_array(vlx, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3, x);
    copy_varlistCompInSubbox_into_array(vlb, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3, b);
    /* solve for var vi in box bi */
    ncols = Blocks_JacobiPrecon.blockdims[blocki];
    //Acol  = Blocks_JacobiPrecon.Mblock[blocki];
    Ap    = Blocks_JacobiPrecon.Ap[blocki];
    Ai    = Blocks_JacobiPrecon.Ai[blocki];
    Ax    = Blocks_JacobiPrecon.Ax[blocki];
    Matrix_BlockJacobi_Solve(Ap,Ai,Ax, x, b, ncols);

    /* set vlx, vlb */
    copy_array_into_varlistCompInSubbox(x, vlx, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3);
    copy_array_into_varlistCompInSubbox(b, vlb, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3);
    /* free arrays */
    free(x); free(b);
  }
}


/* do a linear solve with Jacobi_Preconditioner_from_DiagM, where global
   DiagMinv_JacobiPrecon is set by temporarily switching to finite differencing */
int linSolve_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
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
  int bi, vi, blocki, nblocks;
  int pr = Getv("GridIterators_verbose", "yes");
  int dropbelow = Getd("GridIterators_setABStozero_below");
  int INFO;
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");
  int sbi,sbj,sbk;
  int nsb1 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb1");
  int nsb2 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb2");
  int nsb3 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb3");

  if(pr) printf("linSolve_with_BlockJacobi_precon\n");

  /* allocate memory for blocks in Blocks_JacobiPrecon struct */
  nblocks = (b->n)*(grid->nboxes)*nsb1*nsb2*nsb3;
  Blocks_JacobiPrecon.nblocks = nblocks;
  Blocks_JacobiPrecon.type = 0; /* set type to 0 for now */
  Blocks_JacobiPrecon.blockdims = (int *) calloc(nblocks, sizeof(int));
  Blocks_JacobiPrecon.Mblock 
    = (tSparseVector ***) calloc(nblocks, sizeof(tSparseVector **));
  Blocks_JacobiPrecon.Ap 
    = (LONGINT **) calloc(nblocks, sizeof(LONGINT *));
  Blocks_JacobiPrecon.Ai 
    = (LONGINT **) calloc(nblocks, sizeof(LONGINT *));
  Blocks_JacobiPrecon.Ax 
    = (double **) calloc(nblocks, sizeof(double *));
  
  /* loop over boxes and vars */
  blocki=0;
  forallVarsBoxesAndSubboxes(b, vi,bi, sbi,sbj,sbk, nsb1,nsb2,nsb3)
  {
    tBox *box = grid->box[bi];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i1,i2, j1,j2, k1,k2;
    int col, ncols, j;
    tSparseVector **Acol;
    LONGINT nz;
    LONGINT *Ap;
    LONGINT *Ai;
    double *Ax;

    /* allocate Acol for each block */
    IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
    ncols = (i2-i1)*(j2-j1)*(k2-k1);
    Acol=AllocateSparseVectorArray(ncols);
    if(Acol) { if(pr) printf("allocated %d matrix columns\n", ncols); }
    else       errorexit("no memory for Acol");

    if(use_fd)
    {
      /* save current grid in grid_bak and then convert grid to fin. diff. */
      grid_bak = make_empty_grid(grid->nvariables, 0);
      copy_grid_withoutvars(grid, grid_bak, 0);
      convert_grid_to_fd(grid);
      if(pr) printf("Using finite differencing to set matrix Acol...\n");
    }

    /* set Acol */
    SetMatrixColumns_ForOneVarInOneSubBox_slowly(Acol, vi, bi,
                                                 sbi,sbj,sbk, nsb1,nsb2,nsb3,
                                                 lop, r, x, c1, c2, pr);
    if(pr&&0) 
      for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

    if(use_fd)
    {
      /* restore grid to spectral */
      copy_grid_withoutvars(grid_bak, grid, 0);
      free_grid(grid_bak);
    }

    /* set Ap,Ai,Ax */
    /* count number of entries in sparse matrix */
    nz = 0;
    for(j = 0; j < ncols; j++) nz+=Acol[j]->entries;
    /* allocate memory for Ap,Ai,Ax for matrix */
    allocate_umfpack_dl_matrix(&Ap, &Ai, &Ax, ncols, nz);
    /* set Ap,Ai,Ax and solve */
    set_umfpack_dl_matrix_from_columns(Ap,Ai,Ax, Acol, ncols, dropbelow, 0);

    /* since Acol is now in Ap,Ai,Ax, we free Acol */
    FreeSparseVectorArray(Acol, ncols);
    Acol=NULL;

    /* set each entry in struct */
    Blocks_JacobiPrecon.blockdims[blocki] = ncols;
    Blocks_JacobiPrecon.Mblock[blocki]    = Acol;
    Blocks_JacobiPrecon.Ap[blocki] = Ap;
    Blocks_JacobiPrecon.Ai[blocki] = Ai;
    Blocks_JacobiPrecon.Ax[blocki] = Ax;

    blocki++;
  }
  
  /* solve A x = b with lsolver and BlockJacobi_Preconditioner_from_Blocks */
  INFO = lsolver(x, b, r,c1,c2, itmax,tol,normres, lop, 
                 BlockJacobi_Preconditioner_from_Blocks);

  /* free matrix blocks */
  /* maybe we could save and reuse this matrix for several lin solves,
     if the linear solve succeeds in driving down the residual of lop ???*/
  for(blocki=0; blocki<Blocks_JacobiPrecon.nblocks; blocki++)
  {
    FreeSparseVectorArray(Blocks_JacobiPrecon.Mblock[blocki],
                          Blocks_JacobiPrecon.blockdims[blocki]);
    free(Blocks_JacobiPrecon.Ap[blocki]);
    free(Blocks_JacobiPrecon.Ai[blocki]);
    free(Blocks_JacobiPrecon.Ax[blocki]);
  }
  free(Blocks_JacobiPrecon.blockdims);
  free(Blocks_JacobiPrecon.Mblock);
  free(Blocks_JacobiPrecon.Ap);
  free(Blocks_JacobiPrecon.Ai);
  free(Blocks_JacobiPrecon.Ax);

  return INFO;
}

/* do a linear solve with bicgstab and BlockJacobi_Preconditioner_from_Blocks */
int bicgstab_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("bicgstab_with_BlockJacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2, bicgstab,
                                          itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_gmres_wrapper and BlockJacobi_Preconditioner_from_Blocks */
int templates_gmres_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_gmres_wrapper_with_BlockJacobi_precon: using ");

  /* solve A x = b with GMRES and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2,
                                          templates_gmres_wrapper,
                                          itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicgstab_wrapper and BlockJacobi_Preconditioner_from_Blocks */
int templates_bicgstab_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicgstab_wrapper_with_BlockJacobi_precon: using ");

  /* solve A x = b with bicgstab and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2,
                                          templates_bicgstab_wrapper,
                                          itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_cgs_wrapper and BlockJacobi_Preconditioner_from_Blocks */
int templates_cgs_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_cgs_wrapper_with_BlockJacobi_precon: using ");

  /* solve A x = b with CGS and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2,
                                          templates_cgs_wrapper,
                                          itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_qmr_wrapper and BlockJacobi_Preconditioner_from_Blocks */
int templates_qmr_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_qmr_wrapper_with_BlockJacobi_precon:\n");

  /* solve A x = b with qmr and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2,
                                          templates_qmr_wrapper,
                                          itmax,tol,normres, lop);
  return INFO;
}

/* do a linear solve with templates_bicg_wrapper and BlockJacobi_Preconditioner_from_Blocks */
int templates_bicg_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  int pr = Getv("GridIterators_verbose", "yes");
  int INFO;
  if(pr) printf("templates_bicg_wrapper_with_BlockJacobi_precon:\n");

  /* solve A x = b with bicg and the Precon BlockJacobi_Preconditioner_from_Blocks */
  INFO = linSolve_with_BlockJacobi_precon(x, b, r,c1,c2,
                                          templates_bicg_wrapper,
                                          itmax,tol,normres, lop);
  return INFO;
}
