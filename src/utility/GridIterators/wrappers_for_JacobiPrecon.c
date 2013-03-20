/* wrappers_for_JacobiPrecon.c */
/* Wolfgang Tichy 1/2013 */


#include "sgrid.h"
#include "GridIterators.h"

#ifdef UMFPACK
#include "umfpack.h"
#endif

#ifdef SUITESPARSEQR
#include "SuiteSparseQR_C.h"
#endif


/* blocks of a block diagonal matrix */
typedef struct TBlocks_JacobiPrecon
{
  int nblocks; /* # of blocks */
  int type;    /* 0 means arrays of tSparseVector **, 
                  1 means UMFPACK, 2 means SPQR */
  int *blockdims;           /* array of dims of blocks 0 to nblocks-1 */
  tSparseVector ***Mblock;  /* array of matrices for blocks 0 to nblocks-1 */
  tUMFPACK_A *umfpackA;  /* struct containing all needed for umfpack */
  tSPQR_A *SPQR;         /* struct containing all needed for SPQR*/
} tBlocks_JacobiPrecon;


/* global vars */
/* inverse Diagonal Matrix elements */
double *DiagMinv_JacobiPrecon;

/* blocks of a block diagonal matrix */
tBlocks_JacobiPrecon Blocks_JacobiPrecon;


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

/* solver for one block */
void Matrix_BlockJacobi_Solve(tBlocks_JacobiPrecon Blocks, int i,
                              double *x, double *b, int ncols)
{
  int pr = 0; // Getv("GridIterators_verbose", "yes");

  if(Blocks.type==2)
    SuiteSparseQR_solve_from_tSPQR_A_x_b(Blocks.SPQR[i], x, b, pr);                                             
  else
    umfpack_dl_solve_from_tUMFPACK_A_x_b(Blocks.umfpackA[i], x, b, pr);
}

/* Blocks_JacobiPrecon contains a block diagonal matrix set with
   SetMatrixColumns_ForOneVarInOneSubBox_slowly */
void BlockJacobi_Preconditioner_from_Blocks(tVarList *vlx, tVarList *vlb,
                                            tVarList *vlc1, tVarList *vlc2)
{
  tGrid *grid = vlx->grid;
  int nsb1 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb1");
  int nsb2 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb2");
  int nsb3 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb3");
  int blocki;

  /* loop over vars, boxes and subboxes, i.e. the blocks we have */
  SGRID_TOPLEVEL_Pragma(omp parallel for)
  forallVarsBoxesAndSubboxes_defIndices(vlx, blocki, vi,bi, sbi,sbj,sbk,
                                        nsb1,nsb2,nsb3)
  {
    tBox *box = grid->box[bi];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i1,i2, j1,j2, k1,k2; 
    double *x;
    double *b;
    int ncols;
    /* blocki = sbi + nsb1 * sbj + nsb1*nsb2 * sbk
                + nsb1*nsb2*nsb3 * bi + nsb1*nsb2*nsb3*(grid->nboxes) * vi 
    int blocki = sbi + nsb1 * (sbj + nsb2 * (sbk + 
                 nsb3 * (bi + (grid->nboxes) * vi))); */

    /* allocate mem for x,b */
    IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
    x = dmalloc((i2-i1)*(j2-j1)*(k2-k1));
    b = dmalloc((i2-i1)*(j2-j1)*(k2-k1));

    /* set b */
    copy_varlistCompInSubbox_into_array(vlb, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3, b);
    /* solve for var vi in box bi in subbox sbi,sbj,sbk,
       i.e. solve for x in this block */
    ncols = Blocks_JacobiPrecon.blockdims[blocki];
    //Acol  = Blocks_JacobiPrecon.Mblock[blocki];
    Matrix_BlockJacobi_Solve(Blocks_JacobiPrecon, blocki, x, b, ncols);

    /* set vlx */
    copy_array_into_varlistCompInSubbox(x, vlx, vi, bi,
                                        sbi,sbj,sbk, nsb1,nsb2,nsb3);
    /* free arrays */
    free(x); free(b);
  } End_forallVarsBoxesAndSubboxes_defIndices
}


/* do a linear solve with BlockJacobi_Preconditioner_from_Blocks, where global
   Blocks_JacobiPrecon contains the blocks. The blocks can be set 
   by temporarily switching to finite differencing */
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
  int i, blocki, nblocks;
  int pr = Getv("GridIterators_verbose", "yes");
  int dropbelow = Getd("GridIterators_setABStozero_below");
  int INFO;
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");
  int nsb1 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb1");
  int nsb2 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb2");
  int nsb3 = Getd("GridIterators_Preconditioner_BlockJacobi_nsb3");
  int type;
  
  if(Getv("GridIterators_Preconditioner_type", "SPQR")) /* use SPQR */
    type=2;
  else /* use umfpack */
    type=1;

  if(pr) printf("linSolve_with_BlockJacobi_precon\n");

  /* allocate memory for blocks in Blocks_JacobiPrecon struct */
  nblocks = (b->n)*(grid->nboxes)*nsb1*nsb2*nsb3;
  Blocks_JacobiPrecon.nblocks = nblocks;
  Blocks_JacobiPrecon.type = type;
  Blocks_JacobiPrecon.blockdims = (int *) calloc(nblocks, sizeof(int));
  Blocks_JacobiPrecon.Mblock 
    = (tSparseVector ***) calloc(nblocks, sizeof(tSparseVector **));
  Blocks_JacobiPrecon.umfpackA
    = (tUMFPACK_A *) calloc(nblocks, sizeof(tUMFPACK_A));
  Blocks_JacobiPrecon.SPQR
    = (tSPQR_A *) calloc(nblocks, sizeof(tSPQR_A));

  /* convert grid to fin. diff.? */
  if(use_fd)
  {
    /* save current grid in grid_bak and then convert grid to fin. diff. */
    grid_bak = make_empty_grid(grid->nvariables, 0);
    copy_grid_withoutvars(grid, grid_bak, 0);
    convert_grid_to_fd(grid);
    if(pr) printf(" Using finite differencing to set all matrix blocks.\n");
  }

  if(pr) prTimeIn_s("Time BEFORE setting up blocks: ");

  /* loop over boxes and vars */
#ifndef LEVEL6_Pragmas
  SGRID_TOPLEVEL_Pragma(omp parallel for)
#endif
  forallVarsBoxesAndSubboxes_defIndices(b, blocki, vi,bi, sbi,sbj,sbk,
                                        nsb1,nsb2,nsb3)
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
    tSPQR_A SPQR;

    /* allocate Acol for each block */
    IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
    ncols = (i2-i1)*(j2-j1)*(k2-k1);
    Acol=AllocateSparseVectorArray(ncols);
    if(Acol) { if(0) printf(" allocated %d matrix columns for block%d\n",
                             ncols, blocki); fflush(stdout); }
    else       errorexit(" no memory for Acol");

    /* set Acol */
    if(pr) 
      printf(" setting %d matrix columns for block%d...\n", ncols, blocki);
    fflush(stdout);
    SetMatrixColumns_ForOneVarInOneSubBox_slowly(Acol, vi, bi,
                                                 sbi,sbj,sbk, nsb1,nsb2,nsb3,
                                                 lop, r, x, c1, c2, 0);
    if(pr&&0) prSparseVectorArray(Acol,ncols); 

    /* count number of entries in sparse matrix */
    nz = 0;
    for(j = 0; j < ncols; j++) nz+=Acol[j]->entries;

    /* allocate memory for matrix and set UMFPACK's Ap,Ai,Ax or SPQR's A */
    if(Blocks_JacobiPrecon.type==2) /* if we use SPQR */
    {
#ifdef SUITESPARSEQR
      cholmod_sparse *A;
      allocate_and_init_tSPQR_A_struct(&SPQR, ncols, nz, 0);
      /* make Ap,Ai,Ax point inside SPQR.A */
      A  = (cholmod_sparse *) SPQR.A;
      Ap = (LONGINT *) A->p;
      Ai = (LONGINT *) A->i;
      Ax = (double *)  A->x;
      /* set Ap,Ai,Ax inside SPQR.A */
      set_umfpack_dl_matrix_from_columns(Ap,Ai,Ax, Acol, ncols, dropbelow, 0);
#else
      allocate_and_init_tSPQR_A_struct(&SPQR, ncols, nz, 1);
#endif
      /* now let Ap,Ai,Ax point to NULL again,
         so nothing bad happens if we free them below! */
      Ap = NULL;
      Ai = NULL;
      Ax = NULL;
    }
    else
    {
      SPQR.A = NULL; /* this signals that nothing is in the SPQR struct */
      allocate_umfpack_dl_matrix(&Ap, &Ai, &Ax, ncols, nz);

      /* set Ap,Ai,Ax */
      set_umfpack_dl_matrix_from_columns(Ap,Ai,Ax, Acol, ncols, dropbelow, 0);
    }

    /* since Acol is now in Ap,Ai,Ax or SPQR we free Acol */
    FreeSparseVectorArray(Acol, ncols);
    Acol=NULL;

    /* set each entry in struct */
    Blocks_JacobiPrecon.blockdims[blocki] = ncols;
    Blocks_JacobiPrecon.Mblock[blocki]    = Acol;

#ifdef UMFPACK
    Blocks_JacobiPrecon.umfpackA[blocki].sys = UMFPACK_A;
#endif
    Blocks_JacobiPrecon.umfpackA[blocki].Ap  = Ap;
    Blocks_JacobiPrecon.umfpackA[blocki].Ai  = Ai;
    Blocks_JacobiPrecon.umfpackA[blocki].Ax  = Ax;
    Blocks_JacobiPrecon.umfpackA[blocki].Numeric = NULL; /* it's set below */

    Blocks_JacobiPrecon.SPQR[blocki] = SPQR; /* the QR inside this is set below */

  } End_forallVarsBoxesAndSubboxes_defIndices
  if(pr)
  {
    printf("linSolve_with_BlockJacobi_precon: created %d blocks.\n", nblocks);
    prTimeIn_s("Time AFTER setting up blocks: ");
  }
  if(use_fd)
  {
    /* restore grid to spectral */
    copy_grid_withoutvars(grid_bak, grid, 0);
    free_grid(grid_bak);
  }

  if(pr) prTimeIn_s("Time BEFORE LU or QR factorization: ");

  /* do LU or QR factorization with umfpackA or SPQR */
#ifndef MEMORY_EFFICIENT
  SGRID_TOPLEVEL_Pragma(omp parallel for)
#endif
  for(i=0; i<nblocks; i++)
  {
    if(Blocks_JacobiPrecon.type==2) /* if we use SPQR */
    {
      printf(" SuiteSparseQR_C_factorize_tSPQR_A in block%d...\n", i);
      fflush(stdout);
      SuiteSparseQR_C_factorize_tSPQR_A(&(Blocks_JacobiPrecon.SPQR[i]), 0);
    }
    else
    {
      printf(" umfpack_dl_numeric_from_tUMFPACK_A in block%d...\n", i);
      fflush(stdout);
      umfpack_dl_numeric_from_tUMFPACK_A(&(Blocks_JacobiPrecon.umfpackA[i]),
                                         Blocks_JacobiPrecon.blockdims[i], 0);
    }
  }
  if(pr) prTimeIn_s("Time AFTER LU or QR factorization: ");

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
    free(Blocks_JacobiPrecon.umfpackA[blocki].Ap);
    free(Blocks_JacobiPrecon.umfpackA[blocki].Ai);
    free(Blocks_JacobiPrecon.umfpackA[blocki].Ax);
#ifdef UMFPACK
    umfpack_dl_free_numeric(&(Blocks_JacobiPrecon.umfpackA[blocki].Numeric));
#endif
    free_tSPQR_A_struct(&(Blocks_JacobiPrecon.SPQR[blocki]));
  }
  free(Blocks_JacobiPrecon.blockdims);
  free(Blocks_JacobiPrecon.Mblock);
  free(Blocks_JacobiPrecon.umfpackA);
  free(Blocks_JacobiPrecon.SPQR);

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
