/* SuiteSparseQR_C_interface.c */
/* Wolfgang Tichy 1/2013 */


#include "sgrid.h"
#include "LinAlg.h"

#ifdef SUITESPARSEQR
#include "SuiteSparseQR_C.h"
#define PrintErrorCodesAndExit  \
  { printf("CHOLMOD returned INFO=%d\n", INFO); \
    printf("some common return values:\n"); \
    printf(" CHOLMOD_OK=%d\n", CHOLMOD_OK); \
    printf(" CHOLMOD_NOT_INSTALLED=%d\n", CHOLMOD_NOT_INSTALLED); \
    printf(" CHOLMOD_OUT_OF_MEMORY=%d\n", CHOLMOD_OUT_OF_MEMORY); \
    printf(" CHOLMOD_TOO_LARGE=%d\n", CHOLMOD_TOO_LARGE); \
    printf(" CHOLMOD_INVALID=%d\n", CHOLMOD_INVALID); \
    printf(" CHOLMOD_GPU_PROBLEM= maybe -5\n"); \
    printf(" CHOLMOD_NOT_POSDEF=%d\n", CHOLMOD_NOT_POSDEF); \
    printf(" CHOLMOD_DSMALL=%d\n", CHOLMOD_DSMALL); \
    printf("for more info do:\n"); \
    printf(" grep CHOLMOD_ /usr/include/suitesparse/cholmod_core.h\n"); \
    fflush(stdout); \
    if(INFO<0) errorexiti("CHOLMOD returned INFO=%d", INFO); }
#else
#define PrintErrorCodesAndExit  errorexit("SuiteSparseQR is not compiled in")
#endif

#define CompileSuiteSparseQR \
  errorexit("In order to compile with SuiteSparseQR use\n" \
            "MyConfig with\n" \
            "DFLAGS += -DSUITESPARSEQR\n" \
            "SPECIALINCS += -I/usr/include/suitesparse\n" \
            "SPECIALLIBS += -lspqr -lcholmod\n" \
            "you may also need something like\n" \
            "CLINKER = g++\n" \
            "or\n" \
            "CLINKER = icpc\n")


/***************************************************************************/
/* some helper routines                                                    */
/***************************************************************************/

#ifdef SUITESPARSEQR
/* get a sparse matrix in cholmod_sparse format from the array in Acol */
cholmod_sparse *get_cholmod_sparse_fromAcolumns(tSparseVector **Acol,
                                                LONGINT nrows, LONGINT ncols,
                                                double dropbelow, int pr,
                                                cholmod_common *cc)
{
  cholmod_sparse *A;
  int sorted = 1; /* we sort row indices in each column */
  int packed = 1; /* we use the packed format */
  int stype = 0;  /* we have an unsymmetric matrix */
  int xtype = CHOLMOD_REAL; /* our matrix has real numbers as entries */
  LONGINT j, nz;
  LONGINT *Ap;
  LONGINT *Ai;
  double *Ax;

  /* count the number of entries nz in Acol array */
  nz = 0;
  for(j = 0; j < ncols; j++) nz += Acol[j]->entries;

  /* allocate A */
  A = cholmod_l_allocate_sparse(nrows,ncols, nz, 
                                sorted, packed, stype, xtype, cc);
  /* set relevant arrays inside A */
  Ap = (LONGINT *) A->p;
  Ai = (LONGINT *) A->i;
  Ax = (double *) A->x;
  set_umfpack_dl_matrix_from_columns(Ap,Ai,Ax, Acol,ncols, dropbelow, 0);

  if(!cholmod_l_check_sparse(A, cc))
    errorexit("cholmod_l_check_sparse(A, cc) failed");

  if(pr)
  {
    cholmod_l_print_sparse(A, "A" ,cc);
    printf("get_cholmod_sparse_fromAcolumns: the sparse %ld*%ld matrix "
           "A->p[%ld]=%ld, A->i=%p, A->x=%p is now set!\n",
           nrows, ncols, ncols, Ap[ncols], Ai, Ax);
  }
  return A;
}
#endif

/* allocate and init a tSPQR_A struct, 
   needs to be freed later with free_tSPQR_A_struct */
int allocate_and_init_tSPQR_A_struct(tSPQR_A *SPQR, LONGINT ncols, 
                                     LONGINT nz, int pr)
{
  int sorted = 1; /* we sort row indices in each column */
  int packed = 1; /* we use the packed format */
  int stype = 0;  /* we have an unsymmetric matrix */
#ifdef SUITESPARSEQR
  int xtype = CHOLMOD_REAL; /* our matrix has real numbers as entries */
  cholmod_sparse *A;
  cholmod_common *cc;
  
  /* init cc */
  cc = (cholmod_common *) malloc(sizeof(cholmod_common));
  cholmod_l_start(cc);

  /* get mem for matrix A */
  A = cholmod_l_allocate_sparse(ncols,ncols, nz, sorted,
                                packed, stype, xtype, cc);
  /* set local SPQR struct */
  SPQR->sys      = SPQR_RX_EQUALS_B;
  SPQR->ordering = SPQR_ORDERING_DEFAULT;
  SPQR->tol      = SPQR_DEFAULT_TOL;
  SPQR->A  = (void *) A;
  SPQR->QR = NULL;        /* need to set this later */  
  SPQR->cc = (void *) cc;
#else
  SPQR->sys      = 0;
  SPQR->ordering = 7;
  SPQR->tol      = -2;
  SPQR->A  = NULL;
  SPQR->QR = NULL;
  SPQR->cc = NULL;
#endif
  if(pr)
  { 
    printf("allocate_and_init_tSPQR_A_struct: "
           "A=%p QR=%p cc=%p\n", SPQR->A, SPQR->QR, SPQR->cc);
    fflush(stdout);
  }
  return 0;
}

/* free a tSPQR_A struct */
int free_tSPQR_A_struct(tSPQR_A *SPQR_A)
{
#ifdef SUITESPARSEQR
  cholmod_sparse *A  = (cholmod_sparse *) SPQR_A->A;
  cholmod_common *cc = (cholmod_common *) SPQR_A->cc;
  SuiteSparseQR_C_factorization *QR 
    = (SuiteSparseQR_C_factorization *) SPQR_A->QR;

  if(SPQR_A->A !=NULL)
  {
    /* free everything and finish CHOLMOD */
    SuiteSparseQR_C_free(&QR, cc);
    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
    free(cc);
  }
  return 0;
#else
  return -1;
#endif
}

/***************************************************************************/
/* main solver routines                                                    */
/***************************************************************************/

/* solve A X = B with SuiteSparseQR's SuiteSparseQR_C_backslash_default
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int SuiteSparseQR_solve_fromAcolumns(tSparseVector **Acol,
                                     tVarList *vlx, tVarList *vlb,
                                     double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int INFO=-6662442;
  int nvars=vlx->n;
  int nlines=0;
#ifdef SUITESPARSEQR
  cholmod_common Common, *cc;
  cholmod_sparse *A;
  cholmod_dense *X, *B;

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* start CHOLMOD */
  cc = &Common;
  cholmod_l_start(cc);

  /* allocate and set A */
  A = get_cholmod_sparse_fromAcolumns(Acol, nlines,nlines, dropbelow, pr, cc);

  /* allocate B */
  B = cholmod_l_allocate_dense(nlines,1, nlines, A->xtype, cc);

  /* set B = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlb->n; j++)
      {
        double *bp = box->v[vlb->index[j]];
        double *b = (double *) B->x;
        b[line] = bp[i];
        line++;
      }
  }
  if(!cholmod_l_check_dense(B ,cc))
    errorexit("cholmod_l_check_dense(B, cc) failed");
  if(pr) cholmod_l_print_dense(B, "B", cc);

  if(pr)
  { printf("SuiteSparseQR_solve_fromAcolumns: SuiteSparseQR_C_backslash_default\n"); fflush(stdout); }

  /* Solve, i.e. in matlab notation do: X = A\B */
  X = SuiteSparseQR_C_backslash_default(A, B, cc);
  INFO=cc->status;
  if(pr)
  { 
    cholmod_l_print_dense(X, "X", cc);
    printf(" -> INFO=%d\n", INFO);
    fflush(stdout);
  }
  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = X */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *xp = box->v[vlx->index[j]];
        double *x = (double *) X->x;
        xp[i] = x[line];
        line++;
      }
  }

  /* free everything and finish CHOLMOD */
  cholmod_l_free_sparse(&A, cc);
  cholmod_l_free_dense(&X, cc);
  cholmod_l_free_dense(&B, cc);
  cholmod_l_finish(cc);

#else
  CompileSuiteSparseQR;
#endif
  return INFO;
}


/* get QR factorization and put it into tSPQR_A struct */
int SuiteSparseQR_C_factorize_tSPQR_A(tSPQR_A *SPQR_A, int pr)
{
  int INFO=-1;
#ifdef SUITESPARSEQR
  cholmod_sparse *A  = (cholmod_sparse *) SPQR_A->A;
  cholmod_common *cc = (cholmod_common *) SPQR_A->cc;
  SuiteSparseQR_C_factorization *QR 
                       = (SuiteSparseQR_C_factorization *) SPQR_A->QR;

  if(!cholmod_l_check_sparse(A, cc) || pr)
  {
    printf("SuiteSparseQR_C_factorize_tSPQR_A: input is:\n"
           "A=%p QR=%p cc=%p cc->status=%d\n", A, QR, cc, cc->status);
    cholmod_l_print_sparse(A, "A", cc);
  }

  QR = SuiteSparseQR_C_factorize(SPQR_A->ordering, SPQR_A->tol, A, cc);
  SPQR_A->QR = (void *) QR;
  INFO=cc->status;
  if(pr || INFO!=0)
  { 
    printf("SuiteSparseQR_C_factorize_tSPQR_A: output is:\n"
           "A=%p QR=%p cc=%p cc->status=%d\n", A, QR, cc, cc->status);
    printf("SuiteSparseQR_C_factorize_tSPQR_A -> INFO=%d\n", INFO);
    fflush(stdout);
  }
  if(INFO!=0)
    PrintErrorCodesAndExit;
#else
  CompileSuiteSparseQR;
#endif
  return INFO;
}

/* solve A x = b with umfpack's SuiteSparseQR's SuiteSparseQR_C_solve 
   for a matrix that is already saved in tSPQR_A SPQR_A with SPQR_A.QR already
   set with prior calls to SuiteSparseQR_C_factorize_tSPQR_A
   with x and b in ordinary double arrays */
int SuiteSparseQR_solve_from_tSPQR_A_x_b(tSPQR_A SPQR_A,
                                         double *x, double *b, int pr)
{
  int i;
  int INFO=-6662442;
  int nlines;
#ifdef SUITESPARSEQR
  SuiteSparseQR_C_factorization *QR = (SuiteSparseQR_C_factorization *) SPQR_A.QR;
  cholmod_sparse *A  = (cholmod_sparse *) SPQR_A.A;
  cholmod_common *cc = (cholmod_common *) SPQR_A.cc;
  cholmod_dense *X, *B;

  /* figure out number of lines */
  nlines = A->nrow; 

  /* allocate X, B */
  B = cholmod_l_allocate_dense(nlines,1, nlines, A->xtype, cc);

  /* set B = b */
  for(i=0; i<nlines; i++)
  {
     double *Bx = (double *) B->x;
     Bx[i] = b[i];
  }
  if(!cholmod_l_check_dense(B ,cc))
    errorexit("cholmod_l_check_dense(B, cc) failed");
  if(pr)
  {
    printf("SuiteSparseQR_solve_from_tSPQR_A_x_b: input is:\n"
           "A=%p QR=%p cc=%p cc->status=%d B=%p\n", A, QR, cc, cc->status, B);
    cholmod_l_print_sparse(A, "A", cc);
    cholmod_l_print_dense(B, "B", cc);
  }

  /* Solve, using prviously computed QR */
  X = SuiteSparseQR_C_solve(SPQR_A.sys, QR, B, cc);
  INFO=cc->status;
  if(pr || INFO!=0)
  { 
    printf("SuiteSparseQR_solve_from_tSPQR_A_x_b: output is:\n"
           "A=%p QR=%p cc=%p cc->status=%d B=%p => X=%p\n",
           A, QR, cc, cc->status, B, X);
    cholmod_l_print_dense(X, "X", cc);
    printf(" -> INFO=%d\n", INFO);
    fflush(stdout);
  }
  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set x = X */
  for(i=0; i<nlines; i++)
  {
    double *Xx = (double *) X->x;
    x[i] = Xx[i];
  }

  /* free X, B */
  cholmod_l_free_dense(&X, cc);
  cholmod_l_free_dense(&B, cc);

#else
  CompileSuiteSparseQR;
#endif
  return INFO;
}
