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
  errorexit("SuiteSparseQR_solve_fromAcolumns: in order to compile with SuiteSparseQR use\n"
            "MyConfig with\n"
            "DFLAGS += -DSUITESPARSEQR\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lspqr -lcholmod\n"
            "you may also need something like\n"
            "CLINKER = g++\n"
            "or\n"
            "CLINKER = icpc\n");
#endif
  return INFO;
}
