/* umfpack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"

#ifdef UMFPACK
#include "umfpack.h"
#define PrintErrorCodesAndExit  \
  { printf("umfpack_dl_symbolic returned INFO1=%d\n", INFO1); \
    printf("umfpack_dl_numeric returned INFO2=%d\n", INFO2);  \
    printf("umfpack_dl_solve returned INFO=%d\n", INFO); \
    printf("some common return values:\n"); \
    printf(" UMFPACK_OK=%d\n", UMFPACK_OK); \
    printf(" UMFPACK_WARNING_singular_matrix=%d\n", UMFPACK_WARNING_singular_matrix); \
    printf(" UMFPACK_WARNING_determinant_underflow=%d\n", UMFPACK_WARNING_determinant_underflow); \
    printf(" UMFPACK_WARNING_determinant_overflow=%d\n", UMFPACK_WARNING_determinant_overflow); \
    printf(" UMFPACK_ERROR_out_of_memory=%d\n", UMFPACK_ERROR_out_of_memory); \
    printf(" UMFPACK_ERROR_invalid_Numeric_object=%d\n", UMFPACK_ERROR_invalid_Numeric_object); \
    printf(" UMFPACK_ERROR_invalid_Symbolic_object=%d\n", UMFPACK_ERROR_invalid_Symbolic_object); \
    printf(" UMFPACK_ERROR_invalid_matrix=%d\n", UMFPACK_ERROR_invalid_matrix); \
    printf(" UMFPACK_ERROR_invalid_system=%d\n", UMFPACK_ERROR_invalid_system); \
    printf("for more info do:\n"); \
    printf(" grep UMFPACK_ERROR /usr/include/suitesparse/umfpack.h\n"); \
    fflush(stdout); \
    if(INFO<0) errorexiti("umfpack_dl_solve: dl_solve returned INFO=%d", INFO); }
#else
#define PrintErrorCodesAndExit  errorexit("umfpack is not compiled in")
#endif


/***************************************************************************/
/* some helper routines                                                    */
/***************************************************************************/

/* allocate memory for matrix in umfpack format. Call as:
   allocate_umfpack_dl_matrix(&Ap, &Ai, &Ax, n, nz); */
void allocate_umfpack_dl_matrix(LONGINT **Ap, LONGINT **Ai, double **Ax,
                                LONGINT n, LONGINT nz)
{
  /* allocate memory for matrix */
  *Ap=calloc(n+1, sizeof(LONGINT));
  *Ai=calloc(nz,  sizeof(LONGINT));
  *Ax=calloc(nz,  sizeof(double));
  if(*Ap==NULL || *Ai==NULL || *Ax==NULL)
    errorexit("allocate_umfpack_dl_matrix: out of memory for *Ap, *Ai, *Ax");
}
/* free umfpack matrix allocated by allocate_umfpack_dl_matrix */
void free_umfpack_dl_matrix(LONGINT *Ap, LONGINT *Ai, double *Ax)
{
  free(Ap); free(Ai); free(Ax);
  Ap=NULL;  Ai=NULL;  Ax=NULL;       
}

/* set a matrix in umfpack format (i.e. Ap, Ai, Ax) from Aline */
LONGINT set_umfpack_dl_matrix_from_lines(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                         tSparseVector **Aline, LONGINT nlines,
                                         double dropbelow, int pr)
{
  LONGINT i, j, n;
  double Aij;

  /* set matrix */
  /* All nonzeros are entries, but an entry may be numerically zero. The row
     indices of entries in column j are stored in Ai[Ap[j] ... Ap[j+1]-1].
     The corresponding numerical values are stored in Ax[Ap[j] ... Ap[j+1]-1].
     No duplicate row indices may be present, and the row indices in any
     given column must be sorted in ascending order. The first entry Ap[0]
     must be zero. The total number of entries in the matrix is thus nz =
     Ap[n]. Except for the fact that extra zero entries can be included,
     there is thus a unique compressed column representation of any given
     matrix A. */
  Ap[0] = n = 0;
  for(j = 0; j < nlines; j++)
  {
    Ap[j+1] = Ap[j];
    for(i = 0; i < nlines; i++)
    {
      Aij = GetSparseVectorComponent(Aline[i],j);
      if(fabs(Aij)<=dropbelow) continue;
      Ap[j+1]++;
      Ai[n] = i;
      Ax[n] = Aij;
      n++;
    }
  }
  if(pr)
    printf("set_umfpack_dl_matrix_from_lines: the sparse %ld*%ld matrix "
           "Ap[%ld]=%ld, Ai=%p, Ax=%p is now set!\n",
           nlines, nlines, nlines, Ap[nlines], Ai, Ax);

  if(pr&&0)
  {
    printf("Ax = \n");
    for(i = 0; i < Ap[nlines]; i++)
    {
      printf("%g ", Ax[i]);
    }
    printf("\n");
  }
  /* return number of entries */
  return Ap[nlines];
}

/* set a matrix in umfpack format (i.e. Ap, Ai, Ax) from Acol */
LONGINT set_umfpack_dl_matrix_from_columns(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                           tSparseVector **Acol, LONGINT ncols,
                                           double dropbelow, int pr)
{
  LONGINT i, j, n, ent;
  double Aij;

  /* set matrix */
  /* All nonzeros are entries, but an entry may be numerically zero. The row
     indices of entries in column j are stored in Ai[Ap[j] ... Ap[j+1]-1].
     The corresponding numerical values are stored in Ax[Ap[j] ... Ap[j+1]-1].
     No duplicate row indices may be present, and the row indices in any
     given column must be sorted in ascending order. The first entry Ap[0]
     must be zero. The total number of entries in the matrix is thus nz =
     Ap[n]. Except for the fact that extra zero entries can be included,
     there is thus a unique compressed column representation of any given
     matrix A. */
  Ap[0] = n = 0;
  for(j = 0; j < ncols; j++)
  {
    Ap[j+1] = Ap[j];
    for(ent = 0; ent < Acol[j]->entries; ent++)
    {
      i   = Acol[j]->pos[ent];
      Aij = Acol[j]->val[ent];
      if(fabs(Aij)<=dropbelow) continue;
      Ap[j+1]++;
      Ai[n] = i;
      Ax[n] = Aij;
      n++;
    }
  }
  if(pr)
    printf("set_umfpack_dl_matrix_from_columns: the sparse %ld*%ld matrix "
           "Ap[%ld]=%ld, Ai=%p, Ax=%p is now set!\n",
           ncols, ncols, ncols, Ap[ncols], Ai, Ax);

  if(pr&&0)
  {
    printf("Ax = \n");
    for(i = 0; i < Ap[ncols]; i++)
    {
      printf("%g ", Ax[i]);
    }
    printf("\n");
  }
  /* return number of entries */
  return Ap[ncols];
}


/***************************************************************************/
/* main solver routines                                                    */
/***************************************************************************/

/* solve A x = b with umfpack's umfpack_dl_solve
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int umfpack_dl_solve_fromAcolumns(tSparseVector **Acol,
                                  tVarList *vlx, tVarList *vlb,
                                  double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  LONGINT nz=0;
  double *x;
  double *b;
  LONGINT *Ap;
  LONGINT *Ai;
  double *Ax;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_dl_solve_fromAcolumns: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(j = 0; j < nlines; j++) nz+=Acol[j]->entries;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_dl_solve_fromAcolumns: out of memory for x,b");

  /* allocate memory for matrix */
  allocate_umfpack_dl_matrix(&Ap, &Ai, &Ax, nlines, nz);

  /* set b = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlb->n; j++)
      {
        double *bp = box->v[vlb->index[j]];
        b[line] = bp[i];
        line++;
      }
  }

  /* set matrix */
  set_umfpack_dl_matrix_from_columns(Ap, Ai, Ax, Acol, nlines, dropbelow, pr);
  if(pr) printf("umfpack_dl_solve_fromAcolumns: %ld entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_dl_solve_fromAcolumns: calling umfpack_dl_solve\n"); fflush(stdout); }
  INFO1=umfpack_dl_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_dl_free_symbolic(&Symbolic);
  INFO=umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_dl_free_numeric(&Numeric);
#else
  errorexit("umfpack_dl_solve_fromAcolumns: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_dl_solve_fromAcolumns: umfpack_dl_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_dl_solve_fromAcolumns: setting solution vector\n"); fflush(stdout); }
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *xp = box->v[vlx->index[j]];
        xp[i] = x[line];
        line++;
      }
  }
  if(pr)
  { printf("umfpack_dl_solve_fromAcolumns: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free_umfpack_dl_matrix(Ap, Ai, Ax);

  return INFO;
}


/* solve A x = b with umfpack's umfpack_dl_solve
   for a matrix that is already saved in
   LONGINT *Ap, LONGINT *Ai, double *Ax */
int umfpack_dl_solve_from_Ap_Ai_Ax(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                   tVarList *vlx, tVarList *vlb, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  double *x;
  double *b;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_dl_solve_from_Ap_Ai_Ax: out of memory for x,b");

  /* set b = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlb->n; j++)
      {
        double *bp = box->v[vlb->index[j]];
        b[line] = bp[i];
        line++;
      }
  }

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_dl_solve_from_Ap_Ai_Ax: calling umfpack_dl_solve\n"); fflush(stdout); }
  INFO1=umfpack_dl_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_dl_free_symbolic(&Symbolic);
  INFO=umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_dl_free_numeric(&Numeric);
#else
  errorexit("umfpack_dl_solve_from_Ap_Ai_Ax: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_dl_solve_from_Ap_Ai_Ax: umfpack_dl_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { 
    printf("umfpack_dl_solve_from_Ap_Ai_Ax: setting solution vector\n");
    fflush(stdout);
  }
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *xp = box->v[vlx->index[j]];
        xp[i] = x[line];
        line++;
      }
  }
  if(pr)
  { printf("umfpack_dl_solve_from_Ap_Ai_Ax: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. for x,b */
  free(b);
  free(x);

  return INFO;
}