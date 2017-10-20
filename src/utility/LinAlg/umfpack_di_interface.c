/* umfpack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"

#ifdef UMFPACK
#include "umfpack.h"
#define PrintErrorCodes  \
  { printf("umfpack_di_symbolic returned INFO1=%d\n", INFO1); \
    printf("umfpack_di_numeric returned INFO2=%d\n", INFO2);  \
    printf("umfpack_di_solve returned INFO=%d\n", INFO); \
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
    fflush(stdout); }
#else
#define PrintErrorCodes  errorexit("umfpack is not compiled in")
#endif

#define PrintErrorCodesAndExit  \
  { PrintErrorCodes; \
    if(INFO<0) errorexiti("umfpack_di_solve: di_solve returned INFO=%d", INFO); }

/***************************************************************************/
/* some helper routines                                                    */
/***************************************************************************/

/* allocate memory for matrix in umfpack format. Call as:
   allocate_umfpack_di_matrix(&Ap, &Ai, &Ax, n, nz); */
void allocate_umfpack_di_matrix(int **Ap, int **Ai, double **Ax, int n, int nz)
{
  /* allocate memory for matrix */
  *Ap=calloc(n+1, sizeof(int));
  *Ai=calloc(nz,  sizeof(int));
  *Ax=calloc(nz,  sizeof(double));
  if(*Ap==NULL || *Ai==NULL || *Ax==NULL)
    errorexit("allocate_umfpack_di_matrix: out of memory for *Ap, *Ai, *Ax");
}
/* free umfpack matrix allocated by allocate_umfpack_di_matrix */
void free_umfpack_di_matrix(int *Ap, int *Ai, double *Ax)
{
  free(Ap); free(Ai); free(Ax);
  Ap=NULL;  Ai=NULL;  Ax=NULL;       
}

/* set a matrix in umfpack format (i.e. Ap, Ai, Ax) from Acol */
int set_umfpack_di_matrix_from_lines(int *Ap, int *Ai, double *Ax,
                                     tSparseVector **Aline, int nlines,
                                     double dropbelow, int pr)
{
  int i, j, n, ent;
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
    printf("set_umfpack_di_matrix_from_lines: the sparse %d*%d matrix "
           "Ap[%d]=%d, Ai=%p, Ax=%p is now set!\n",
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
int set_umfpack_di_matrix_from_columns(int *Ap, int *Ai, double *Ax,
                                       tSparseVector **Acol, int ncols,
                                       double dropbelow, int pr)
{
  int i, j, n, ent;
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
    printf("set_umfpack_di_matrix_from_columns: the sparse %d*%d matrix "
           "Ap[%d]=%d, Ai=%p, Ax=%p is now set!\n",
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

/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse line vectors that was written 
   by SetMatrixLines_slowly */
int umfpack_di_solve_fromAlines(tSparseVector **Aline, tVarList *vlx,
                                tVarList *vlb, double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j,n;
  double Aij;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  int nz=0;
  double *x;
  double *b;
  int *Ap;
  int *Ai;
  double *Ax;
#ifdef UMFPACK
  double Info[UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
#endif
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_di_solve_fromAlines: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(i = 0; i < nlines; i++)
    for(j = 0; j < nlines; j++)
      if(GetSparseVectorComponent(Aline[i],j)!=0.0) nz++;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_di_solve_fromAlines: out of memory for x,b");

  /* allocate memory for matrix */
  allocate_umfpack_di_matrix(&Ap, &Ai, &Ax, nlines, nz);

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
  set_umfpack_di_matrix_from_lines(Ap, Ai, Ax, Aline, nlines, dropbelow, pr);
  if(pr) printf("umfpack_di_solve_fromAlines: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  INFO1=INFO2=INFO=-6662442;
#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_di_solve_fromAlines: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, Info);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
  if(pr)
  {
    printf("  Info[UMFPACK_RCOND]=%g\n", Info[UMFPACK_RCOND]);
    printf("umfpack_di_solve_fromAlines: umfpack_di_solve -> INFO=%d\n", INFO);
    fflush(stdout);
  }
#else
  errorexit("umfpack_di_solve_fromAlines: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif

  if(INFO!=0) PrintErrorCodesAndExit;
  if(INFO1!=0 || INFO2!=0) PrintErrorCodes;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_di_solve_fromAlines: setting solution vector\n"); fflush(stdout); }
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
  { printf("umfpack_di_solve_fromAlines: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free_umfpack_di_matrix(Ap, Ai, Ax);

  return translate_SuiteSparse_warnings(INFO);
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse line vectors that was written 
   by SetMatrixLines_forSortedVars_slowly */
int umfpack_di_solve_forSortedVars_fromAlines(tSparseVector **Aline, 
				              tVarList *vlx, tVarList *vlb,
                  		              double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j,n;
  double Aij;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  int nz=0;
  double *x;
  double *b;
  int *Ap;
  int *Ai;
  double *Ax;
#ifdef UMFPACK
  double Info[UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
#endif
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_di_solve_forSortedVars_fromAlines: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(i = 0; i < nlines; i++)
    for(j = 0; j < nlines; j++)
      if(GetSparseVectorComponent(Aline[i],j)!=0.0) nz++;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_di_solve_forSortedVars_fromAlines: out of memory for x,b");

  /* allocate memory for matrix */
  allocate_umfpack_di_matrix(&Ap, &Ai, &Ax, nlines, nz);

  /* set b = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlb->n; j++)
    {
      double *bp = box->v[vlb->index[j]];
      forallpoints(box,i)
      {
        b[line] = bp[i];
        line++;
      }
    }
  }

  /* set matrix */
  set_umfpack_di_matrix_from_lines(Ap, Ai, Ax, Aline, nlines, dropbelow, pr);
  if(pr) printf("umfpack_di_solve_forSortedVars_fromAlines: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  INFO1=INFO2=INFO=-6662442;
#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_di_solve_forSortedVars_fromAlines: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, Info);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
  if(pr)
  {
    printf("  Info[UMFPACK_RCOND]=%g\n", Info[UMFPACK_RCOND]);
    printf("umfpack_di_solve_forSortedVars_fromAlines: umfpack_di_solve -> INFO=%d\n", INFO);
    fflush(stdout);
  }
#else
  errorexit("umfpack_di_solve_forSortedVars_fromAlines: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif

  if(INFO!=0) PrintErrorCodesAndExit;
  if(INFO1!=0 || INFO2!=0) PrintErrorCodes;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_di_solve_forSortedVars_fromAlines: setting solution vector\n"); fflush(stdout); }
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlx->n; j++)
    {
      double *xp = box->v[vlx->index[j]];
      forallpoints(box,i)
      {
        xp[i] = x[line];
        line++;
      }
    } 
  } 
  if(pr) { printf("umfpack_di_solve_forSortedVars_fromAlines: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free_umfpack_di_matrix(Ap, Ai, Ax);

  return translate_SuiteSparse_warnings(INFO);
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int umfpack_di_solve_fromAcolumns(tSparseVector **Acol,
                                  tVarList *vlx, tVarList *vlb,
                                  double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  int nz=0;
  double *x;
  double *b;
  int *Ap;
  int *Ai;
  double *Ax;
#ifdef UMFPACK
  double Info[UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
#endif
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_di_solve_fromAcolumns: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(j = 0; j < nlines; j++) nz+=Acol[j]->entries;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_di_solve_fromAcolumns: out of memory for x,b");

  /* allocate memory for matrix */
  allocate_umfpack_di_matrix(&Ap, &Ai, &Ax, nlines, nz);

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
  set_umfpack_di_matrix_from_columns(Ap, Ai, Ax, Acol, nlines, dropbelow, pr);
  if(pr) printf("umfpack_di_solve_fromAcolumns: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  INFO1=INFO2=INFO=-6662442;
#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_di_solve_fromAcolumns: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, Info);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
  if(pr)
  {
    printf("  Info[UMFPACK_RCOND]=%g\n", Info[UMFPACK_RCOND]);
    printf("umfpack_di_solve_fromAcolumns: umfpack_di_solve -> INFO=%d\n", INFO);
    fflush(stdout);
  }
#else
  errorexit("umfpack_di_solve_fromAcolumns: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif

  if(INFO!=0) PrintErrorCodesAndExit;
  if(INFO1!=0 || INFO2!=0) PrintErrorCodes;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_di_solve_fromAcolumns: setting solution vector\n"); fflush(stdout); }
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
  { printf("umfpack_di_solve_fromAcolumns: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free_umfpack_di_matrix(Ap, Ai, Ax);

  return translate_SuiteSparse_warnings(INFO);
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int umfpack_di_solve_forSortedVars_fromAcolumns(tSparseVector **Acol,
      tVarList *vlx, tVarList *vlb,
      double dropbelow, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  int nz=0;
  double *x;
  double *b;
  int *Ap;
  int *Ai;
  double *Ax;
#ifdef UMFPACK
  double Info[UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
#endif
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_di_solve_forSortedVars_fromAcolumns: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(j = 0; j < nlines; j++) nz+=Acol[j]->entries;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_di_solve_forSortedVars_fromAcolumns: out of memory for x,b");

  /* allocate memory for matrix */
  allocate_umfpack_di_matrix(&Ap, &Ai, &Ax, nlines, nz);

  /* set b = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlb->n; j++)
    {
      double *bp = box->v[vlb->index[j]];
      forallpoints(box,i)
      {
        b[line] = bp[i];
        line++;
      }
    }
  }

  /* set matrix */
  set_umfpack_di_matrix_from_columns(Ap, Ai, Ax, Acol, nlines, dropbelow, pr);
  if(pr) printf("umfpack_di_solve_forSortedVars_fromAcolumns: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  INFO1=INFO2=INFO=-6662442;
#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_di_solve_forSortedVars_fromAcolumns: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, Info);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
  if(pr)
  {
    printf("  Info[UMFPACK_RCOND]=%g\n", Info[UMFPACK_RCOND]);
    printf("umfpack_di_solve_forSortedVars_fromAcolumns: umfpack_di_solve -> INFO=%d\n", INFO);
    fflush(stdout);
  }
#else
  errorexit("umfpack_di_solve_forSortedVars_fromAcolumns: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif

  if(INFO!=0) PrintErrorCodesAndExit;
  if(INFO1!=0 || INFO2!=0) PrintErrorCodes;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_di_solve_forSortedVars_fromAcolumns: setting solution vector\n"); fflush(stdout); }
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlx->n; j++)
    {
      double *xp = box->v[vlx->index[j]];
      forallpoints(box,i)
      {
        xp[i] = x[line];
        line++;
      }
    }
  }
  if(pr)
  { printf("umfpack_di_solve_forSortedVars_fromAcolumns: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free_umfpack_di_matrix(Ap, Ai, Ax);

  return translate_SuiteSparse_warnings(INFO);
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix that is already saved in
   int *Ap, int *Ai, double *Ax */
int umfpack_di_solve_from_Ap_Ai_Ax(int *Ap, int *Ai, double *Ax,
                                   tVarList *vlx, tVarList *vlb, int pr)
{
  tGrid *grid = vlx->grid;
  int i,j;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  double *x;
  double *b;
#ifdef UMFPACK
  double Info[UMFPACK_INFO];
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
#endif
  int INFO, INFO1, INFO2;

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));
  if(b==NULL || x==NULL)
    errorexit("umfpack_di_solve_from_Ap_Ai_Ax: out of memory for x,b");

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

  INFO1=INFO2=INFO=-6662442;
#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_di_solve_from_Ap_Ai_Ax: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, Info);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
  if(pr)
  {
    printf("  Info[UMFPACK_RCOND]=%g\n", Info[UMFPACK_RCOND]);
    printf("umfpack_di_solve_from_Ap_Ai_Ax: umfpack_di_solve -> INFO=%d\n", INFO);
    fflush(stdout);
  }
#else
  errorexit("umfpack_di_solve_from_Ap_Ai_Ax: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif

  if(INFO!=0) PrintErrorCodesAndExit;
  if(INFO1!=0 || INFO2!=0) PrintErrorCodes;

  /* set vlx = x */
  if(pr)
  { 
    printf("umfpack_di_solve_from_Ap_Ai_Ax: setting solution vector\n");
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
  { printf("umfpack_di_solve_from_Ap_Ai_Ax: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. for x,b */
  free(b);
  free(x);

  return translate_SuiteSparse_warnings(INFO);
}
