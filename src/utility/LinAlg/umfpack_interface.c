/* umfpack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"

#ifdef UMFPACK
#include "umfpack.h"
#endif


#define PrintErrorCodesAndExit  \
  { printf("umfpack_di_symbolic returned INFO1=%d\n", INFO1); \
    printf("umfpack_di_numeric returned INFO2=%d\n", INFO2);  \
    printf("umfpack_di_solve returned INFO=%d\n", INFO); \
    printf("some common error codes:\n"); \
    printf(" UMFPACK_ERROR_out_of_memory=%d\n", UMFPACK_ERROR_out_of_memory); \
    printf(" UMFPACK_ERROR_invalid_Numeric_object=%d\n", UMFPACK_ERROR_invalid_Numeric_object); \
    printf(" UMFPACK_ERROR_invalid_Symbolic_object=%d\n", UMFPACK_ERROR_invalid_Symbolic_object); \
    printf(" UMFPACK_ERROR_invalid_matrix=%d\n", UMFPACK_ERROR_invalid_matrix); \
    printf(" UMFPACK_ERROR_invalid_system=%d\n", UMFPACK_ERROR_invalid_system); \
    printf("for more info do:\n"); \
    printf(" grep UMFPACK_ERROR /usr/include/suitesparse/umfpack.h\n"); \
    fflush(stdout); \
    errorexiti("umfpack_di_solve: di_solve returned INFO=%d", INFO); }


/***************************************************************************/
/* some helper routines                                                    */
/***************************************************************************/

/* set a matrix in umfpack format (i.e. Ap, Ai, Ax) from Acol */
int set_umfpack_matrix_from_lines(int *Ap, int *Ai, double *Ax,
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
    printf("set_umfpack_matrix_from_lines: the sparse %d*%d matrix "
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
int set_umfpack_matrix_from_columns(int *Ap, int *Ai, double *Ax,
                                    tSparseVector **Acol, int nlines,
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
    printf("set_umfpack_matrix_from_columns: the sparse %d*%d matrix "
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


/***************************************************************************/
/* main solver routines                                                    */
/***************************************************************************/

/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse line vectors that was written 
   by SetMatrixLines_slowly */
int umfpack_solve(tSparseVector **Aline, tVarList *vlx, tVarList *vlb,
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
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_solve: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(i = 0; i < nlines; i++)
    for(j = 0; j < nlines; j++)
      if(GetSparseVectorComponent(Aline[i],j)!=0.0) nz++;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));

  /* allocate memory matrix */
  Ap=calloc(nlines+1, sizeof(double));
  Ai=calloc(nz, sizeof(int));
  Ax=calloc(nz, sizeof(double));
  if(b==NULL || x==NULL || Ap==NULL || Ai==NULL || Ax==NULL)
    errorexit("umfpack_solve: out of memory for x,b, Ap, Ai, Ax");

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
  if(pr) printf("umfpack_solve: the sparse %d*%d matrix "
                "Ap[%d]=%d, Ai=%p, Ax=%p is now set!\n",
                nlines, nlines, nlines, Ap[nlines], Ai, Ax);
  if(pr) printf("umfpack_solve: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  if(pr&&0)
  {
    printf("Ax = \n");
    for(i = 0; i < nz; i++)
    {
      printf("%g ", Ax[i]);
    }
    printf("\n");
  }

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_solve: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
#else
  errorexit("umfpack_solve: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_solve: umfpack_di_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_solve: setting solution vector\n"); fflush(stdout); }
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
  { printf("umfpack_solve: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free(Ap);
  free(Ai);
  free(Ax);

  return INFO;
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse line vectors that was written 
   by SetMatrixLines_forSortedVars_slowly */
int umfpack_solve_forSortedVars(tSparseVector **Aline, 
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
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_solve_forSortedVars: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(i = 0; i < nlines; i++)
    for(j = 0; j < nlines; j++)
      if(GetSparseVectorComponent(Aline[i],j)!=0.0) nz++;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));

  /* allocate memory matrix */
  Ap=calloc(nlines+1, sizeof(double));
  Ai=calloc(nz, sizeof(int));
  Ax=calloc(nz, sizeof(double));
  if(b==NULL || x==NULL || Ap==NULL || Ai==NULL || Ax==NULL)
    errorexit("umfpack_solve_forSortedVars: out of memory for x,b, Ap, Ai, Ax");

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
  if(pr) printf("umfpack_solve_forSortedVars: the sparse %d*%d matrix "
                "Ap[%d]=%d, Ai=%p, Ax=%p is now set!\n",
                nlines, nlines, nlines, Ap[nlines], Ai, Ax);
  if(pr) printf("umfpack_solve_forSortedVars: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

  if(pr&&0)
  {
    printf("Ax = \n");
    for(i = 0; i < nz; i++)
    {
      printf("%g ", Ax[i]);
    }
    printf("\n");
  }

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_solve_forSortedVars: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
#else
  errorexit("umfpack_solve_forSortedVars: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_solve_forSortedVars: umfpack_di_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_solve_forSortedVars: setting solution vector\n"); fflush(stdout); }
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
  if(pr) { printf("umfpack_solve_forSortedVars: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free(Ap);
  free(Ai);
  free(Ax);

  return INFO;
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int umfpack_solve_fromAcolumns(tSparseVector **Acol,
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
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_solve_fromAcolumns: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(j = 0; j < nlines; j++) nz+=Acol[j]->entries;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));

  /* allocate memory matrix */
  Ap=calloc(nlines+1, sizeof(double));
  Ai=calloc(nz, sizeof(int));
  Ax=calloc(nz, sizeof(double));
  if(b==NULL || x==NULL || Ap==NULL || Ai==NULL || Ax==NULL)
    errorexit("umfpack_solve_fromAcolumns: out of memory for x,b, Ap, Ai, Ax");

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
  set_umfpack_matrix_from_columns(Ap, Ai, Ax, Acol, nlines, dropbelow, pr);
  if(pr) printf("umfpack_solve_fromAcolumns: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_solve_fromAcolumns: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
#else
  errorexit("umfpack_solve_fromAcolumns: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_solve_fromAcolumns: umfpack_di_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_solve_fromAcolumns: setting solution vector\n"); fflush(stdout); }
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
  { printf("umfpack_solve_fromAcolumns: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free(Ap);
  free(Ai);
  free(Ax);

  return INFO;
}


/* solve A x = b with umfpack's umfpack_di_solve
   for a matrix made up of sparse column vectors that was written 
   by SetMatrixColumns_slowly */
int umfpack_solve_forSortedVars_fromAcolumns(tSparseVector **Acol,
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
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int INFO, INFO1, INFO2;

  if(pr) { printf("umfpack_solve_forSortedVars_fromAcolumns: setting sparse matrix\n"); fflush(stdout); }

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* count number of entries in sparse matrix */
  for(j = 0; j < nlines; j++) nz+=Acol[j]->entries;

  /* allocate memory for b and x */
  b=calloc(nlines, sizeof(double));
  x=calloc(nlines, sizeof(double));

  /* allocate memory matrix */
  Ap=calloc(nlines+1, sizeof(double));
  Ai=calloc(nz, sizeof(int));
  Ax=calloc(nz, sizeof(double));
  if(b==NULL || x==NULL || Ap==NULL || Ai==NULL || Ax==NULL)
    errorexit("umfpack_solve_forSortedVars_fromAcolumns: out of memory for x,b, Ap, Ai, Ax");

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
  set_umfpack_matrix_from_columns(Ap, Ai, Ax, Acol, nlines, dropbelow, pr);
  if(pr) printf("umfpack_solve_forSortedVars_fromAcolumns: %d entries of magnitude <= %g were dropped\n",
                nz-Ap[nlines], dropbelow);

#ifdef UMFPACK
  /* call umfpack routine */
  if(pr)
  { printf("umfpack_solve_forSortedVars_fromAcolumns: calling umfpack_di_solve\n"); fflush(stdout); }
  INFO1=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO2=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
#else
  errorexit("umfpack_solve_forSortedVars_fromAcolumns: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/suitesparse\n"
            "SPECIALLIBS += -lumfpack -lamd -lblas");
#endif
  if(pr)
  { 
    printf("umfpack_solve_forSortedVars_fromAcolumns: umfpack_di_solve -> INFO=%d\n", INFO); 
    fflush(stdout);
  }

  if(INFO!=0)
    PrintErrorCodesAndExit;

  /* set vlx = x */
  if(pr)
  { printf("umfpack_solve_forSortedVars_fromAcolumns: setting solution vector\n"); fflush(stdout); }
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
  { printf("umfpack_solve_forSortedVars_fromAcolumns: vector vlx=%p is now set!\n", vlx); fflush(stdout);}

  /* free mem. */
  free(b);
  free(x);
  free(Ap);
  free(Ai);
  free(Ax);

  return INFO;
}
