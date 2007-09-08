/* umfpack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"

#ifdef UMFPACK
#include "umfpack.h"
#endif


/* solve A x = b with umfpack's umfpack_di_solve */
int umfpack_solve(tSparseVector **Aline, tVarList *vlx, tVarList *vlb, int pr)
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
  int INFO;

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
      if(Aij==0.0) continue;
      Ap[j+1]++;
      Ai[n] = i;
      Ax[n] = Aij;
      n++;
    }
  }
  if(pr) printf("umfpack_solve: the sparse %d*%d matrix "
                "Ap[%d]=%d, Ai=%p, Ax=%p is now set!\n",
                nlines, nlines, nlines, Ap[nlines], Ai, Ax);

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
  if(pr) printf("umfpack_solve: calling umfpack_di_solve\n");
  INFO=umfpack_di_symbolic(nlines, nlines, Ap, Ai, Ax, &Symbolic, null, null);
  INFO=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
  umfpack_di_free_symbolic(&Symbolic);
  INFO=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
  umfpack_di_free_numeric(&Numeric);
#else
  errorexit("umfpack_solve: in order to compile with umfpack use MyConfig with\n"
            "DFLAGS += -DUMFPACK\n"
            "SPECIALINCS += -I/usr/include/ufsparse\n"
            "SPECIALLIBS += -lumfpack -lamd");
#endif
  if(pr) printf("umfpack_solve: umfpack_di_solve -> INFO=%d\n", INFO);

  if(INFO!=0)
  {
    printf("INFO=%d\n", INFO);
    errorexiti("umfpack_di_solve: di_solve returned INFO=%d", INFO);
  }

  /* set vlx = x */
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

  /* free mem. */
  free(b);
  free(x);
  free(Ap);
  free(Ai);
  free(Ax);

  return INFO;
}
