/* lapack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"


/* solve A x = b with lapack's dgesv */
int lapack_dgesv(tSparseVector **Aline, tVarList *vlx, tVarList *vlb, int pr)
{
  tGrid *grid = vlx->grid;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  double *xb;
  double *AT;
  int i,j, *IPIV, INFO, NRHS=1;

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* allocate memory for b, which will be overwritten with x */
  xb=calloc(nlines, sizeof(double));

  /* allocate memory matrix AT = A^T*/
  AT=calloc(nlines*nlines, sizeof(double));

  /* allocate memory for IPIV */
  IPIV=calloc(nlines, sizeof(int));

  /* set xb = vlb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlb->n; j++)
      {
        double *b = box->v[vlb->index[j]];
        xb[line] = b[i];
        line++;
      }
  }

  /* set AT = A^T :  AT[j+nlines*i]=A_{ji}; */
  for(i = 0; i < nlines; i++)
    for(j = 0; j < nlines; j++)
      AT[j+nlines*i] = GetSparseVectorComponent(Aline[j],i);

  if(pr) printf("lapack_dgesv: the %d*%d matrix AT=%p is now set!\n",
                nlines, nlines, AT);

  if(pr&&0)
  {
    printf("AT[i][j]:\n");
    for(i = 0; i < nlines; i++)
    {
      printf("AT[%d][j] = ", i);
      for(j = 0; j < nlines; j++)  printf("%g ", AT[j+nlines*i]);
      printf("\n");
    }
  }

#ifdef LAPACK
  /* call lapack routine */
  if(pr) printf("lapack_dgesv: calling lapack's dgesv\n");
  dgesv_(&nlines, &NRHS, AT, &nlines, IPIV, xb, &nlines, &INFO);
#else
  errorexit("lapack_dgesv: in order to compile with lapack use MyConfig with\n"
            "DFLAGS += -DLAPACK\n"
            "SPECIALLIBS += -llapack");
#endif
  if(pr) printf("lapack_dgesv: dgesv -> INFO=%d\n", INFO);

  if(INFO!=0)
  {
    printf("INFO=%d\n", INFO);
    printf("  if INFO<0: INFO =-i, the i-th argument had an illegal value\n");
    printf("  if INFO>0: INFO = i, U(i,i) is exactly zero.\n"
           "             The factorization has been completed, "
           "but the factor U is\n"
           "             exactly singular, "
           "so the solution could not be computed.\n");
    errorexiti("lapack_dgesv: dgesv returned INFO=%d", INFO);
  }

  /* set vlx = xb */
  line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *x = box->v[vlx->index[j]];
        x[i] = xb[line];
        line++;
      }
  }

  free(xb);
  free(AT);
  free(IPIV);

  return INFO;
}
