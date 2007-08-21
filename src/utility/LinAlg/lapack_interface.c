/* lapack_interface.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"


/* solve A x = b with lapack's dgesv */
int lapack_dgesv(tSparseVector **Aline, tVarList *vlx, tVarList *vlb)
{
  tGrid *grid = vlx->grid;
  int bi, line;
  int nvars=vlx->n;
  int nlines=0;
  double *xb;
  double *AT;
  int i,j, c1,c2, *pivot, ok;

  /* figure out number of lines */
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;

  /* allocate memory for b, which will be overwritten with x */
  xb=calloc(nlines, sizeof(double));

  /* allocate memory matrix AT = A^T*/
  AT=calloc(nlines*nlines, sizeof(double));

  /* allocate memory for pivot */
  pivot=calloc(nlines, sizeof(int));

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

#ifdef LAPACK
  /* call lapack routine */
  dgesv_(&c1, &c2, AT, &c1, pivot, xb, &c1, &ok);
#else
  errorexit("lapack_dgesv: in order to compile with lapack use MyConfig with\n"
            "DFLAGS += -DLAPACK\n"
            "SPECIALLIBS += -llapack");
#endif
  
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
  free(pivot);

  return ok;
}
