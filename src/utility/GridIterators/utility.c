/* utility.c */
/* Wolfgang Tichy 8/2008 */

#include "sgrid.h"
#include "GridIterators.h"


/* dot product of to varlists on grid */
double GridDotProduct(tVarList *vlu, tVarList *vlv)
{
  int i, j, b;
  double sum = 0;
  tGrid *grid = vlu->grid;

  for (j = 0; j < vlu->n; j++)
    forallboxes(grid,b)
    {
      tBox *box = grid->box[b];
      double *u = box->v[vlu->index[j]];
      double *v = box->v[vlv->index[j]];

      forallpoints(box,i)  sum += u[i] * v[i];
    }
  return sum;
}


/* L2 Norm of varlist on grid */
double GridL2Norm(tVarList *vlu)
{
  double sum=0;
  tGrid *grid = vlu->grid;
  int b, i, j, n=0;

  for(j = 0; j < vlu->n; j++)
    forallboxes(grid,b)
    {
      tBox *box = grid->box[b];
      double *u  = box->v[vlu->index[j]];

      forallpoints(box,i) { sum += u[i]*u[i];  n++; }
    }
  return sqrt(sum/n);
}


/* L2 Norm of varlist in box b */
double BoxL2Norm(int b, tVarList *vlu)
{
  double sum=0;
  tGrid *grid = vlu->grid;
  tBox *box = grid->box[b];
  int i, j, n=0;

  for(j = 0; j < vlu->n; j++)
  {
    double *u = box->v[vlu->index[j]];
    forallpoints(box,i) { sum += u[i]*u[i];  n++; }
  }
  return sqrt(sum/n);
}

/* L2 Norm of single variable (with index iu) in box */
double varBoxL2Norm(tBox *box, int iu)
{
  double res;
  tVarList *vlu = vlalloc(box->grid);
  vlpush(vlu, iu);
  res = BoxL2Norm(box->b, vlu);
  vlfree(vlu);
  return res;
}
