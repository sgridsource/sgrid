/* utility.c */
/* Wolfgang Tichy 8/2008 */

#include "sgrid.h"
#include "GridIterators.h"



double dot(tVarList *vlu, tVarList *vlv)
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



double norm2(tVarList *vlu)
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
