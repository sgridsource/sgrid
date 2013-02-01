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
  vlpushone(vlu, iu);
  res = BoxL2Norm(box->b, vlu);
  vlfree(vlu);
  return res;
}


/*******************************************************/
/* copy between ordinary C-arrays and varlists on grid */
/******************************************************/
/* set x = vlx compatible with SetMatrixLines_slowly */
void copy_varlist_into_array(tVarList *vlx, double *x)
{
  tGrid *grid = vlx->grid;
  int bi;
  int line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *px = box->v[vlx->index[j]];
        x[line] = px[i];
        line++;
      }
  }
}

/* set vlx = x compatible with SetMatrixLines_slowly */
void copy_array_into_varlist(double *x, tVarList *vlx)
{
  tGrid *grid = vlx->grid;
  int bi;
  int line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *px = box->v[vlx->index[j]];
        px[i] = x[line];
        line++;
      }
  }
}

/* set x = vlx compatible with SetMatrixLines_forSortedVars_slowly */
void copy_varlist_into_array_forSortedVars(tVarList *vlx, double *x)
{
  tGrid *grid = vlx->grid;
  int bi;
  int line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlx->n; j++)
      forallpoints(box,i)
      {
        double *px = box->v[vlx->index[j]];
        x[line] = px[i];
        line++;
      }
  }
}

/* set vlx = x compatible with SetMatrixLines_forSortedVars_slowly */
void copy_array_into_varlist_forSortedVars(double *x, tVarList *vlx)
{
  tGrid *grid = vlx->grid;
  int bi;
  int line = 0;
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int i,j;

    for(j = 0; j < vlx->n; j++)
      forallpoints(box,i)
      {
        double *px = box->v[vlx->index[j]];
        px[i] = x[line];
        line++;
      }
  }
}

/* copy var list vlx into array xa. The outer loop is over vars in vlx. */
/* xa = vlx */
void copy_vl_into_array_outervlloop(tVarList *vlx, double *xa)
{
  tGrid *grid = vlx->grid;
  int j, i, boxi, n;
  
  for(n=0,j=0; j<vlx->n; j++)
    forallboxes(grid,boxi)
    {
      tBox *box = grid->box[boxi];      
      double *px  = vlldataptr(vlx,  box, j);
      forallpoints(box,i) 
      {
        xa[n] = px[i];
        n++;
      }
    }
}

/* copy array xa into var list vlx. The outer loop is over vars in vlx. */
/* vlx = xa */
void copy_array_into_vl_outervlloop(tVarList *vlx, double *xa)
{
  tGrid *grid = vlx->grid;
  int j, i, boxi, n;
  
  for(n=0,j=0; j<vlx->n; j++)
    forallboxes(grid,boxi)
    {
      tBox *box = grid->box[boxi];      
      double *px  = vlldataptr(vlx,  box, j);
      forallpoints(box,i) 
      {
        px[i] = xa[n];
        n++;
      }
    }
}
