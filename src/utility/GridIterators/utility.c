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
    int i;

    forallpoints(box,i)
    {
      int j;
      for(j = 0; j < vlx->n; j++)
      {
        double *px = box->v[vlx->index[j]];
        x[line] = px[i];
        line++;
      }
    }
  }
}
/* same as above but with thread safe: int line */
void copy_varlist_into_array_LEVEL3(tVarList *vlx, double *x)
{
  tGrid *grid = vlx->grid;
  int nvar = vlx->n;
  int bi;
  
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;
    int lb = nvar*nnodes*bi;
    int i;

    SGRID_LEVEL3_Pragma(omp parallel for)
    forallpoints(box,i)
    {
      int li = nvar*i + lb;
      int j;
      for(j = 0; j < nvar; j++)
      {
        double *px = box->v[vlx->index[j]];
        int line = j + li;
        x[line] = px[i];
      }
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
    int i;

    forallpoints(box,i)
    {
      int j;
      for(j = 0; j < vlx->n; j++)
      {
        double *px = box->v[vlx->index[j]];
        px[i] = x[line];
        line++;
      }
    }
  }
}
/* same as above but with thread safe: int line */
void copy_array_into_varlist_LEVEL3(double *x, tVarList *vlx)
{
  tGrid *grid = vlx->grid;
  int nvar = vlx->n;
  int bi;
  
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;
    int lb = nvar*nnodes*bi;
    int i;

    SGRID_LEVEL3_Pragma(omp parallel for)
    forallpoints(box,i)
    {
      int li = nvar*i + lb;
      int j;
      for(j = 0; j < nvar; j++)
      {
        double *px = box->v[vlx->index[j]];
        int line = j + li;
        px[i] = x[line];
      }
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
void copy_array_into_vl_outervlloop(double *xa, tVarList *vlx)
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

/* copy varlistComp In Subbox_into_an array: x = vlx of comp in subbox */
void copy_varlistCompInSubbox_into_array(tVarList *vlx, int vi, int bi,
                                         int sbi, int sbj, int sbk,
                                         int nsb1, int nsb2, int nsb3,
                                         double *x)
{
  tGrid *grid = vlx->grid;
  tBox *box = grid->box[bi];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i1,i2, j1,j2, k1,k2; 
  int i,j,k;
        
  IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);

  for(k=k1; k<k2; k++)
  for(j=j1; j<j2; j++)
  for(i=i1; i<i2; i++)
  {
    double *xp = box->v[vlx->index[vi]];
    int ijk = Index(i,j,k); 
    int ix  = Ind_n1n2( (i-i1),(j-j1),(k-k1), (i2-i1),(j2-j1) );

    x[ix] = xp[ijk];
  }
}

/* copy array into varlistComp In Subbox: vlx of comp in subbox = x */
void copy_array_into_varlistCompInSubbox(double *x, 
                                         tVarList *vlx, int vi, int bi,
                                         int sbi, int sbj, int sbk,
                                         int nsb1, int nsb2, int nsb3)
{
  tGrid *grid = vlx->grid;
  tBox *box = grid->box[bi];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i1,i2, j1,j2, k1,k2; 
  int i,j,k;
        
  IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);

  for(k=k1; k<k2; k++)
  for(j=j1; j<j2; j++)
  for(i=i1; i<i2; i++)
  {
    double *xp = box->v[vlx->index[vi]];
    int ijk = Index(i,j,k); 
    int ix  = Ind_n1n2( (i-i1),(j-j1),(k-k1), (i2-i1),(j2-j1) );

    xp[ijk] = x[ix];
  }
}
