/* WTsolver.c */
/* Wolfgang Tichy 12/2007 */

#include "sgrid.h"
#include "GridIterators.h"


/* call all relevant routines to set up an umfpack call using
   umfpack_solve_forSortedVars_fromAcolumns */
int WTsolver(tVarList *vlx, tVarList *vlb, 
             tVarList *r, tVarList *c1,tVarList *c2,
	     int itmax, double tol, double *normres,
	     void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = vlb->grid;
  int bi;
  int pr = Getv("GridIterators_verbose", "yes");
  int line, nlines, nvars, INFO;
  double *x;
  double *b;
  tSparseVector **Aline;

  if(pr) printf("WTsolver:\n");

  /* allocate Aline */
  nvars=vlx->n; 
  nlines = 0;
  forallboxes(grid,bi)  nlines+=(grid->box[bi]->nnodes)*nvars;
  Aline = calloc(nlines, sizeof(*Aline));
  if(Aline)  if(pr) printf("allocated %d matrix lines\n", nlines);
  else       errorexit("no memory for Aline");
  for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();

  /* set Aline */                
  SetMatrixLines_slowly(Aline, lop, r, vlx, c1, c2, pr);
  if(pr&&0) 
    for(line=0; line<nlines; line++) prSparseVector(Aline[line]);

  /* allocate x,b */
  x = calloc(nlines, sizeof(double));
  b = calloc(nlines, sizeof(double));
  if(x==NULL || b==NULL) errorexit("WTsolver: no memory for x,b");
     
  /* set x = vlx , b = vlb */
  copy_varlist_into_array(vlx, x);
  copy_varlist_into_array(vlb, b);

  /* solve A x = b with WTiterator */
  INFO=WTiterator(Aline, nlines, x, b, itmax, tol, normres);
                 
  /* set vlx = x */
  copy_array_into_varlist(x, vlx);

  /* free matrix Aline and x,b */
  for(line=0; line<nlines; line++)  FreeSparseVector(Aline[line]);
  free(Aline);
  free(x);
  free(b);

  return INFO;
}


int WTiterator(tSparseVector **Aline, int nlines,
               double *x, double *b, int itmax, double tol, double *normres)
{
  int i,j,ent, it;
  int pr = Getv("GridIterators_verbose", "yes");
  double *mm;
  double *f;
  double *dx;
  double Aij, newres, sfac,ss;

  mm = calloc(nlines, sizeof(double));
  f  = calloc(nlines, sizeof(double));
  dx = calloc(nlines, sizeof(double));
  if(mm==NULL || f==NULL || dx==NULL)
    errorexit("WTiterator: no memory for mm,f,dx");

  /* mm_j = sum_i A_ij A_ij  = norm2 of column vectors */
  for(i=0; i<nlines; i++)
    for(ent = 0; ent < Aline[j]->entries; ent++)
    {
       j   = Aline[i]->pos[ent];
       Aij = Aline[j]->val[ent];
       mm[j] += Aij*Aij;
    }

  /* start iterations */
  for(it=0; it<itmax; it++)
  { 
    /* compute the residual f_i */
    /* f_i = sum_j A_ij x_j - b_i */
    SparseMatrixLines_times_vector(Aline, nlines, x, f);
    for(i=0; i<nlines; i++)  f[i] -= b[i];
  
    /* dx_i = - sum_j A_ij f_j / mm_j */
    for(i=0; i<nlines; i++)  f[i] = -f[i] / mm[i];
    SparseMatrixLines_times_vector(Aline, nlines, f, dx);

    ss = 1.0;
    sfac = ss;
    do /* do step x_i = x_i + s dx_i and check if it reduces residue */
    {
      if(ss!=1.0 && pr) printf("WTiterator: backtracking to ss=%g\n", ss);

      /* x_i = x_i + s dx_i */
      for(i=0; i<nlines; i++)  x[i] = x[i] + sfac*dx[i];

      /* compute norm of residual f_i */
      SparseMatrixLines_times_vector(Aline, nlines, x, f);
      for(i=0; i<nlines; i++)  f[i] -= b[i];
      newres = scalarproduct_vectors(f,f, nlines);

      ss = ss*0.9;
      sfac = -sfac + ss*0.9;
    } while(newres > *normres);

    *normres = newres;
    if(pr) { printf("WTiterator: %d  %g\n", it, *normres); fflush(stdout); }
    if(*normres <= tol)  break;
  } /* end iterations loop */

  free(mm); 
  free(f);
  return it;
}

/* set x = vlx */
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

/* set vlx = x */
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
