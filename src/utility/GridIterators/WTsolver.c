/* WTsolver.c */
/* Wolfgang Tichy 12/2007 */

#include "sgrid.h"
#include "GridIterators.h"


/* prepare matrix and vectors and then call my own iterator */
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
  if(Aline) { if(pr) printf("allocated %d matrix lines\n", nlines); }
  else        errorexit("no memory for Aline");
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
  int i, it;
  int pr = Getv("GridIterators_verbose", "yes");
  double *mm;
  double *f;
  double *dx;
  double mag_dx, newres, sfac,s0;

  mm = calloc(nlines, sizeof(double));
  f  = calloc(nlines, sizeof(double));
  dx = calloc(nlines, sizeof(double));
  if(mm==NULL || f==NULL || dx==NULL)
    errorexit("WTiterator: no memory for mm,f,dx");

//nlines=2;
//Aline[0]->entries=2;
//Aline[1]->entries=2;
//
//Aline[0]->pos[0]=0;
//Aline[0]->pos[1]=1;
//Aline[1]->pos[0]=0;
//Aline[1]->pos[1]=1;
//
//Aline[0]->val[0]=1;
//Aline[0]->val[1]=2;
//Aline[1]->val[0]=1;
//Aline[1]->val[1]=-1;
//
//b[0]=3;
//b[1]=0;
//
//x[0]=-22;
//x[1]=41;

  /* compute the residual f_i */
  /* f_i = sum_j A_ij x_j - b_i */
  SparseMatrixLines_times_vector(Aline, nlines, x, f);
  for(i=0; i<nlines; i++)  f[i] -= b[i];
  *normres = sqrt(scalarproduct_vectors(f,f, nlines)/nlines);
  if(pr) 
  { 
    printf("WTiterator: %d  res=%g\n", 0, *normres); fflush(stdout);
  }
  if(*normres <= tol) return 0;

  /* start iterations */
  for(it=1; it<=itmax; it++)
  { 
    /* S := (x^T A^T - b^T)(A x - b)/2
       dS/dx_j = [A^T A x]_j - [A^T b]_j */
    /* compute dx_j = - [A^T A x]_j + [A^T b]_j */
    SparseMatrixLines_times_vector(Aline, nlines, x, mm);
    SparseMatrixLinesTranspose_times_vector(Aline, nlines, mm, f);
    for(i=0; i<nlines; i++)  dx[i] = -f[i];
    SparseMatrixLinesTranspose_times_vector(Aline, nlines, b, mm);
    for(i=0; i<nlines; i++)  dx[i] += mm[i];
    mag_dx = sqrt(scalarproduct_vectors(dx,dx, nlines));
    if(mag_dx!=0.0) for(i=0; i<nlines; i++) dx[i] = *normres * dx[i]/mag_dx;
//printf("WTiterator: %d  %g   ", it, *normres);
//printf("(x,y)=(%g,%g) (dx,dy)=(%g,%g)\n",x[0],x[1], dx[0],dx[1]);
    
    sfac = s0 = 1.0;
    do /* do step x_i = x_i + s dx_i and check if it reduces residue */
    {
      /* x_i = x_i + s dx_i */
      for(i=0; i<nlines; i++)  x[i] = x[i] + sfac*dx[i];

      /* compute norm of residual f_i */
      SparseMatrixLines_times_vector(Aline, nlines, x, f);
      for(i=0; i<nlines; i++)  f[i] -= b[i];
      newres = sqrt(scalarproduct_vectors(f,f, nlines)/nlines);

      if(newres <= *normres) break; 
      
      /* step was bad, go back and reduce sfac */
      for(i=0; i<nlines; i++)  x[i] = x[i] - sfac*dx[i];
      sfac = sfac*0.1;
      if(pr&&0) printf("WTiterator: newres=%g backtracking to sfac=%g\n", 
                    newres, sfac);
    } while(sfac>0);

    *normres = newres;
    if(pr)
    { 
      printf("WTiterator: %d  res=%g  sfac=%g\n", it, *normres, sfac);
      fflush(stdout);
    }

//printf("WTiterator: %d  %g   ", it, *normres);
//printf("(x,y)=(%g,%g) (dx,dy)=(%g,%g)\n",x[0],x[1], dx[0],dx[1]);
    if(*normres <= tol)  break;
  } /* end iterations loop */

  free(mm); 
  free(f);
  free(dx);
  return it;
}


/* this one works only if nlines = 2, bad!!! */
int WTiterator1(tSparseVector **Aline, int nlines,
               double *x, double *b, int itmax, double tol, double *normres)
{
  int i,ent, it;
  int pr = Getv("GridIterators_verbose", "yes");
  double *mm;
  double *f;
  double *dx;
  double Aij, newres, sfac,s0;

  mm = calloc(nlines, sizeof(double));
  f  = calloc(nlines, sizeof(double));
  dx = calloc(nlines, sizeof(double));
  if(mm==NULL || f==NULL || dx==NULL)
    errorexit("WTiterator: no memory for mm,f,dx");
//
//nlines=2;
//Aline[0]->entries=2;
//Aline[1]->entries=2;
//
//Aline[0]->pos[0]=0;
//Aline[0]->pos[1]=1;
//Aline[1]->pos[0]=0;
//Aline[1]->pos[1]=1;
//
//Aline[0]->val[0]=1;
//Aline[0]->val[1]=2;
//Aline[1]->val[0]=1;
//Aline[1]->val[1]=-1;
//
//b[0]=3333;
//b[1]=0;
//
//x[0]=-22;
//x[1]=41;

  /* mm_i = sum_i A_ij A_ij  = norm2 of row vectors */
  for(i=0; i<nlines; i++)
    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       /* j   = Aline[i]->pos[ent]; */
       Aij = Aline[i]->val[ent];
       mm[i] += Aij*Aij;
    }

  /* compute the residual f_i */
  /* f_i = sum_j A_ij x_j - b_i */
  SparseMatrixLines_times_vector(Aline, nlines, x, f);
  for(i=0; i<nlines; i++)  f[i] -= b[i];
  *normres = sqrt(scalarproduct_vectors(f,f, nlines)/nlines);
  if(pr) { printf("WTiterator: %d  %g\n", 0, *normres); fflush(stdout); }

  /* start iterations */
  for(it=1; it<=itmax; it++)
  { 
    /* dx_j = - sum_k f_k A_kj / mm_k */
    for(i=0; i<nlines; i++)  f[i] = -f[i] / mm[i];
    vector_times_SparseMatrixLines(f, Aline, nlines, dx);

    sfac = s0 = 1.0;
    do /* do step x_i = x_i + s dx_i and check if it reduces residue */
    {
      /* x_i = x_i + s dx_i */
      for(i=0; i<nlines; i++)  x[i] = x[i] + sfac*dx[i];

      /* compute norm of residual f_i */
      SparseMatrixLines_times_vector(Aline, nlines, x, f);
      for(i=0; i<nlines; i++)  f[i] -= b[i];
      newres = sqrt(scalarproduct_vectors(f,f, nlines)/nlines);
//printf("%d  %g      ", it-1, *normres);
//break;
      if(newres <= *normres) break; 
      
      /* step was bad, go back and reduce sfac */
      for(i=0; i<nlines; i++)  x[i] = x[i] - sfac*dx[i];
      sfac = sfac*0.9;
      if(sfac>0 && sfac<1e-6) sfac = -s0;
      if(pr) printf("WTiterator: newres=%g backtracking to sfac=%g\n", 
                    newres, sfac);
    } while(sfac>0 || sfac<-1e-6);

    *normres = newres;
    if(pr) { printf("WTiterator: %d  %g\n", it, *normres); fflush(stdout); }
//printf("(x,y)=(%g,%g) (dx,dy)=(%g,%g)\n",x[0],x[1], dx[0],dx[1]);

    if(*normres <= tol)  break;
  } /* end iterations loop */

  free(mm); 
  free(f);
  free(dx);
  return it;
}
