/* wrappers_for_templates.c */
/* Wolfgang Tichy 8/2007 */

/* from "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates/				*/

#include "sgrid.h"
#include "GridIterators.h"
#include "wrappers_for_templates.h"

/* global vars in this file */
void (*lop_fortemplates)(tVarList *, tVarList *, tVarList *, tVarList *);
void (*precon_fortemplates)(tVarList *, tVarList *, tVarList *, tVarList *);
tVarList *r_fortemplates;
tVarList *x_fortemplates;
tVarList *c1_fortemplates;
tVarList *c2_fortemplates;
int dim_fortemplates;


/* copy var list vlx into array xa. The outer loop is over vars in vlx. */
/* xa = vlx */
int copy_vl_into_array_outervlloop(tVarList *vlx, double *xa, int dim)
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
        if(n>=dim) return dim;
      }
    }
  return n;
}

/* copy array xa into var list vlx. The outer loop is over vars in vlx. */
/* vlx = xa */
int copy_array_into_vl_outervlloop(tVarList *vlx, double *xa, int dim)
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
        if(n>=dim) return dim;
      }
    }
  return n;
}


/* y := alpha*A*x + beta*y, where A is matrix defined by lop */
int matvec(double *alpha, double *x, double *beta, double *y)
{
  /* compute r = A*x */
  copy_array_into_vl_outervlloop(x_fortemplates, x, dim_fortemplates);
  lop_fortemplates(r_fortemplates, x_fortemplates, c1_fortemplates, c2_fortemplates);

  /* set x_fortemplates = y */
  copy_array_into_vl_outervlloop(x_fortemplates, y, dim_fortemplates);

  /* set r=alpha*r + beta*x_fortemplates */
  vladd(r_fortemplates, *alpha, r_fortemplates, *beta, x_fortemplates);

  /* copy r into y */
  copy_vl_into_array_outervlloop(r_fortemplates, y, dim_fortemplates);
  
  return 0;
}

/* Precon: solves M*x = b for x */
int psolve(double *x, double *b)
{
  /* solve M*x = b for x, solution is x_fortemplates, b is in r_fortemplates */
  copy_array_into_vl_outervlloop(r_fortemplates, b, dim_fortemplates);
  precon_fortemplates(x_fortemplates, r_fortemplates, c1_fortemplates, c2_fortemplates);

  /* copy x_fortemplates into x */
  copy_vl_into_array_outervlloop(x_fortemplates, x, dim_fortemplates);
  
  return 0;
}


/* call GMRES from templates */
int templates_gmres_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  int N; /* dim of matrix */
  double *B;
  double *X;
  int RESTRT;
  double *WORK;		int LDW;
  double *H;		int LDH;
  int ITER;
  double RESID;
  int INFO;

  /* set int vars */
  N = 0 ;
  for(j = 0; j < b->n; j++)
    forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  RESTRT = 1000;
  LDW = (N) + 1; 
  LDH = (RESTRT+1) + 1;
  ITER = itmax;
  RESID = tol;
  if(pr) printf("templates_gmres_wrapper: itmax=%d tol=%g "
                "N=%d RESTRT=%d LDW=%d LDH=%d\n",
                itmax, tol, N, RESTRT, LDW, LDH);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*(RESTRT+4), sizeof(double));
  H = (double *) calloc(LDH*(RESTRT+2), sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL || H==NULL)
    errorexit("templates_gmres_wrapper: out of memory for X, B, WORK, H");


  /* setup global vars and functions needed in matvec and psolve */
  lop_fortemplates	= lop;
  precon_fortemplates	= precon;
  r_fortemplates	= r;
  x_fortemplates	= x;
  c1_fortemplates	= c1;
  c2_fortemplates	= c2;
  dim_fortemplates	= N;

  /* setup local B and X */
  copy_vl_into_array_outervlloop(b, B, N);
  copy_vl_into_array_outervlloop(x, X, N);


#ifdef TEMPLATES
  /* call gmres from templates */
  gmres_(&N, B, X, &RESTRT, WORK, &LDW, H, &LDH, &ITER, &RESID,
          matvec, psolve, &INFO);
#else
  errorexit("templates_gmres_wrapper: to compile with templates "
            "use MyConfig with\n"
            "DFLAGS += -DTEMPLATES\n"
            "SPECIALLIBS += -l~/Packages/dctemplates/libiteratortemplatates.a "
            "-l~/Packages/dctemplates/F2CLIBS/libF77.a "
            "-l~/Packages/dctemplates/F2CLIBS/libI77.a");
#endif

  /* read out vlx and normres */
  copy_array_into_vl_outervlloop(x, X, N);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);
  free(H);

  if(pr) printf("templates_gmres_wrapper: ITER=%d RESID=%g INFO=%d\n",
                ITER, RESID, INFO);

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}


/* call bicgstab_ from templates */
int templates_bicgstab_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  int N; /* dim of matrix */
  double *B;
  double *X;
  double *WORK;		int LDW;
  int ITER;
  double RESID;
  int INFO;

  /* set int vars */
  N = 0 ;
  for(j = 0; j < b->n; j++)
    forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = (N) + 1; 
  ITER = itmax;
  RESID = tol;
  if(pr) printf("templates_bicgstab_wrapper: itmax=%d tol=%g "
                "N=%d LDW=%d\n", itmax, tol, N, LDW);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*7, sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_bicgstab_wrapper: out of memory for X, B, WORK");


  /* setup global vars and functions needed in matvec and psolve */
  lop_fortemplates	= lop;
  precon_fortemplates	= precon;
  r_fortemplates	= r;
  x_fortemplates	= x;
  c1_fortemplates	= c1;
  c2_fortemplates	= c2;
  dim_fortemplates	= N;

  /* setup local B and X */
  copy_vl_into_array_outervlloop(b, B, N);
  copy_vl_into_array_outervlloop(x, X, N);


#ifdef TEMPLATES
  /* call bicgstab from templates */
  bicgstab_(&N, B, X, WORK, &LDW, &ITER, &RESID, matvec, psolve, &INFO);
#else
  errorexit("templates_bicgstab_wrapper: to compile with templates "
            "use MyConfig with\n"
            "DFLAGS += -DTEMPLATES\n"
            "SPECIALLIBS += -l~/Packages/dctemplates/libiteratortemplatates.a "
            "-l~/Packages/dctemplates/F2CLIBS/libF77.a "
            "-l~/Packages/dctemplates/F2CLIBS/libI77.a");
#endif

  /* read out vlx and normres */
  copy_array_into_vl_outervlloop(x, X, N);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  if(pr) printf("templates_bicgstab_wrapper: ITER=%d RESID=%g INFO=%d\n",
                ITER, RESID, INFO);

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}
