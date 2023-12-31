/* Newton_linesrch.c */
/* Wolfgang Tichy 5/2013 */

#include "sgrid.h"
#include "GridIterators.h"


#define COPY_ARRAY_INTO_VL copy_array_into_varlist
#define COPY_VL_INTO_ARRAY copy_varlist_into_array



/* structure we use to store varlists and funcs */
typedef struct tNEWTONARGS {
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2);
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2);
  tVarList *vlFu;
  tVarList *vlu;
  tVarList *vlc1;
  tVarList *vlc2;
  tVarList *vlres;
  tVarList *vlJdu;
  tVarList *vldu;
  tVarList *vld1;
  tVarList *vld2;
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  void (*Precon)(tVarList *Hinv_v, tVarList *v, tVarList *, tVarList *);
  int pr;
} tNewtonArgs;


/* vector function for WT_newton */
void F_x_for_Newton_linesrch(int n, double *x, double *Fx, void *par)
{
  tNewtonArgs *p =  (tNewtonArgs *) par;

  /* copy x into vlu. NOTE: x[1...n] not x[0...n-1] */
  COPY_ARRAY_INTO_VL(x+1, p->vlu);

  /* call func Fu */
  p->Fu(p->vlFu, p->vlu, p->vlc1, p->vlc2);

  /* set Fx array */
  COPY_VL_INTO_ARRAY(p->vlFu, Fx+1);
}

/* Jacobian of vector function times dx (evaluated at x) */
void J_x_dx_for_Newton_linesrch(int n, double *dx, double *Jdx,
                                double *x, void *par)
{
  tNewtonArgs *p =  (tNewtonArgs *) par;

  errorexit("J_x_dx_for_Newton_linesrch should never be called");

  /* copy dx,x into vldu,vlu. NOTE: x[1...n] not x[0...n-1] */
  COPY_ARRAY_INTO_VL(dx+1, p->vldu);
  COPY_ARRAY_INTO_VL(x+1,  p->vlu);

  /* call func Fu */
  p->Jdu(p->vlJdu, p->vldu, p->vld1, p->vld2);

  /* set Jdx array */
  COPY_VL_INTO_ARRAY(p->vlJdu, Jdx+1);
}

/* linear solver */
int linSol_for_Newton_linesrch(int n, double *b, double *dx,
           void (*J_dx)(int, double *, double *, double *, void *),
           int (*precon)(int n, double *b, double *dx, double *x, void *par),
           double *x, void *par, int itmax, double tol)
{
  tNewtonArgs *p =  (tNewtonArgs *) par;
  int ret;
  double res;

  /* copy x,b into vldu,vlFu. NOTE: dx[1...n] not dx[0...n-1] */
  COPY_ARRAY_INTO_VL(dx+1, p->vldu);
  COPY_ARRAY_INTO_VL(b+1,  p->vlFu);

  /* get current norm */
  res = norm2(p->vlFu);
  if(p->pr) printf("Newton_linesrch:  Newton residual = %e\n", res);

  /* call linear solver */
  res=0;
  if(p->pr)
  {
    prdivider(0);
    printf("linSol_for_Newton_linesrch: itmax=%d tol=%e\n", itmax,
    tol);
  }
  ret = p->linSolver(p->vldu, p->vlFu, p->vlres, p->vld1, p->vld2,
                     itmax, tol, &res, p->Jdu, p->Precon);
  if(p->pr)
    printf("linSol_for_Newton_linesrch: ret=%d res=%e\n", ret, res);

  /* set dx array */
  COPY_VL_INTO_ARRAY(p->vldu, dx+1);
  //printf("dx[666+1]=%g\n", dx[666+1]);
  return ret;
}

/* this precon should not be needed because the linear solver above
   calls p->Precon, and not this one */
int precon_for_Newton_linesrch(int n, double *b, double *dx,
                               double *x, void *par)
{
  errorexit("precon_for_Newton_linesrch should never be called");
  return -1;
}



/* interface that takes same args as Newton but then uses 
int WT_newton(double *x, int n, int *check,
        void (*F_x)(int, double *x, double *Fx, void *par),
        void (*J_x_dx)(int, double *dx, double *Jdx, double *x, void *par),
        void *par, int MAXITS, double TOLF,
        int (*linSol)(int n, double *b, double *dx,
            void (*Jdx)(int, double *, double *, double *, void *),
            int (*precon)(int n, double *b, double *dx, double *x, void *par),
            double *x, void *par, int itmax, double tol),
        int (*precon)(int n, double *b, double *dx, double *x, void *par),
        int linitmax, double lintolfac);
   to actually find the roots */
int Newton_linesrch(
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlu, tVarList *vlFu, tVarList *vlc1, tVarList *vlc2,
  int itmax, double tol, double *normres, int pr,
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  void (*Precon)(tVarList *Hinv_v, tVarList *v, tVarList *, tVarList *),
  tVarList *vldu, tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int linSolv_itmax, double linSolv_tolFac, double linSolv_tol)
{
  tGrid *grid = vlFu->grid;
  int i, n, inewton, check;
  double *x;
  tNewtonArgs Nargs[1];
  void *par;
  double res;

  /* save all args */
  Nargs->Fu   = Fu;
  Nargs->Jdu  = Jdu;
  Nargs->vlFu = vlFu;
  Nargs->vlu  = vlu;
  Nargs->vlc1 = vlc1;
  Nargs->vlc2 = vlc2;
  Nargs->vlres= vlres;
  Nargs->vldu = vldu;
  Nargs->vld1 = vld1;
  Nargs->vld2 = vld2;
  Nargs->linSolver = linSolver;
  Nargs->Precon = Precon;
  Nargs->pr   = pr;
  par = (void *) Nargs;

  /* set number of unknows WT_newton solves for */
  n = 0;
  forallboxes(grid,i)  n += grid->box[i]->nnodes;
  n = (vlFu->n) * n;

  /* make array for x. NOTE: x[1...n] not x[0...n-1] */
  x = (double *) calloc(n+1, sizeof(double));

  /* copy vlu into x array. NOTE: x[1...n] not x[0...n-1] */
  COPY_VL_INTO_ARRAY(vlu, x+1);

  /* solve with Newton-Raphson iterations: */
  if(pr) printf("Newton_linesrch:  starting Newton iterations, itmax=%d tol=%g\n",
                itmax, tol);

  /* call my Newton solver that deals with arrays */
  inewton = WT_newton(x, n, &check, F_x_for_Newton_linesrch,
                      J_x_dx_for_Newton_linesrch, par, itmax, tol,
                      linSol_for_Newton_linesrch, precon_for_Newton_linesrch,
                      linSolv_itmax, linSolv_tolFac);
  /* copy x array into vlu. NOTE: x[1...n] not x[0...n-1] */
  COPY_ARRAY_INTO_VL(x+1, vlu);

  /* get residual */
  Fu(vlFu, vlu, vlc1, vlc2);
  res = norm2(vlFu);
  if(pr) printf("Newton_linesrch: after %d steps"
                "  Newton residual = %e\n", inewton, res);

  /* warn if we didn't converge */
  if (inewton >= itmax)
  {
    printf("Newton_linesrch warning: *** Too many Newton steps! ");
    if(res <= tol) printf("*** \n");
    else	   printf("Tolerance goal not reached! *** \n");
  }
  if(inewton<0 || check)
    printf("Newton_linesrch warning: inewton=%d  check=%d\n", inewton, check);

  /* mem for x */
  free(x);
  return inewton;
}
