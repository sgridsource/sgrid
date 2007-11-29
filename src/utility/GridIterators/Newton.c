/* Newton.c */
/* Wolfgang Tichy 7/2003 */

#include "sgrid.h"
#include "GridIterators.h"


/* Newton Raphson solver, which solves F(u)=0 on a grid.
   It takes the function F(u), its linearization J(du), some var lists and
   pars as well as a linear Solver and a Preconditioner as arguments.      */
int Newton(
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlu, tVarList *vlFu, tVarList *vlc1, tVarList *vlc2,
  int itmax, double tol, double *normres, int pr,
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  void (*linPrecon)(tVarList *Hinv_v, tVarList *v, tVarList *, tVarList *),
  tVarList *vldu, tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int linSolv_itmax, double linSolv_tolFac, double linSolv_tol)
{
  tGrid *grid = vlFu->grid;
  int i, j, inewton, b;
  double res, lin_normres = 0;
  int lin_its = 0;
  
  /* solve with Newton-Raphson iterations: */
  if(pr) printf("Newton:  starting Newton iterations ...\n");
  for (inewton = 0; inewton < itmax; inewton++)
  {
    /* compute vlFu = F(u) */
    Fu(vlFu, vlu, vlc1, vlc2);
    res = norm2(vlFu);
    *normres = res;

    if(pr)
    {
      printf("Newton step %d: linSolves=%d res=%.4e"
             "  Newton residual = %.4e\n",
             inewton, lin_its, lin_normres, *normres);
      fflush(stdout);
    }
    
    if (*normres <= tol) break;

    /* solve linear equation */
    lin_its=linSolver(vldu, vlFu, vlres, vld1, vld2, 
                      linSolv_itmax, 
                      max2( (*normres)*linSolv_tolFac, linSolv_tol ),
                      &lin_normres,
	              Jdu, linPrecon);
    /* if(pr) printf("Newton: after linSolver: %e\n",lin_normres); */

    /* do Newton step: u^{n+1} = u^{n} - du */
    for(j = 0; j < vlu->n; j++)
      forallboxes(grid,b)
      {
        tBox *box = grid->box[b];
	double *u  = box->v[vlu ->index[j]]; 
	double *du = box->v[vldu->index[j]]; 

        forallpoints(box,i)
        {
          u[i] -= du[i];  /* do Newton step: u^{n+1} = u^{n} - du */
          du[i] = 0;      /* reset du to zero */
        }
      }

    /* sync vlu. sync is not needed if du is synced */
    /* bampi_vlsynchronize(vlu); */

    /* boundary conditions for vlu are set in the function
       Fu(vlFu, vlu, vlc1, vlc2);  supplied by the user.         */
  } 

  /* warn if we didn't converge */
  if (inewton >= itmax)
  {
    printf("Newton warning: *** Too many Newton steps! ");
    printf("Tolerance goal not reached! *** \n");
    Fu(vlFu, vlu, vlc1, vlc2);
    res = norm2(vlFu);
    *normres = res;
    printf("Newton: Residual after %d Newton steps:"
           "  Newton residual = %e\n", inewton, *normres);
  }

  return inewton;
}
