/* Newton.c */
/* Wolfgang Tichy 7/2003 */

#include "sgrid.h"
#include "GridIterators.h"


/* structure we can use to store varlists and funcs */
typedef struct tNEWTONSTEPVARS {
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2);
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2);
  tVarList *vlFu;
  tVarList *vlu;
  tVarList *vlc1;
  tVarList *vlc2;
  tVarList *vlres;
  tVarList *vldu;
  tVarList *vld1;
  tVarList *vld2;
  tVarList *vltemp;
} tNewtonStepVars;


/* function pointer that is called before the linear solver,
   if it does not point to NULL */
int (*Newton_do_before_linSolver)(tGrid *g);

/* funcs */
void do_partial_Newton_step(tVarList *vlu, double lambda, tVarList *vldu);
void do_Newton_step(tVarList *vlu, tVarList *vldu, double oldres,
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlFu,  tVarList *vlc1, tVarList *vlc2,
  tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int pr);




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
  if(pr) printf("Newton:  starting Newton iterations, itmax=%d tol=%g\n",
                itmax, tol);
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

    /* user defined function that is called before the linear solver,
       if it does not point to NULL */
    if(Newton_do_before_linSolver != NULL)
      Newton_do_before_linSolver(grid);

    /* solve linear equation */
    lin_its=linSolver(vldu, vlFu, vlres, vld1, vld2, 
                      linSolv_itmax, 
                      max2( (*normres)*linSolv_tolFac, linSolv_tol ),
                      &lin_normres,
	              Jdu, linPrecon);
    /* if(pr) printf("Newton: after linSolver: %e\n",lin_normres); */

    /* do Newton step: u^{n+1} = u^{n} - du */
    /* do_partial_Newton_step(vlu, 1.0, vldu); */
    do_Newton_step(vlu, vldu, res, Fu, Jdu,
                   vlFu,  vlc1, vlc2, vlres, vld1, vld2, pr);

    /* sync vlu. sync is not needed if du is synced */
    /* bampi_vlsynchronize(vlu); */

    /* boundary conditions for vlu are set in the function
       Fu(vlFu, vlu, vlc1, vlc2);  supplied by the user.         */
  } 

  /* warn if we didn't converge */
  if (inewton >= itmax)
  {
    Fu(vlFu, vlu, vlc1, vlc2);
    res = norm2(vlFu);
    *normres = res;
    printf("Newton warning: *** Too many Newton steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
    printf("Newton: Residual after %d Newton steps:"
           "  Newton residual = %e\n", inewton, *normres);
  }

  return inewton;
}


/* partial Newton step */
/* do Newton step: u^{n+1} = u^{n} - lambda*du */
/* Note: if lambda=1 this is a full Newton step. */
void do_partial_Newton_step(tVarList *vlu, double lambda, tVarList *vldu)
{
  tGrid *grid = vlu->grid;
  int i, j, b;
  
  for(j = 0; j < vlu->n; j++)
    forallboxes(grid,b)
    {
      tBox *box = grid->box[b];
      double *u  = box->v[vlu ->index[j]]; 
      double *du = box->v[vldu->index[j]]; 

      forallpoints(box,i)
        u[i] -= lambda*du[i]; /* do Newton step: u^{n+1} = u^{n} - lambda*du */
    }
}

/* function for brent, which gives Newton residual as a function of 
   lambda only */
double Newton_residual_of_lambda(double lambda, void *p)
{
  tNewtonStepVars *pars;
  pars = (tNewtonStepVars *) p;

  /* set vlu = vltemp */  
  vlcopy(pars->vlu, pars->vltemp);
  do_partial_Newton_step(pars->vlu, lambda, pars->vldu);
  pars->Fu(pars->vlFu, pars->vlu, pars->vlc1, pars->vlc2);
  return norm2(pars->vlFu);
}

/* do one Newton step and with backtracking if needed */
void do_Newton_step(tVarList *vlu, tVarList *vldu, double oldres,
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlFu,  tVarList *vlc1, tVarList *vlc2,
  tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int pr)
{
  double cl = 1.0;

  /* do full Newton step */
  do_partial_Newton_step(vlu, cl, vldu);

  /* do we use backtracking? */
  if(Getv("GridIterators_Newtonstep", "backtrack"))
  {
    void *p;
    tNewtonStepVars pars[1]; /* create tNewtonStepVars *pars but with memory */
    double al = 0.0;
    double bl = 1e-6;
    double tol = 0.01;
    double lmin, lambda;
    double fx;  
    tVarList *vltemp;
    double resb, resc;

    /* res of full Newton step */
    Fu(vlFu, vlu, vlc1, vlc2);
    resc = norm2(vlFu);
    if(resc<oldres)
    {
      vlsetconstant(vldu, 0.0); 
      return; /* stick with full Newton step if it decreases res */
    }

    /* do Newton step only to lambda = bl */
    do_partial_Newton_step(vlu, bl-cl, vldu);
    Fu(vlFu, vlu, vlc1, vlc2);
    resb = norm2(vlFu);
    if(resb>=oldres)
    { /* errorexit("do_Newton_step: (res at bl) > (res at 0)"); */
      lambda = 1e-4;
      if(pr)
      {
        printf("do_Newton_step: Warning: "
               "(res at lambda=bl=%g) > (res at lambda=0)\n", bl);
        printf("do_Newton_step: trying to escape local min. with "
               "lambda=%g\n", lambda);
      }
      do_partial_Newton_step(vlu, lambda-bl, vldu);
      vlsetconstant(vldu, 0.0);
      return; /* do Newton step by amount lambda, even though res is worse */
    }
    /* save old vlu in vltemp */
    vltemp = AddDuplicateEnable(vlu, "_GridIterators_Newtonstep_temp");
    do_partial_Newton_step(vlu, -bl, vldu); /* go back to lambda=0 */
    vlcopy(vltemp, vlu);
            
    /* set pars for Newton_residual_of_lambda */
    pars->Fu   = Fu;
    pars->Jdu  = Jdu;
    pars->vlFu = vlFu;
    pars->vlu  = vlu;
    pars->vlc1 = vlc1;
    pars->vlc2 = vlc2;
    pars->vlres= vlres;
    pars->vldu = vldu;
    pars->vld1 = vld1;
    pars->vld2 = vld2;
    pars->vltemp = vltemp;
    p = (void *) pars;

    /* call brent with Newton_residual_of_lambda to find the lambda=lmin 
       where the res in min */
    fx = brent_with_pointer_to_pars(al, bl, cl, Newton_residual_of_lambda,
                                    tol, &lmin, p);
    if(lmin < bl) lambda = bl;
    else          lambda = lmin;
    if(pr) printf("do_Newton_step:  backtracking to lambda=%g\n", lambda);
    
    /* do reduced Newton step */
    vlcopy(vlu, vltemp);
    do_partial_Newton_step(vlu, lambda, vldu);

    /* free varlist vltemp */
    vlfree(vltemp);
  }
  /* reset du to zero */
  vlsetconstant(vldu, 0.0);
}

/* Newton startup initialization */
int Init_Newton(tGrid *grid)
{
  /* set function pointer to NULL */
  Newton_do_before_linSolver = NULL;

  return 0;
}
