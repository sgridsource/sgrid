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
int (*Newton_do_before_linSolver)(tGrid *grid, tNewtonResults inp,
                                  tNewtonArgPointers out);

/* funcs */
void do_partial_Newton_step(tVarList *vlu, double lambda, tVarList *vldu);
void do_random_Newton_step(tVarList *vlu, double eps, tVarList *vldu);
double do_Newton_step(tVarList *vlu, tVarList *vldu, double oldres,
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
  int inewton;
  double res, lin_normres = 0;
  int lin_its = 0;
  double lambda = 0.0;
  double ms = Getd("GridIterators_Newton_minstep");
  tNewtonResults before_lin_in;
  tNewtonArgPointers before_lin_out;
  tLinSolver linSolver1;
  int linSolverOK, linSolverPos;

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
    before_lin_in.lambda = lambda;
    before_lin_out.linSolver = &linSolver;
    if(Newton_do_before_linSolver != NULL)
      Newton_do_before_linSolver(grid, before_lin_in, before_lin_out);

    /* Call some linear solver and do Newton step as needed.
       We can switch the linear solver if step fails. */
    linSolver1 = linSolver; /* try first with linSolver passed in to Newton */
    linSolverPos = 0;       /* position of next linSolver in GridIterators_Newton_linSolvers */
    do /* loop while !linSolverOK */
    {
      /* solve linear equation */
      lin_its=linSolver1(vldu, vlFu, vlres, vld1, vld2, linSolv_itmax, 
                         max2( (*normres)*linSolv_tolFac, linSolv_tol ),
                         &lin_normres, Jdu, linPrecon);
      linSolverOK=1;
      /* if(pr) printf("Newton: after linSolver: %e\n",lin_normres); */
      if(lin_its>=0) /* linSolver1 found a solution without problems */
      {
        /* do Newton step: u^{n+1} = u^{n} - du */
        /* do_partial_Newton_step(vlu, 1.0, vldu); */
        lambda = do_Newton_step(vlu, vldu, res, Fu, Jdu,
                                vlFu,  vlc1, vlc2, vlres, vld1, vld2, pr);
      }
      else /* linSolver1 returned error code */
      {
        /* reset du to zero */
        vlsetconstant(vldu, 0.0);
        /* call Jdu with vldu=0,
           which hopefully resets all derivs of vldu to zero */
        if(Getv("GridIterators_Newton_EndOfStep", "Jdu"))
          Jdu(vlres, vldu, vld1, vld2);

        /* signal that we stayed at previous vlu */
        lambda = -1.0;
        printf("Newton error: *** Linear solver failed and returned %d ***\n",
               lin_its);
        printf("  signaling that solution is unchanged by setting lambda=%g\n",
               lambda);
      }

      /* do we seem to be caught in a local min, i.e. close to lambda=-1 ? */
      if(fabs(lambda+1.0)<ms)
      {
        printf("lambda=%g : we may be stuck in a local minimum!\n", lambda);

        /* do random Newton step if we seem to be caught in a local min */
        if(Getv("GridIterators_Newton_atlocalMin","escapeMin"))
        {
          double eps = Getd("GridIterators_Newton_randomstepsize");
          printf(" ==> taking random Newton step of size eps=%g\n", eps);
          do_random_Newton_step(vlu, eps, vldu);
        }
        /* quit iterations by setting inewton = itmax */
        if(Getv("GridIterators_Newton_atlocalMin","quit"))
        {
          inewton = itmax;
          printf(" ==> quitting (by going to step %d)\n", inewton+1);
        }
        /* try another linSolver */
        if(Getv("GridIterators_Newton_atlocalMin","AltLinSolver"))
        {
          char *linSolvers = Gets("GridIterators_Newton_linSolvers");
          char *word = cmalloc(strlen(linSolvers)+1);
          int p0 = linSolverPos;
          linSolverPos = sscan_word_at_p(linSolvers, linSolverPos, word);
          if(linSolverPos != EOF)
          {
            printf(" ==> retrying with %s\n", word);
            printf("     (char %d-%d in GridIterators_Newton_linSolvers)\n",
                   p0, linSolverPos);
            linSolverOK=0;
            linSolver1 = get_linSolver_by_name(word);
          }
          free(word);
        }
      }
    } while(!linSolverOK);

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

/* random Newton step: u^{n+1} = u^{n} + uav * eps*r,
   where r is random in [-1,1] and eps<1 */
void do_random_Newton_step(tVarList *vlu, double eps, tVarList *vldu)
{
  tGrid *grid = vlu->grid;
  int i, j, b;
  double uav = norm2(vlu);
  char str[1000];

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];

    for(j = 0; j < vlu->n; j++)
    {
      double *u  = box->v[vlu ->index[j]];
      /* double *du = box->v[vldu->index[j]]; */

      forallpoints(box,i)
        u[i] += uav*eps*(RND()-0.5)*2.0; /* u^{n+1} = u^{n} + uav*eps*r */
    }
    /* make sure vlu is unique */
    snprintf(str, 999, "box%d_Coordinates", box->b);
    if(Getv(str, "SphericalDF"))  reset_doubleCoveredPoints(vlu);
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
double do_Newton_step(tVarList *vlu, tVarList *vldu, double oldres,
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlFu,  tVarList *vlc1, tVarList *vlc2,
  tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int pr)
{
  double lambda = 0.0;

  /* do full Newton step */
  do_partial_Newton_step(vlu, 1.0, vldu);
/*
if(strstr(VarName(vlu->index[0]), "BNSdata_Psi"))
{
tGrid *grid=vlu->grid;
//prvarlist(vldu);
//printvar(grid, "BNSdata_Psi_l");
grid->time = -66;
write_grid(grid);
exit(5);
}
*/
  /* do we use backtracking? */
  if(Getv("GridIterators_Newtonstep", "backtrack"))
  {
    void *p;
    tNewtonStepVars pars[1]; /* create tNewtonStepVars *pars but with memory */
    double al = -1.1;
    double bl = -0.5;
    double cl = -0.0;
    double tol = 0.0001;
    double lmin;
    double fx, fa,fb,fc;
    tVarList *vltemp;
    double resc;

    /* res of full Newton step */
    Fu(vlFu, vlu, vlc1, vlc2);
    resc = norm2(vlFu);

    /* backtrack only if Newton step increases res or optimal is on */
    if(resc>=oldres || Getv("GridIterators_Newtonstep", "optimal"))
    {
      /* save old vlu in vltemp */
      vltemp = AddDuplicateEnable(vlu, "_GridIterators_Newtonstep_temp");
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

      /* bracket the min */
      mnbrak_with_pointer_to_pars(&al, &bl, &cl, &fa, &fb, &fc,
                                  Newton_residual_of_lambda, p);

      /* call brent with Newton_residual_of_lambda to find the lambda=lmin 
         where the res in min */
      fx = brent_with_pointer_to_pars(al, bl, cl, Newton_residual_of_lambda,
                                      tol, &lmin, p);
      /* if(lmin < bl) lambda = bl;
         else          lambda = lmin; */
      lambda = lmin;
      if(pr)
      {
        printf("do_Newton_step:  mnbrak -> (al,bl,cl)=(%g,%g,%g), "
               "(fa,fb,fc)=(%.3g,%.3g,%.3g)\n", al,bl,cl, fa,fb,fc);
        printf("do_Newton_step:  brent tol=%g -> lmin=%g, fx=%.3g\n",
               tol, lmin, fx);
        printf("do_Newton_step:  backtracking by lambda=%g\n", lambda);
      }

      /* do reduced Newton step */
      vlcopy(vlu, vltemp);
      do_partial_Newton_step(vlu, lambda, vldu);

      /* free varlist vltemp */
      vlfree(vltemp);
    }
  }
  /* reset du to zero */
  vlsetconstant(vldu, 0.0);

  /* call Jdu with vldu=0, 
     which hopefully resets all derivs of vldu to zero */
  if(Getv("GridIterators_Newton_EndOfStep", "Jdu"))
    Jdu(vlres, vldu, vld1, vld2);

  return lambda;
}

/* Newton startup initialization */
int Init_Newton(tGrid *grid)
{
  /* set function pointer to NULL */
  Newton_do_before_linSolver = NULL;

  return 0;
}

/* get linear solver function pointer from string */
tLinSolver get_linSolver_by_name(char *name)
{
  tLinSolver linear_solver=NULL;

  /* choose linear solver */
  if(strcmp(name, "bicgstab")==0)
    linear_solver=bicgstab;
  else if(strcmp(name, "bicgstab_with_fd_UMFPACK_precon")==0)
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
  else if(strcmp(name, "LAPACK_dgesv_wrapper")==0)
    linear_solver=LAPACK_dgesv_wrapper;
  else if(strcmp(name, "templates_gmres_wrapper")==0)
    linear_solver=templates_gmres_wrapper;
  else if(strcmp(name, "templates_gmres_wrapper_with_fd_UMFPACK_precon")==0)
    linear_solver=templates_gmres_wrapper_with_fd_UMFPACK_precon;
  else if(strcmp(name, "templates_bicgstab_wrapper")==0)
    linear_solver=templates_bicgstab_wrapper;
  else if(strcmp(name, "templates_bicgstab_wrapper_with_fd_UMFPACK_precon")==0)
    linear_solver=templates_bicgstab_wrapper_with_fd_UMFPACK_precon;
  else if(strcmp(name, "templates_cgs_wrapper")==0)
    linear_solver=templates_cgs_wrapper;
  else if(strcmp(name, "templates_cgs_wrapper_with_fd_UMFPACK_precon")==0)
    linear_solver=templates_cgs_wrapper_with_fd_UMFPACK_precon;
  else if(strcmp(name, "UMFPACK_solve_wrapper")==0)
    linear_solver=UMFPACK_solve_wrapper;
  else if(strcmp(name, "SuiteSparseQR_solve_wrapper")==0)
    linear_solver=SuiteSparseQR_solve_wrapper;
  else if(strcmp(name, "UMFPACK_solve_forSortedVars_wrapper")==0)
    linear_solver=UMFPACK_solve_forSortedVars_wrapper;
  else if(strcmp(name, "bicgstab_with_Jacobi_precon")==0)
    linear_solver=bicgstab_with_Jacobi_precon;
  else if(strcmp(name, "templates_gmres_wrapper_with_Jacobi_precon")==0)
    linear_solver=templates_gmres_wrapper_with_Jacobi_precon;
  else if(strcmp(name, "templates_gmres_wrapper_with_BlockJacobi_precon")==0)
    linear_solver=templates_gmres_wrapper_with_BlockJacobi_precon;
  else if(strcmp(name, "templates_bicgstab_wrapper_with_Jacobi_precon")==0)
    linear_solver=templates_bicgstab_wrapper_with_Jacobi_precon;
  else if(strcmp(name, "templates_cgs_wrapper_with_Jacobi_precon")==0)
    linear_solver=templates_cgs_wrapper_with_Jacobi_precon;
  else if(strcmp(name, "templates_qmr_wrapper_with_Jacobi_precon")==0)
    linear_solver=templates_qmr_wrapper_with_Jacobi_precon;
  else if(strcmp(name, "templates_bicg_wrapper_with_Jacobi_precon")==0)
    linear_solver=templates_bicg_wrapper_with_Jacobi_precon;
  else if(strcmp(name, "templates_bicgstab_wrapper_with_BlockJacobi_precon")==0)
    linear_solver=templates_bicgstab_wrapper_with_BlockJacobi_precon;
  else if(strcmp(name, "ZIB_gmres_wrapper_with_BlockJacobi_precon")==0)
    linear_solver=ZIB_gmres_wrapper_with_BlockJacobi_precon;
  else if(strcmp(name, "ZIB_gbit_wrapper_with_BlockJacobi_precon")==0)
    linear_solver=ZIB_gbit_wrapper_with_BlockJacobi_precon;
  else if(strcmp(name, "ZIB_pcg_wrapper_with_BlockJacobi_precon")==0)
    linear_solver=ZIB_pcg_wrapper_with_BlockJacobi_precon;
  else if(strcmp(name, "templates_qmr_wrapper")==0)
    linear_solver=templates_qmr_wrapper;
  else if(strcmp(name, "templates_bicg_wrapper")==0)
    linear_solver=templates_bicg_wrapper;
  else if(strcmp(name, "bicgstab_with_SOR_precon")==0)
    linear_solver=bicgstab_with_SOR_precon;
  else if(strcmp(name, "templates_gmres_wrapper_with_SOR_precon")==0)
    linear_solver=templates_gmres_wrapper_with_SOR_precon;
  else if(strcmp(name, "templates_bicgstab_wrapper_with_SOR_precon")==0)
    linear_solver=templates_bicgstab_wrapper_with_SOR_precon;
  else if(strcmp(name, "templates_cgs_wrapper_with_SOR_precon")==0)
    linear_solver=templates_cgs_wrapper_with_SOR_precon;
  else
    errorexit("InitRealisticBBH: unknown RealisticBBH_linSolver");

  return linear_solver;
}
