/* Poisson2.c */
/* Wolfgang Tichy 2016 */

/* Example module, where we solve two Poisson equations for Psi and Chi:
   The eqs are in F_Poisson2 and the linearized versions are in J_Poisson2 .
   The domain is a spherical shell. Psi and Chi are equal to different
   constants on the inner and outer shell boundary. */

#include "sgrid.h"
#include "Poisson2.h"



void Poisson2_set_BCs(tVarList *VLFu, tVarList *VLu, tVarList *VLuAll,
                      int nonlin);
void Precon_fd_Poisson2_UMFPACK(tVarList *vlx, tVarList *vlb,
                               tVarList *vlc1, tVarList *vlc2);


/* initialize Poisson2 */
int Poisson2_startup(tGrid *grid)
{
  int b;
  printf("Initializing Poisson2:\n");

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("Poisson2_Psi", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Poisson2_Chi", 0, 1, 1.0);

  /* enable all vars */
  enablevar(grid, Ind("Poisson2_Psi"));
  enablevar(grid, Ind("Poisson2_rh1"));
  enablevar(grid, Ind("Poisson2_Psix"));
  enablevar(grid, Ind("Poisson2_Psixx"));
  enablevar(grid, Ind("Poisson2_Err_Psi"));
  enablevar(grid, Ind("Poisson2_Chi"));
  enablevar(grid, Ind("Poisson2_rh2"));
  enablevar(grid, Ind("Poisson2_Chix"));
  enablevar(grid, Ind("Poisson2_Chixx"));
  enablevar(grid, Ind("Poisson2_Err_Chi"));
  enablevar(grid, Ind("Poisson2_temp1"));
  enablevar(grid, Ind("Poisson2_temp2"));
  enablevar(grid, Ind("Poisson2_temp3"));
  enablevar(grid, Ind("Poisson2_temp4"));
  enablevar(grid, Ind("Poisson2_temp5"));
  enablevar(grid, Ind("Poisson2_temp6"));
  
  /* set initial values in all in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *Psi = box->v[Ind("Poisson2_Psi")];
    double *Chi = box->v[Ind("Poisson2_Chi")];
    double *rh1 = box->v[Ind("Poisson2_rh1")];
    double *rh2 = box->v[Ind("Poisson2_rh2")];

    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      /* set Psi and Chi */
      if(Getv("Poisson2_grid", "SphericalDF"))
      {
        Psi[i] = x*y*z; // 1.0/sqrt(x*x + y*y + z*z);
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
      }
    }
  }

  return 0;
}


/* Solve the Equations */
int Poisson2_solve(tGrid *grid)
{
  int    itmax = Geti("Poisson2_itmax");
  double tol   = Getd("Poisson2_tol");
  int    linSolver_itmax  = Geti("Poisson2_linSolver_itmax");
  double linSolver_tolFac = Getd("Poisson2_linSolver_tolFac");
  double linSolver_tol    = Getd("Poisson2_linSolver_tol");
  double normresnonlin;
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  void (*Precon)(tVarList *vlx, tVarList *vlb, tVarList *vlc1, tVarList *vlc2);
  tVarList *VLu, *VLFu, *VLuAll;
  tVarList *VLlu, *VLJlu, *VLluAll;
  tVarList *VLrhs;

  /* allocate varlists */
  VLu  = vlalloc(grid);
  VLuAll= vlalloc(grid);

  /* add Poisson2_Psi and Poisson2_Chi to VLu */
  vlpush(VLu, Ind("Poisson2_Psi"));
  vlpush(VLu, Ind("Poisson2_Chi"));

  /* add Poisson2_Psi and Poisson2_Chi and their derivs to VLuAll */
  vlpush(VLuAll, Ind("Poisson2_Psi"));
  vlpush(VLuAll, Ind("Poisson2_Psix"));
  vlpush(VLuAll, Ind("Poisson2_Psixx"));
  vlpush(VLuAll, Ind("Poisson2_Chi"));
  vlpush(VLuAll, Ind("Poisson2_Chix"));
  vlpush(VLuAll, Ind("Poisson2_Chixx"));
  /* add Poisson2_rh1 and Poisson2_rh2 to VLuAll */
  vlpush(VLuAll, Ind("Poisson2_rh1"));
  vlpush(VLuAll, Ind("Poisson2_rh2"));

  /* allocate memory */
  enablevarlist(VLu);
  enablevarlist(VLuAll);

  /* now duplicate VLu to get VLFu */
  VLFu = AddDuplicateEnable(VLu, "_Err");

  /* now duplicate VLJlu, VLu for linarized Eqs. */
  VLJlu   = AddDuplicateEnable(VLFu, "_l");
  VLlu    = AddDuplicateEnable(VLu,  "_l");
  VLluAll = AddDuplicateEnable(VLuAll, "_l");

  /* choose linear solver */
  if(Getv("Poisson2_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("Poisson2_linSolver", "bicgstab_with_fd_UMFPACK_precon"))
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
  else if(Getv("Poisson2_linSolver", "templates_GMRES_with_fd_UMFPACK_precon"))
    linear_solver=templates_gmres_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson2_linSolver", "templates_BICGSTAB_with_fd_UMFPACK_precon"))
    linear_solver=templates_bicgstab_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson2_linSolver", "templates_CGS_with_fd_UMFPACK_precon"))
    linear_solver=templates_cgs_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson2_linSolver", "LAPACK"))
    linear_solver=LAPACK_dgesv_wrapper;
  else if(Getv("Poisson2_linSolver", "templates_GMRES"))
    linear_solver=templates_gmres_wrapper;
  else if(Getv("Poisson2_linSolver", "templates_BICGSTAB"))
    linear_solver=templates_bicgstab_wrapper;
  else if(Getv("Poisson2_linSolver", "templates_CGS"))
    linear_solver=templates_cgs_wrapper;
  else if(Getv("Poisson2_linSolver", "UMFPACK"))
    linear_solver=UMFPACK_solve_wrapper;
  else if(Getv("Poisson2_linSolver", "templates_GMRES_with_Jacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_Jacobi_precon;
  else if(Getv("Poisson2_linSolver", "templates_GMRES_with_SOR_precon"))
    linear_solver=templates_gmres_wrapper_with_SOR_precon;
  else if(Getv("Poisson2_linSolver", "SOR"))
    linear_solver=SOR_Iterator;
  else if(Getv("Poisson2_linSolver", "templates_GMRES_with_BlockJacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson2_linSolver", "ZIB_GMRES_with_BlockJacobi_precon"))
    linear_solver=ZIB_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson2_linSolver", "ZIB_GBIT_with_BlockJacobi_precon"))
    linear_solver=ZIB_gbit_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson2_linSolver", "ZIB_PCG_with_BlockJacobi_precon"))
    linear_solver=ZIB_pcg_wrapper_with_BlockJacobi_precon;
  else
    errorexit("Poisson2_solve: unknown Poisson2_linSolver");

  /* choose preconditioner */
  if(Getv("Poisson2_linSolver_Precon", "I"))
    Precon=Preconditioner_I;
  else if(Getv("Poisson2_linSolver_Precon", "fd_UMFPACK"))
    Precon=Precon_fd_Poisson2_UMFPACK;
  else if(Getv("Poisson2_linSolver_Precon", "templates"))
    Precon=templates_Preconditioner_for_templates_solver;
  else
    errorexit("Poisson2_solve: unknown Poisson2_linSolver_Precon");

  /* call Newton solver */
  Newton(F_Poisson2, J_Poisson2, VLu, VLFu, VLuAll, VLluAll,
         itmax, tol, &normresnonlin, 1,
         linear_solver, Precon, VLlu, VLJlu, VLuAll, VLluAll,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);

  /* free varlists */     
  VLDisableFree(VLlu);
  VLDisableFree(VLluAll);
  VLDisableFree(VLJlu);     
  vlfree(VLu);
  vlfree(VLuAll);
  vlfree(VLFu);
        
  return 0;
}


/* compute absolute error in ANALYSIS */
int Poisson2_analyze(tGrid *grid)
{
  int b;
  printf("Poisson2: computing absolute error\n");
  
  /* set initial data in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *Psi    = box->v[Ind("Poisson2_Psi")];
    double *PsiErr = box->v[Ind("Poisson2_Psi_Err")];
    double *Chi    = box->v[Ind("Poisson2_Chi")];
    double *ChiErr = box->v[Ind("Poisson2_Chi_Err")];


    /* subtract true values */
    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      if(Getv("Poisson2_grid", "SphericalDF"))
      {
        PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
        ChiErr[i] = Chi[i]-2.0/sqrt(x*x + y*y + z*z);
      }
    }
  }
  return 0;
}



/* evaluate Poisson eqn for VLu */
void F_Poisson2(tVarList *VLFu, tVarList *VLu,
               tVarList *VLuAll, tVarList *VLluAll)
{
  tGrid *grid = VLu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[VLFu->index[0]];
    double *Psi  = box->v[VLu->index[0]];
    double *Psixx = box->v[VLuAll->index[4]];
    double *Psiyy = box->v[VLuAll->index[7]];
    double *Psizz = box->v[VLuAll->index[9]];
    double *FChi = box->v[VLFu->index[1]];
    double *Chi  = box->v[VLu->index[1]];
    double *Chixx = box->v[VLuAll->index[14]];
    double *Chiyy = box->v[VLuAll->index[17]];
    double *Chizz = box->v[VLuAll->index[19]];
    double *rh1  = box->v[VLuAll->index[20]];
    double *rh2  = box->v[VLuAll->index[21]];
    int i;

    /* compute the derivs */
    D_and_DD_of_S(box, VLu->index[0], VLuAll->index[1], VLuAll->index[4]);
    D_and_DD_of_S(box, VLu->index[1], VLuAll->index[11], VLuAll->index[14]);

    /* Poisson eqs */
    forallpoints(box, i)
    {
      FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i] - rh1[i];
      FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i] - rh2[i];
    }
  }

  /* BCs */
  Poisson2_set_BCs(VLFu, VLu, VLuAll, 1);
}


/* evaluate linearized eqns */
void J_Poisson2(tVarList *VLJlu, tVarList *VLlu,
               tVarList *VLuAll, tVarList *VLluAll)
{
  tGrid *grid = VLlu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *JlPsi = box->v[VLJlu->index[0]];
    double *lPsi  = box->v[VLlu->index[0]];
    double *lPsixx = box->v[VLluAll->index[4]];
    double *lPsiyy = box->v[VLluAll->index[7]];
    double *lPsizz = box->v[VLluAll->index[9]];
    double *JlChi = box->v[VLJlu->index[1]];
    double *lChi  = box->v[VLlu->index[1]];
    double *lChixx = box->v[VLluAll->index[14]];
    double *lChiyy = box->v[VLluAll->index[17]];
    double *lChizz = box->v[VLluAll->index[19]];
    double *rh1  = box->v[VLluAll->index[20]];
    double *rh2  = box->v[VLluAll->index[21]];
    int i;

    /* compute the derivs */
    D_and_DD_of_S(box, VLlu->index[0], VLluAll->index[1], VLluAll->index[4]);
    D_and_DD_of_S(box, VLlu->index[1], VLluAll->index[11], VLluAll->index[14]);

    /* Poisson eqs */
    forallpoints(box, i)
    {
      JlPsi[i] = lPsixx[i] + lPsiyy[i] + lPsizz[i];
      JlChi[i] = lChixx[i] + lChiyy[i] + lChizz[i];
    }
  }

  /* BCs */
  Poisson2_set_BCs(VLJlu, VLlu, VLluAll, 0);
}



/* switch to fin. diff. and use UMFPACK as lin solver as precon */
void Precon_fd_Poisson2_UMFPACK(tVarList *vlx, tVarList *vlb,
                               tVarList *vlc1, tVarList *vlc2)
{
  tVarList *vlr = AddDuplicateEnable(vlx, "_Precon_r");
  tGrid *grid = vlx->grid;
  tGrid *grid_bak=make_empty_grid(grid->nvariables, 0);
  int INFO;
  int itmax = 1;  /* UMPACK does not iterate */
  double tol = 0; /* UMPACK has no tol, it solves explicitly */
  double normres;

  /* save current grid in grid_bak and then convert grid to fin. diff. */
  copy_grid_withoutvars(grid, grid_bak, 0);
  convert_grid_to_fd(grid);

  /* when we call J_Poisson2, we evaluate lin. eqns, but now with fd */
  INFO = UMFPACK_solve_wrapper(vlx, vlb, vlr, vlc1, vlc2, 
                               itmax, tol, &normres,
                               J_Poisson2, Preconditioner_I);
  /* restore grid */
  copy_grid_withoutvars(grid_bak, grid, 0);
  free_grid(grid_bak);
  VLDisableFree(vlr);
}


/* set BCs for a varlist */
void Poisson2_set_BCs(tVarList *VLFu, tVarList *VLu, tVarList *VLuAll,
                      int nonlin)
{
  tGrid *grid = VLu->grid;
  int b;
  int vind;

  for(vind=0; vind<VLu->n; vind++)
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[VLFu->index[vind]];
      double *Psi  = box->v[VLu->index[vind]];
      double *Psix = box->v[VLuAll->index[9*vind]];
      double *Psiy = box->v[VLuAll->index[9*vind+1]];
      double *Psiz = box->v[VLuAll->index[9*vind+2]];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;

      /* BCs */
      if(Getv("Poisson2_grid", "SphericalDF"))
      {
        forplane1(i,j,k, n1,n2,n3, 0)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*(vind+1)*nonlin;

        forplane1(i,j,k, n1,n2,n3, n1-1)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*(vind+1)*nonlin;
      }

    } /* end forallboxes */
}
