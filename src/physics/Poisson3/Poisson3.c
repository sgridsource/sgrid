/* Poisson3.c */
/* Wolfgang Tichy 2017 */

/* Example module, where we solve two Poisson equations for Psi and Chi:
   The eqs are in F_Poisson3 and the linearized versions are in J_Poisson3 .
   The domain is a spherical shell. Psi and Chi are equal to different
   constants on the inner and outer shell boundary. */

#include "sgrid.h"
#include "Poisson3.h"



void Poisson3_set_BCs(tVarList *VLFu, tVarList *VLu, tVarList *VLuAll,
                      int nonlin);
void Precon_fd_Poisson3_UMFPACK(tVarList *vlx, tVarList *vlb,
                               tVarList *vlc1, tVarList *vlc2);


int Poisson3_initboxes(tGrid *grid)
{
  printf("Initializing Poisson3 boxes:\n");
  if(Getv("Poisson3_grid","2starcubes"))
  {
    /* test cubed spheres */
    two_full_cubes_touching_at_x0(grid, 0, 0.25, 0.05,0.15, 0.1,0.2);
  }
  return 0;
}

/* initialize Poisson3 */
int Poisson3_startup(tGrid *grid)
{
  int b;
  int pr=1;
  printf("Initializing Poisson3:\n");

  if(Getv("Poisson3_grid","2starcubes"))
  {
    /* set bface structures */
    set_touching_bfaces_of_boxes_with_same_facepoints(grid, 0, grid->nboxes);
    set_all_bfaces_with_ob_minus1_to_outerbound(grid, 0, grid->nboxes);
    if(pr) forallboxes(grid, b) printbfaces(grid->box[b]);
  }
  else
  {
    /* set bface structures that contain info about box boundaries */
    Coordinates_set_bfaces(grid);
  }

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("Poisson3_Psi", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Poisson3_Chi", 1, 1, 1.0);

  /* enable all vars */
  enablevar(grid, Ind("Poisson3_Psi"));
  enablevar(grid, Ind("Poisson3_rh1"));
  enablevar(grid, Ind("Poisson3_Psix"));
  enablevar(grid, Ind("Poisson3_Psixx"));
  enablevar(grid, Ind("Poisson3_Err_Psi"));
  enablevar(grid, Ind("Poisson3_Chi"));
  enablevar(grid, Ind("Poisson3_rh2"));
  enablevar(grid, Ind("Poisson3_Chix"));
  enablevar(grid, Ind("Poisson3_Chixx"));
  enablevar(grid, Ind("Poisson3_Err_Chi"));
  enablevar(grid, Ind("Poisson3_temp1"));
  enablevar(grid, Ind("Poisson3_temp2"));
  enablevar(grid, Ind("Poisson3_temp3"));
  enablevar(grid, Ind("Poisson3_temp4"));
  enablevar(grid, Ind("Poisson3_temp5"));
  enablevar(grid, Ind("Poisson3_temp6"));
  
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
    double *Psi = box->v[Ind("Poisson3_Psi")];
    double *Chi = box->v[Ind("Poisson3_Chi")];
    double *rh1 = box->v[Ind("Poisson3_rh1")];
    double *rh2 = box->v[Ind("Poisson3_rh2")];
    double sqrtPI3 = sqrt(PI*PI*PI);

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
      /* Soln of "Lapl Psi = rh1"  is "Psi = 1.0/(4*PI*r)*erf(r)" */
      Psi[i] = 0.0; // 1.0/sqrt(x*x + y*y + z*z);
      Chi[i] = 0.9;
      rh1[i] = -exp(-(x*x + y*y + z*z))/sqrtPI3;
      rh2[i] = 0.0;
      //Psi[i] = 0.5*b + (b+1)*x;
      //Chi[i] = b - 2*(b+1)*x;
    }
  }

  return 0;
}


/* Solve the Equations */
int Poisson3_solve(tGrid *grid)
{
  int    itmax = Geti("Poisson3_itmax");
  double tol   = Getd("Poisson3_tol");
  int    linSolver_itmax  = Geti("Poisson3_linSolver_itmax");
  double linSolver_tolFac = Getd("Poisson3_linSolver_tolFac");
  double linSolver_tol    = Getd("Poisson3_linSolver_tol");
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

  /* add Poisson3_Psi and Poisson3_Chi to VLu */
  vlpush(VLu, Ind("Poisson3_Psi"));
  vlpush(VLu, Ind("Poisson3_Chi"));

  /* add Poisson3_Psi and Poisson3_Chi and their derivs to VLuAll */
  vlpush(VLuAll, Ind("Poisson3_Psi"));
  vlpush(VLuAll, Ind("Poisson3_Psix"));
  vlpush(VLuAll, Ind("Poisson3_Psixx"));
  vlpush(VLuAll, Ind("Poisson3_Chi"));
  vlpush(VLuAll, Ind("Poisson3_Chix"));
  vlpush(VLuAll, Ind("Poisson3_Chixx"));
  /* add Poisson3_rh1 and Poisson3_rh2 to VLuAll */
  vlpush(VLuAll, Ind("Poisson3_rh1"));
  vlpush(VLuAll, Ind("Poisson3_rh2"));

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
  if(Getv("Poisson3_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("Poisson3_linSolver", "bicgstab_with_fd_UMFPACK_precon"))
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
  else if(Getv("Poisson3_linSolver", "templates_GMRES_with_fd_UMFPACK_precon"))
    linear_solver=templates_gmres_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson3_linSolver", "templates_BICGSTAB_with_fd_UMFPACK_precon"))
    linear_solver=templates_bicgstab_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson3_linSolver", "templates_CGS_with_fd_UMFPACK_precon"))
    linear_solver=templates_cgs_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson3_linSolver", "LAPACK"))
    linear_solver=LAPACK_dgesv_wrapper;
  else if(Getv("Poisson3_linSolver", "templates_GMRES"))
    linear_solver=templates_gmres_wrapper;
  else if(Getv("Poisson3_linSolver", "templates_BICGSTAB"))
    linear_solver=templates_bicgstab_wrapper;
  else if(Getv("Poisson3_linSolver", "templates_CGS"))
    linear_solver=templates_cgs_wrapper;
  else if(Getv("Poisson3_linSolver", "UMFPACK"))
    linear_solver=UMFPACK_solve_wrapper;
  else if(Getv("Poisson3_linSolver", "templates_GMRES_with_Jacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_Jacobi_precon;
  else if(Getv("Poisson3_linSolver", "templates_GMRES_with_SOR_precon"))
    linear_solver=templates_gmres_wrapper_with_SOR_precon;
  else if(Getv("Poisson3_linSolver", "SOR"))
    linear_solver=SOR_Iterator;
  else if(Getv("Poisson3_linSolver", "templates_GMRES_with_BlockJacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson3_linSolver", "ZIB_GMRES_with_BlockJacobi_precon"))
    linear_solver=ZIB_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson3_linSolver", "ZIB_GBIT_with_BlockJacobi_precon"))
    linear_solver=ZIB_gbit_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson3_linSolver", "ZIB_PCG_with_BlockJacobi_precon"))
    linear_solver=ZIB_pcg_wrapper_with_BlockJacobi_precon;
  else
    errorexit("Poisson3_solve: unknown Poisson3_linSolver");

  /* choose preconditioner */
  if(Getv("Poisson3_linSolver_Precon", "I"))
    Precon=Preconditioner_I;
  else if(Getv("Poisson3_linSolver_Precon", "fd_UMFPACK"))
    Precon=Precon_fd_Poisson3_UMFPACK;
  else if(Getv("Poisson3_linSolver_Precon", "templates"))
    Precon=templates_Preconditioner_for_templates_solver;
  else
    errorexit("Poisson3_solve: unknown Poisson3_linSolver_Precon");
/*
F_Poisson3(VLFu, VLu, VLuAll, VLluAll);
grid->time=0.5;
write_grid(grid);
exit(99);
*/
  /* call Newton solver */
  Newton(F_Poisson3, J_Poisson3, VLu, VLFu, VLuAll, VLluAll,
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
int Poisson3_analyze(tGrid *grid)
{
  int b;
  printf("Poisson3: computing absolute error\n");
  
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
    double *Psi    = box->v[Ind("Poisson3_Psi")];
    double *ErrPsi = box->v[Ind("Poisson3_Err_Psi")];
    double *Chi    = box->v[Ind("Poisson3_Chi")];
    double *ErrChi = box->v[Ind("Poisson3_Err_Chi")];

    /* subtract true values */
    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];
      double r;

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      r = sqrt(x*x + y*y + z*z);

      /* Soln of "Lapl Psi = rh1"  with rh1=-exp(-r*r)/sqrtPI3;
         is "Psi = 1.0/(4*PI*r)*erf(r)" */
      ErrPsi[i] = Psi[i] - 1.0/(4*PI*r)*erf(r);
      ErrChi[i] = Chi[i] - 1.0;
    }
  }
  return 0;
}



/* evaluate Poisson eqn for VLu */
void F_Poisson3(tVarList *VLFu, tVarList *VLu,
                tVarList *VLuAll, tVarList *VLluAll)
{
  tGrid *grid = VLu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[VLFu->index[0]];
    double *Psi  = box->v[VLu->index[0]];
    double *Psix = box->v[VLuAll->index[1]];
    double *Psiy = box->v[VLuAll->index[2]];
    double *Psiz = box->v[VLuAll->index[3]];
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
      FPsi[i] = Chi[i]*Psixx[i] + Psiyy[i] + Chi[i]*Chi[i]*Psizz[i] - rh1[i];
      FChi[i] = Chixx[i] + Psix[i]*Chiyy[i] + Chizz[i] - rh2[i];
FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i] - rh1[i];
FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i] - rh2[i];
    }
  }

  /* BCs */
  Poisson3_set_BCs(VLFu, VLu, VLuAll, 1);
}


/* evaluate linearized eqns */
void J_Poisson3(tVarList *VLJlu, tVarList *VLlu,
                tVarList *VLuAll, tVarList *VLluAll)
{
  tGrid *grid = VLlu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *Psi  = box->v[VLuAll->index[0]];
    double *Psix = box->v[VLuAll->index[1]];
    double *Psiy = box->v[VLuAll->index[2]];
    double *Psiz = box->v[VLuAll->index[3]];
    double *Psixx = box->v[VLuAll->index[4]];
    double *Psiyy = box->v[VLuAll->index[7]];
    double *Psizz = box->v[VLuAll->index[9]];
    double *Chi  = box->v[VLuAll->index[10]];
    double *Chixx = box->v[VLuAll->index[14]];
    double *Chiyy = box->v[VLuAll->index[17]];
    double *Chizz = box->v[VLuAll->index[19]];
    double *JlPsi = box->v[VLJlu->index[0]];
    double *lPsi  = box->v[VLlu->index[0]];
    double *lPsix = box->v[VLluAll->index[1]];
    double *lPsiy = box->v[VLluAll->index[2]];
    double *lPsiz = box->v[VLluAll->index[3]];
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
      JlPsi[i] = Chi[i]*lPsixx[i] + lChi[i]*Psixx[i]  +  lPsiyy[i]  +
                 Chi[i]*Chi[i]*lPsizz[i] + 2.0*lChi[i]*Psizz[i];
      JlChi[i] = lChixx[i] + Psix[i]*lChiyy[i] + lPsix[i]*Chiyy[i] + lChizz[i];
JlPsi[i] = lPsixx[i] + lPsiyy[i] + lPsizz[i];
JlChi[i] = lChixx[i] + lChiyy[i] + lChizz[i];
    }
  }

  /* BCs */
  Poisson3_set_BCs(VLJlu, VLlu, VLluAll, 0);
}



/* switch to fin. diff. and use UMFPACK as lin solver as precon */
void Precon_fd_Poisson3_UMFPACK(tVarList *vlx, tVarList *vlb,
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

  /* when we call J_Poisson3, we evaluate lin. eqns, but now with fd */
  INFO = UMFPACK_solve_wrapper(vlx, vlb, vlr, vlc1, vlc2, 
                               itmax, tol, &normres,
                               J_Poisson3, Preconditioner_I);
  /* restore grid */
  copy_grid_withoutvars(grid_bak, grid, 0);
  free_grid(grid_bak);
  VLDisableFree(vlr);
}
