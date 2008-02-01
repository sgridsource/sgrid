/* BNSdata.c */
/* Wolfgang Tichy 2007 */


#include "sgrid.h"
#include "BNSdata.h"

#define Power pow


/* global var lists */
tVarList *vlu, *vlFu, *vluDerivs;
tVarList *vldu, *vlJdu, *vlduDerivs;


/* functions in this file */
void set_BNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);
void compute_ABphi_from_xyz(tBox *box, double *A, double *B, double *phi,
                            double x, double y, double z);
void make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(tGrid *grid,
     tVarList **vlw,  tVarList **vlwDerivs,  tVarList **vlFw, 
     tVarList **vldw, tVarList **vldwDerivs, tVarList **vlJdw, char *Name);
void free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(
     tVarList *vlw,  tVarList *vlwDerivs,  tVarList *vlFw,
     tVarList *vldw, tVarList *vldwDerivs, tVarList *vlJdw);
int BNS_Eqn_Iterator(tGrid *grid, int itmax, double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);


/* initialize BNSdata */
int BNSdata_startup(tGrid *grid)
{
  int b;
/*
  double A         = Getd("BNSdata_A");
*/
  printf("Initializing BNSdata:\n");

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("BNSdata_Psi",   1, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_Bx",    0, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_alphaP",1, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_Sigma", 0, 1, 1.0);

  /* enable all BNSdata vars */
  enablevar(grid, Ind("BNSdata_Psi"));
  enablevar(grid, Ind("BNSdata_Psix"));
  enablevar(grid, Ind("BNSdata_Psixx"));
  enablevar(grid, Ind("BNSdata_Bx"));
  enablevar(grid, Ind("BNSdata_Bxx"));
  enablevar(grid, Ind("BNSdata_Bxxx"));
  enablevar(grid, Ind("BNSdata_alphaP"));
  enablevar(grid, Ind("BNSdata_alphaPx"));
  enablevar(grid, Ind("BNSdata_alphaPxx"));
  enablevar(grid, Ind("BNSdata_Sigma"));
  enablevar(grid, Ind("BNSdata_Sigmax"));
  enablevar(grid, Ind("BNSdata_Sigmaxx"));
  enablevar(grid, Ind("BNSdata_vRSx"));
  enablevar(grid, Ind("BNSdata_q"));
  enablevar(grid, Ind("BNSdata_vRSxx"));
  enablevar(grid, Ind("BNSdata_qx"));
  enablevar(grid, Ind("BNSdata_temp1"));
  enablevar(grid, Ind("BNSdata_temp2"));
  enablevar(grid, Ind("BNSdata_temp3"));
  enablevar(grid, Ind("BNSdata_temp4"));

  /* enable some lapse and shift of ADMvars */
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));
  
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
    double *Psi    = box->v[Ind("BNSdata_Psi")];
    double *alphaP = box->v[Ind("BNSdata_alphaP")];
    double *q      = box->v[Ind("BNSdata_q")];

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
      /* set Psi and alphaP */
      if(Getv("BNSdata_grid", "SphericalDF"))
      {
        Psi[i] = 1.0;
        alphaP[i] = 1.0;
        q[i]= 1.0;
      }
      else if(Getv("BNSdata_grid", "AnsorgNS") || 
              Getv("BNSdata_grid", "4ABphi_2xyz"))
      {
        double xmax1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,0.0,0.0);
        double xmin1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,1.0,0.0);
        double xmax2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,1.0,0.0);
        double xmin2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,0.0,0.0);
        double R1  = 0.5*(xmax1-xmin1);
        double R2  = 0.5*(xmax2-xmin2);

        Psi[i] = 1.0;
        alphaP[i] = 1.0;
        q[i]= 0.0;

        if(Getv("BNSdata_guess", "test"))
        {
          double xc1 = 0.5*(xmax1+xmin1);
          double xc2 = 0.5*(xmax2+xmin2);

          if(b==0||b==5)
          {
            q[i] = 0.01*( 1.0 - ((x-xc1)*(x-xc1) + y*y + z*z) )/(R1*R1);
          }
          if(b==3||b==4)
          {
            q[i] = 0.01*( 1.0 - ((x-xc2)*(x-xc2) + y*y + z*z) )/(R2*R2);
          }
        }
      }
    }
  }

  return 0;
}


/* Solve the Equations */
int BNSdata_solve(tGrid *grid)
{
  int    Newton_itmax = Geti("BNSdata_Newton_itmax");
  double Newton_tol   = Getd("BNSdata_Newton_tol");
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  tVarList *vldummy;

  /* choose linear solver */
  if(Getv("BNSdata_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("BNSdata_linSolver", "LAPACK"))
    linear_solver=LAPACK_dgesv_wrapper;
  else if(Getv("BNSdata_linSolver", "templates_GMRES"))
    linear_solver=templates_gmres_wrapper;
  else if(Getv("BNSdata_linSolver", "templates_BICGSTAB"))
    linear_solver=templates_bicgstab_wrapper;
  else if(Getv("BNSdata_linSolver", "templates_CGS"))
    linear_solver=templates_cgs_wrapper;
  else if(Getv("BNSdata_linSolver", "UMFPACK"))
    linear_solver=UMFPACK_solve_wrapper;
  else if(Getv("BNSdata_linSolver", "UMFPACK_forSortedVars"))
    linear_solver=UMFPACK_solve_forSortedVars_wrapper;
  else if(Getv("BNSdata_linSolver", "WTsolver"))
    linear_solver=WTsolver;
  else
    errorexit("BNSdata_solve: unknown BNSdata_linSolver");

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vluDerivs= vlalloc(grid);

  /* add all vars to vlu */
  vlpush(vlu, Ind("BNSdata_Psi"));
  vlpush(vlu, Ind("BNSdata_Bx"));
  vlpush(vlu, Ind("BNSdata_alphaP"));
  vlpush(vlu, Ind("BNSdata_Sigma"));

  /* add derivs to vluDerivs */
  vlpush(vluDerivs, Ind("BNSdata_Psix"));
  vlpush(vluDerivs, Ind("BNSdata_Psixx"));
  vlpush(vluDerivs, Ind("BNSdata_Bxx"));
  vlpush(vluDerivs, Ind("BNSdata_Bxxx"));
  vlpush(vluDerivs, Ind("BNSdata_alphaPx"));
  vlpush(vluDerivs, Ind("BNSdata_alphaPxx"));
  vlpush(vluDerivs, Ind("BNSdata_Sigmax"));
  vlpush(vluDerivs, Ind("BNSdata_Sigmaxx"));

  /* enable vlu, vluDerivs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 

  /* now duplicate vlu to get vlFu */  
  vlFu = AddDuplicateEnable(vlu, "_Err");

  /* now duplicate vlFu, vlu and vluDerivs for linarized Eqs. */
  vlJdu      = AddDuplicateEnable(vlFu, "_l");
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlduDerivs = AddDuplicateEnable(vluDerivs, "_l");

  /* How we solve the coupled ell. eqns */
  if(Getv("BNSdata_EllSolver_method", "allatonce"))
  { /* solve the coupled ell. eqns all together */
// remove this later:
//Setd("GridIterators_setABStozero_below", 1e-12); // remove later
//vlFu->n = vlu->n = vlJdu->n = vldu->n = 1;
//Yo(1);
//J_BNSdata(vlJdu, vldu, vlduDerivs, vlu);
Yo(2);
F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
printf("calling write_grid(grid)\n");
write_grid(grid);
//exit(11);

    /* call Newton solver */
    vldummy = vlJdu;
    Newton(F_BNSdata, J_BNSdata, vlu, vlFu, vluDerivs, vldummy,
           Newton_itmax, Newton_tol, &normresnonlin, 1,
           linear_solver, Preconditioner_I, vldu, vlJdu, vlduDerivs, vlu,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
  }
  else if(Getv("BNSdata_EllSolver_method", "sequential"))
  { /* solve the coupled ell. eqns one after an other */
    BNS_Eqn_Iterator(grid, 2, 1e-7, &normresnonlin, linear_solver, 1);
  }
  else
    errorexit("BNSdata_solve: unknown BNSdata_EllSolver_method");

  /* free varlists */     
  VLDisableFree(vldu);
  VLDisableFree(vlduDerivs);
  VLDisableFree(vlJdu);     
  vlfree(vlu);
  vlfree(vluDerivs);
  vlfree(vlFu);
        
  return 0;
}

int setBNSdata(tGrid *grid)
{
  /* call solver */
  BNSdata_solve(grid);

  /* enable all ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("rho"));
  enablevar(grid, Ind("jx"));
  enablevar(grid, Ind("Sxx"));

  /* set all ADM vars */
  setADMvars(grid);

  return 0;
}

/* compute absolute error in ANALYSIS */
int BNSdata_analyze(tGrid *grid)
{
  int b;
/*
  double A         = Getd("BNSdata_A");
*/
  printf("BNSdata: computing absolute error\n");
  
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
    double *Psi    = box->v[Ind("BNSdata_Psi")];
    double *PsiErr = box->v[Ind("BNSdata_Psi_Err")];
    double *alphaP    = box->v[Ind("BNSdata_alphaP")];
    double *alphaPErr = box->v[Ind("BNSdata_alphaP_Err")];

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
      if(Getv("BNSdata_grid", "SphericalDF"))
      {
        PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
        alphaPErr[i] = alphaP[i]-2.0/sqrt(x*x + y*y + z*z);
      }
      else if(Getv("BNSdata_grid", "AnsorgNS") || 
              Getv("BNSdata_grid", "4ABphi_2xyz"))
      {
        double xmax1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,0.0,0.0);
        double xmin1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,1.0,0.0);
        double xmax2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,1.0,0.0);
        double xmin2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,0.0,0.0);
        double R1  = 0.5*(xmax1-xmin1);
        double xc1 = 0.5*(xmax1+xmin1);
        double R2  = 0.5*(xmax2-xmin2);
        double xc2 = 0.5*(xmax2+xmin2);

        if(b==1||b==2||b==3||b==4)
        {
          PsiErr[i] = Psi[i]-1.0/sqrt((x-xc1)*(x-xc1) + y*y + z*z);
        }
        if(b==0||b==5)
        {
          PsiErr[i] = 
            Psi[i]-(-0.5*((x-xc1)*(x-xc1)+y*y+z*z)/(R1*R1*R1)+ 1/R1 + 0.5/R1);
        }
        if(b==1||b==2||b==0||b==5)
        {
          alphaPErr[i] = alphaP[i]-2.0/sqrt((x-xc2)*(x-xc2) + y*y + z*z);
        }
        if(b==3||b==4)
        {
          alphaPErr[i] = 
            alphaP[i]-(-((x-xc2)*(x-xc2)+y*y+z*z)/(R2*R2*R2) + 2/R2 + 1/R2);
        }

      }
    }
  }
  return 0;
}


/* evaluate BNSdata eqns for vlu */
void F_BNSdata(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2)
{
  BNS_CTS(vlFu,vlu,  vlc2,vlc2,vluDerivs, 1);
                   /* ^----^----^--------not used by BNS_CTS if nonlin=1 */
  /* BCs */
  set_BNSdata_BCs(vlFu, vlu, vluDerivs, 1);
}

/* evaluate linearized BNSdata eqns */
void J_BNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  BNS_CTS(vlJdu,vlu,  vlJdu,vldu,vlduDerivs, 0);
        /* ^--not used by BNS_CTS if nonlin=0 */
  /* BCs */
  set_BNSdata_BCs(vlJdu, vldu, vlduDerivs, 0);
}


/* evaluate eqn for a SINGLE one comp. var vlw */
void F_oneComp(tVarList *vlFw, tVarList *vlw,
               tVarList *vlwDerivs, tVarList *vlc2)
{
  /* Note: vlFu,vlu contains vlFw,vlw */
  BNS_CTS(vlFu,vlu,  vlJdu,vldu,vlduDerivs, 1);
                   /* ^-----^----^--------not used by BNS_CTS if nonlin=1 */
  /* BCs */
  set_BNSdata_BCs(vlFw, vlw, vlwDerivs, 1);
}

/* evaluate linearized eqn for a SINGLE one comp. var for vldw */
void J_oneComp(tVarList *vlJdw, tVarList *vldw,
               tVarList *vldwDerivs, tVarList *vlw)
{
  /* Note: vlJdu,vldu contains vlJdw,vldw */
  BNS_CTS(vlFu,vlu,  vlJdu,vldu,vlduDerivs, 0);
        /* ^--not used by BNS_CTS if nonlin=0 */
  /* BCs */
  set_BNSdata_BCs(vlJdw, vldw, vldwDerivs, 0);
}


/* set BCs for a varlist */
/* NOTE: this works only for a varlist made up of scalars!!!
         because to compute the varlist index of the derivs we
         use stuff like vind*9 */
void set_BNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  int vind;
  int vindDerivs=0;

  if( grid->box[0]->n1 != grid->box[1]->n1 ||
      grid->box[3]->n1 != grid->box[2]->n1 ||
      grid->box[1]->n1 != grid->box[2]->n1 ) 
    errorexit("all n1 in boxes0-3 must be the same because we currently use "
              "lines like:\n"
              "FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];\n"
              "where Psi = grid->box[0]->v[vi], P = grid->box[1]->v[vi]");
                      
  for(vind=0; vind<vlu->n; vind++)
  {
    int ncomp = VarNComponents(vlu->index[vind]);
    double PsiFarLimit = VarFarLimit(vlu->index[vind])*nonlin;

    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[vlFu->index[vind]];
      double *Psi  = box->v[vlu->index[vind]];
      double *Psix = box->v[vluDerivs->index[vindDerivs]];
      double *Psiy = box->v[vluDerivs->index[vindDerivs+1]];
      double *Psiz = box->v[vluDerivs->index[vindDerivs+2]];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;

      /* BCs */
      if(Getv("BNSdata_grid", "SphericalDF"))
      {
        forplane1(i,j,k, n1,n2,n3, 0)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*(vind+1)*nonlin;

        forplane1(i,j,k, n1,n2,n3, n1-1)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*(vind+1)*nonlin;
      }
      else if (Getv("BNSdata_grid", "AnsorgNS"))
      {
        double *P;
        double *dP[4];
        double *BM;

        /* special rho=0 case??? */
        if(b==0 || b==1 || b==2 || b==3)
        {
          int pl;
          char str[1000];
          snprintf(str, 999, "box%d_basis2", b);
          if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* loop over rho=0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                  FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                {
                  double Psi_rho = Psiy[Index(i,j,k)];
                  double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                  FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                }
              }
          }
          /* same as before, but also interpolate to rho=0 */
          else if(Getv("BNSdata_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
          /* same as before, but do it only in box0/3 at A,B=1,0 and A,B=1,1 */
          else if(Getv("BNSdata_regularization", 
                       "regularity_on_axis_at_center") && (b==0 || b==3) )
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              i=n1-1;   /* do it only at A=1 */
              j=pl;
              for(k=0; k<n3; k++)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
        } /* end: special rho=0 case??? */

        if(b==0)  /* in box0 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box0 and box1 */
          P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box0 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==3)  /* in box3 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box3 and box2 */
          P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box3 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==1)
        {
          BM = (double *) calloc(max3(grid->box[0]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[0], 1, BM, 0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box0 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[0]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box1_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else if(b==2)
        {
          BM = (double *) calloc(max3(grid->box[3]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* values at A=1 are equal in box1 and box2 */
          P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[3], 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box3 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[3]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box2_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else errorexiti("b=%d should be impossible!", b);
      } /* end: else if (Getv("BNSdata_grid", "AnsorgNS")) */
      else if (Getv("BNSdata_grid", "4ABphi_2xyz"))
      {
        double *P;
        double *dP[4];
        double *X, *Y, *Z;
        double *Pcoeffs;
        double Pinterp;
        double x,y,z;

        /* special rho=0 case??? */
        if(b==0 || b==1 || b==2 || b==3)
        {
          int pl;
          char str[1000];
          snprintf(str, 999, "box%d_basis2", b);
          if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* loop over rho=0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_ijk = u_ij0 (not u_phi_phi=0) */
                  FPsi[Index(i,j,k)] = Psi[Index(i,j,k)]-Psi[Index(i,j,0)];
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                {
                  double Psi_rho = Psiy[Index(i,j,k)];
                  double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                  FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                }
              }
          }
          /* same as before, but also interpolate to rho=0 */
          else if(Getv("BNSdata_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
        } /* end: special rho=0 case??? */

        if(b==0)  /* in box0 */
        {
          /* values at A=0 are equal in box0 and box1 */
          P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* values at A=Amin are interpolated from box5 */
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[5]->v[vlu->index[vind]]; /* values in box5 */
          Pcoeffs = grid->box[5]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[5], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
          {
            int ind=Index(i,j,k);
            x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            Pinterp = spec_interpolate(grid->box[5], Pcoeffs, x,y,z);
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==3)  /* in box3 */
        {
          /* values at A=0 are equal in box3 and box2 */
          P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* values at A=Amin are interpolated from box4 */
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[4]->v[vlu->index[vind]]; /* values in box4 */
          Pcoeffs = grid->box[4]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[4], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
          {
            int ind=Index(i,j,k);
            x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
            Pinterp = spec_interpolate(grid->box[4], Pcoeffs, x,y,z);
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==1)  /* in box1 */
        {
          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }
          
          /* Psi=0 at infinity */
          if(Getv("box1_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *BM = (double *) calloc(n2, sizeof(double));
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
            free(BM);
          }
        }
        else if(b==2)  /* in box2 */
        {
          /* values at A=1 are equal in box1 and box2 */
          P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];

          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }

          /* Psi=0 at infinity */
          if(Getv("box2_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *BM = (double *) calloc(n2, sizeof(double));
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
            free(BM);
          }
        }
        else if(b==5)  /* in box5 */
        {
          /* values at border are interpolated from box0 */
          double A,B,phi;
          int pl, k_phi;
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[0]->v[vlu->index[vind]]; /* values in box0 */
          Pcoeffs = grid->box[0]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[0], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int i0;
            int ind=Index(pl,0,0);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);

              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int i0;
            int ind=Index(0,pl,0);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int i0;
            int ind=Index(0,0,pl);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else if(b==4)  /* in box4 */
        {
          /* values at border are interpolated from box3 */
          double A,B,phi;
          int pl, k_phi;
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[3]->v[vlu->index[vind]]; /* values in box3 */
          Pcoeffs = grid->box[3]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[3], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int i0;
            int ind=Index(pl,0,0);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int i0;
            int ind=Index(0,pl,0);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int i0;
            int ind=Index(0,0,pl);
            phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
                                      X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else errorexiti("b=%d should be impossible!", b);
      } /* end: else if (Getv("BNSdata_grid", "4ABphi_2xyz")) */

    } /* end forallboxes */
    /* increase index for derivs */
    vindDerivs += 3;
    if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}

/* compute A,B,phi from x,y,z */
void compute_ABphi_from_xyz(tBox *box, double *A, double *B, double *phi,
                            double x, double y, double z)
{
  double ph=Arg(y,z);
  if(ph<0) ph = 2.0*PI + ph;
  if(*A==0.0) *A+=1e-10;
  if(*A==1.0) *A-=1e-10;
  if(*B==0.0) *B+=1e-10;
  if(*B==1.0) *B-=1e-10;
  *phi = ph;
  XYZ_of_xyz(box, A,B,phi, x,y,z);
  *phi = ph;
  if(*A<0.0) *A=0.0;
  if(*A>1.0) *A=1.0;
  if(*B<0.0) *B=0.0;
  if(*B>1.0) *B=1.0;
}

/* make var lists that contain VarComp Name, its derivs, its errors,
   and the linearized var and its derivs and errors */
void make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(tGrid *grid,
     tVarList **vlw,  tVarList **vlwDerivs,  tVarList **vlFw, 
     tVarList **vldw, tVarList **vldwDerivs, tVarList **vlJdw, char *Name)
{
  char *str;
  
  str = (char *) calloc(strlen(Name)+20, sizeof(char) );

  /* allocate varlists */
  *vlw       = vlalloc(grid);
  *vlwDerivs = vlalloc(grid);
  *vlFw      = vlalloc(grid);
  *vldw      = vlalloc(grid);
  *vldwDerivs= vlalloc(grid);
  *vlJdw     = vlalloc(grid);

  /* add Name to vlw, ... */
  sprintf(str, "%s", Name);        vlpushone(*vlw,       Ind(str));
  sprintf(str, "%sx", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sy", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sz", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxx", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxy", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%syy", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%syz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%szz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%s_Err", Name);    vlpushone(*vlFw,      Ind(str));
  sprintf(str, "%s_l", Name);      vlpushone(*vldw,       Ind(str));
  sprintf(str, "%sx_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sy_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sz_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxx_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxy_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%syy_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%syz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%szz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%s_Err_l", Name);  vlpushone(*vlJdw,      Ind(str));
  free(str);
}

/* the var lists vlw, vlwDerivs, vlFw, vldw, vldwDerivs, vlJdw */
void free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(
     tVarList *vlw,  tVarList *vlwDerivs,  tVarList *vlFw,
     tVarList *vldw, tVarList *vldwDerivs, tVarList *vlJdw)
{
  vlfree(vlw);
  vlfree(vlwDerivs);
  vlfree(vlFw);
  vlfree(vldw);
  vlfree(vldwDerivs);
  vlfree(vlJdw);
}          

/* solve the coupled ell. eqns one after an other, and iterate */
int BNS_Eqn_Iterator(tGrid *grid, int itmax, double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  int    Newton_itmax = Geti("BNSdata_Newton_itmax");
  double Newton_tol   = Getd("BNSdata_Newton_tol");
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int it;
  int prN=0;

  if(pr) printf("BNS_Eqn_Iterator:  starting iterations ...\n");

  for (it = 0; it < itmax; it++)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* make new vlw, ... for Psi */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Psi");
    /* call Newton solver for Psi */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_Psi:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Bx */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Bx");
    /* call Newton solver for Bx */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_Bx:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for By */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_By");
    /* call Newton solver for By */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_By:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Bz */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Bz");
    /* call Newton solver for Bz */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_Bz:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for alphaP */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_alphaP");
    /* call Newton solver for alphaP */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_alphaP:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Sigma */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Sigma");
    /* call Newton solver for Sigma */
    prdivider(1);
    printf("Solving elliptic Eqn for BNSdata_Sigma:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* compute error vlFu = F(u) */
    F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
    *normres = GridL2Norm(vlFu);
    if(pr)
    {
      printf("BNS_Eqn_Iterator step %d residual = %.4e\n", it, *normres);
      fflush(stdout);
    }
    if (*normres <= tol) break;
  }
  /* warn if we didn't converge */
  if (it >= itmax)
  {
    printf("BNS_Eqn_Iterator: *** Too many steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
  }
  return it;
}
