/* Poisson.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "Poisson.h"

#define Power pow



void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);
void ABphi_of_xyz(tBox *box, double *A, double *B, double *phi,
                  double x, double y, double z);
void convert_grid_to_fd_old(tGrid *grid);
void Precon_fd_Poisson_UMFPACK(tVarList *vlx, tVarList *vlb,
                               tVarList *vlc1, tVarList *vlc2);


/* initialize Poisson */
int Poisson_startup(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
*/
  printf("Initializing Poisson:\n");

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("Poisson_Psi", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Poisson_Chi", 0, 1, 1.0);

  /* enable all vars */
  enablevar(grid, Ind("Poisson_Psi"));
  enablevar(grid, Ind("Poisson_rh1"));
  enablevar(grid, Ind("Poisson_Psix"));
  enablevar(grid, Ind("Poisson_Psixx"));
  enablevar(grid, Ind("Poisson_Err_Psi"));
  enablevar(grid, Ind("Poisson_Chi"));
  enablevar(grid, Ind("Poisson_rh2"));
  enablevar(grid, Ind("Poisson_Chix"));
  enablevar(grid, Ind("Poisson_Chixx"));
  enablevar(grid, Ind("Poisson_Err_Chi"));
  enablevar(grid, Ind("Poisson_temp1"));
  enablevar(grid, Ind("Poisson_temp2"));
  enablevar(grid, Ind("Poisson_temp3"));
  enablevar(grid, Ind("Poisson_temp4"));
  enablevar(grid, Ind("Poisson_temp5"));
  enablevar(grid, Ind("Poisson_temp6"));
  
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
    double *Psi = box->v[Ind("Poisson_Psi")];
    double *Chi = box->v[Ind("Poisson_Chi")];
    double *rh1 = box->v[Ind("Poisson_rh1")];
    double *rh2 = box->v[Ind("Poisson_rh2")];

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
      if(Getv("Poisson_grid", "SphericalDF"))
      {
        Psi[i] = x*y*z; // 1.0/sqrt(x*x + y*y + z*z);
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
      }
      else if(Getv("Poisson_grid", "AnsorgNS") || 
              Getv("Poisson_grid", "4ABphi_2xyz"))
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

        Psi[i] = 0.0;
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
        if(b==0 || b==5)  rh1[i] = -3.0/(R1*R1*R1);
        if(b==3 || b==4)  rh2[i] = -6.0/(R2*R2*R2);

        if(Getv("Poisson_guess", "exact"))
        {
          double xc1 = 0.5*(xmax1+xmin1);
          double xc2 = 0.5*(xmax2+xmin2);

          if(b==1||b==2||b==3||b==4)
          {
            Psi[i] = 1.0/sqrt((x-xc1)*(x-xc1) + y*y + z*z);
          }
          if(b==0||b==5)
          {
            Psi[i] = 
              (-0.5*((x-xc1)*(x-xc1)+y*y+z*z)/(R1*R1*R1)+ 1/R1 + 0.5/R1);
          }
          if(b==1||b==2||b==0||b==5)
          {
            Chi[i] = 2.0/sqrt((x-xc2)*(x-xc2) + y*y + z*z);
          }
          if(b==3||b==4)
          {
            Chi[i] = 
              (-((x-xc2)*(x-xc2)+y*y+z*z)/(R2*R2*R2) + 2/R2 + 1/R2);
          }
        }

/*
        Psi[i] = b*b;
        Psi[i] = 1*exp(-0.5*(
                 (x-1)*(x-1)/(1*1) + (y-0)*(y-0)/(1.5*1.5)
               + (z-0.5)*(z-0.5)/(2*2)   ) );
*/
        //Psi[i] = 1.0/sqrt(x*x + y*y + z*z+1);
        //Psi[i] = 1.0/sqrt(x*x + (y-1.5)*(y-1.5) + (z-0.2)*(z-0.2)+1);
      }
    }
  }

  return 0;
}


/* Solve the Equations */
int Poisson_solve(tGrid *grid)
{
  int    itmax = Geti("Poisson_itmax");
  double tol   = Getd("Poisson_tol");
  int    linSolver_itmax  = Geti("Poisson_linSolver_itmax");
  double linSolver_tolFac = Getd("Poisson_linSolver_tolFac");
  double linSolver_tol    = Getd("Poisson_linSolver_tol");
  double normresnonlin;
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  void (*Precon)(tVarList *vlx, tVarList *vlb, tVarList *vlc1, tVarList *vlc2);
  tVarList *vlu, *vlFu, *vluDerivs;
  tVarList *vldu, *vlr, *vlduDerivs;
  tVarList *vlrhs;

  /* Set GridIterators_setABStozero_below to 1e-12 for compatibility */
  Setd("GridIterators_setABStozero_below", 1e-12); // remove later

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vlFu = vlalloc(grid);
  vluDerivs= vlalloc(grid);
  vlrhs = vlalloc(grid);

  /* add Poisson_Psi and Poisson_Chi to vlu */
  vlpush(vlu, Ind("Poisson_Psi"));
  vlpush(vlu, Ind("Poisson_Chi"));

  /* add derivs of Poisson_Psi and Poisson_Chi to vluDerivs */
  vlpush(vluDerivs, Ind("Poisson_Psix"));
  vlpush(vluDerivs, Ind("Poisson_Psixx"));
  vlpush(vluDerivs, Ind("Poisson_Chix"));
  vlpush(vluDerivs, Ind("Poisson_Chixx"));

  /* add Poisson_rh1 and Poisson_rh2 to vlrhs */
  vlpush(vlrhs, Ind("Poisson_rh1"));
  vlpush(vlrhs, Ind("Poisson_rh2"));

  /* add Poisson_Err_Psi and Poisson_Err_Chi to vlFu */
  vlpush(vlFu, Ind("Poisson_Err_Psi"));
  vlpush(vlFu, Ind("Poisson_Err_Chi"));
    
  /* enable vlu, vluDerivs, vlrhs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 
  enablevarlist(vlrhs);
  enablevarlist(vlFu);

  /* now duplicate vlu and vluDerivs for linarized Eqs. */
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlr        = AddDuplicateEnable(vlFu, "_l");
  vlduDerivs = vluDerivs; /* maybe: vlduDerivs=AddDuplicateEnable(vluDerivs, "_l"); */

  /* choose linear solver */
  if(Getv("Poisson_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("Poisson_linSolver", "bicgstab_with_fd_UMFPACK_precon"))
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
  else if(Getv("Poisson_linSolver", "templates_GMRES_with_fd_UMFPACK_precon"))
    linear_solver=templates_gmres_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson_linSolver", "templates_BICGSTAB_with_fd_UMFPACK_precon"))
    linear_solver=templates_bicgstab_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson_linSolver", "templates_CGS_with_fd_UMFPACK_precon"))
    linear_solver=templates_cgs_wrapper_with_fd_UMFPACK_precon;
  else if(Getv("Poisson_linSolver", "LAPACK"))
    linear_solver=LAPACK_dgesv_wrapper;
  else if(Getv("Poisson_linSolver", "templates_GMRES"))
    linear_solver=templates_gmres_wrapper;
  else if(Getv("Poisson_linSolver", "templates_BICGSTAB"))
    linear_solver=templates_bicgstab_wrapper;
  else if(Getv("Poisson_linSolver", "templates_CGS"))
    linear_solver=templates_cgs_wrapper;
  else if(Getv("Poisson_linSolver", "UMFPACK"))
    linear_solver=UMFPACK_solve_wrapper;
  else if(Getv("Poisson_linSolver", "templates_GMRES_with_Jacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_Jacobi_precon;
  else if(Getv("Poisson_linSolver", "templates_GMRES_with_SOR_precon"))
    linear_solver=templates_gmres_wrapper_with_SOR_precon;
  else if(Getv("Poisson_linSolver", "SOR"))
    linear_solver=SOR_Iterator;
  else if(Getv("Poisson_linSolver", "templates_GMRES_with_BlockJacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson_linSolver", "ZIB_GMRES_with_BlockJacobi_precon"))
    linear_solver=ZIB_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson_linSolver", "ZIB_GBIT_with_BlockJacobi_precon"))
    linear_solver=ZIB_gbit_wrapper_with_BlockJacobi_precon;
  else if(Getv("Poisson_linSolver", "ZIB_PCG_with_BlockJacobi_precon"))
    linear_solver=ZIB_pcg_wrapper_with_BlockJacobi_precon;
  else
    errorexit("Poisson_solve: unknown Poisson_linSolver");

  /* choose preconditioner */
  if(Getv("Poisson_linSolver_Precon", "I"))
    Precon=Preconditioner_I;
  else if(Getv("Poisson_linSolver_Precon", "fd_UMFPACK"))
    Precon=Precon_fd_Poisson_UMFPACK;
  else if(Getv("Poisson_linSolver_Precon", "templates"))
    Precon=templates_Preconditioner_for_templates_solver;
  else
    errorexit("Poisson_solve: unknown Poisson_linSolver_Precon");

// remove this later:
tGrid *grid_bak=make_empty_grid(grid->nvariables, 1);
copy_grid_withoutvars(grid, grid_bak, 1);
convert_grid_to_fd(grid);
int b;
for(b=0;b<grid->nboxes;b++)
{
printf("box%d D\n",b);
printmatrix(grid->box[b]->D1, grid->box[b]->n1);
printmatrix(grid->box[b]->D2, grid->box[b]->n2);
printmatrix(grid->box[b]->D3, grid->box[b]->n3);
printf("box%d Mcoeffs\n",b);
printmatrix(grid->box[b]->Mcoeffs1, grid->box[b]->n1);
printmatrix(grid->box[b]->Mcoeffs2, grid->box[b]->n2);
printmatrix(grid->box[b]->Mcoeffs3, grid->box[b]->n3);
printf("box%d Meval\n",b);
printmatrix(grid->box[b]->Meval1, grid->box[b]->n1);
printmatrix(grid->box[b]->Meval2, grid->box[b]->n2);
printmatrix(grid->box[b]->Meval3, grid->box[b]->n3);
}
//printmatrix(grid->box[0]->D3,     grid->box[0]->n3);
//printmatrix(grid->box[1]->D1, 	  grid->box[1]->n1);
//printmatrix(grid->box[2]->Mcoeffs2,     grid->box[2]->n2);
//printmatrix(grid->box[4]->D3, 	  grid->box[4]->n3);
//printmatrix(grid->box[3]->Meval3,     grid->box[3]->n3);

copy_grid_withoutvars(grid_bak, grid, 1);
free_grid(grid_bak);


//vlFu->n = vlu->n = vlr->n = vldu->n = 1;

F_Poisson(vlFu, vlu, vluDerivs, vlrhs);
printf("calling write_grid(grid)\n");
write_grid(grid);
//exit(11);

  /* call Newton solver */
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, vlrhs,
         itmax, tol, &normresnonlin, 1,
         linear_solver, Precon, vldu, vlr, vlduDerivs, vlrhs,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);

  /* free varlists */     
  VLDisableFree(vldu);
  /* VLDisableFree(vlduDerivs); */
  VLDisableFree(vlr);     
  vlfree(vlu);
  vlfree(vluDerivs);
  vlfree(vlrhs);
  vlfree(vlFu);
        
  return 0;
}


/* compute absolute error in ANALYSIS */
int Poisson_analyze(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
*/
  printf("Poisson: computing absolute error\n");
  
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
    double *Psi    = box->v[Ind("Poisson_Psi")];
    double *PsiErr = box->v[Ind("Poisson_Err_Psi")];
    double *Chi    = box->v[Ind("Poisson_Chi")];
    double *ChiErr = box->v[Ind("Poisson_Err_Chi")];


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
      if(Getv("Poisson_grid", "SphericalDF"))
      {
        //PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
        //ChiErr[i] = Chi[i]-2.0/sqrt(x*x + y*y + z*z);
      }
      else if(Getv("Poisson_grid", "AnsorgNS") || 
              Getv("Poisson_grid", "4ABphi_2xyz"))
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

//if(i==0) printf("R1=%g xc1=%g\n", R1, xc1);

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
          ChiErr[i] = Chi[i]-2.0/sqrt((x-xc2)*(x-xc2) + y*y + z*z);
        }
        if(b==3||b==4)
        {
          ChiErr[i] = 
            Chi[i]-(-((x-xc2)*(x-xc2)+y*y+z*z)/(R2*R2*R2) + 2/R2 + 1/R2);
        }

      }
    }
  }
  return 0;
}


void filterSing(tBox *box, double *u)
{
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  int filt2=n2-1;
  int filt3=n3-1; // n3-1; // n3/2;
  int i,j,k;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int N, double X)=NULL;
  static int linelen=0;
  static double *uline=NULL;
  static double *ufline=NULL;
  int m3;

  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints, 
                            &filter_coeffs, &basisfunc);
  initfiltermatrix(box->F2, n2+1-filt2, n2, 
                   get_coeffs, filter_coeffs, eval_onPoints);

  get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints, 
                            &filter_coeffs, &basisfunc);
  initfiltermatrix(box->F3, n3+1-filt3, n3, 
                   get_coeffs, filter_coeffs, eval_onPoints);

  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    ufline = (double*) realloc(ufline, linelen * sizeof(double));
  }

  i=box->n1-1;
//  /* filter along phi-direc if A=1 for B=Bmin and B=Bmax */
//  for(j=0;j<n2;j=j+n2-1)
  /* filter along phi-direc if A=1 for all B */
  for(j=0;j<n2;j++)
  {
    get_memline(u, uline,  3, i, j, n1, n2, n3);
    matrix_times_vector(box->F3 , uline, ufline, box->n3);
    put_memline(u, ufline, 3, i, j, n1, n2, n3);
  }
  /* filter along B-direc if A=1 for all phi */
  for(k=0;k<n3;k++)
  {
    get_memline(u, uline,  2, i, k, n1, n2, n3);
    matrix_times_vector(box->F2, uline, ufline, box->n2);
    put_memline(u, ufline, 2, i, k, n1, n2, n3);
  }
  /* filter along phi-direc if A=1 for B=Bmin and B=Bmax */
  //  for(j=0;j<n2;j=j+n2-1)

}                




/* evaluate Poisson eqn for vlu */
void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlrhs)
{
  tGrid *grid = vlu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[vlFu->index[0]];
    // double *Psi  = box->v[vlu->index[0]];
    double *rh1  = box->v[vlrhs->index[0]];
    double *Psixx = box->v[vluDerivs->index[3]];
    double *Psiyy = box->v[vluDerivs->index[6]];
    double *Psizz = box->v[vluDerivs->index[8]];
    double *FChi = box->v[vlFu->index[1]];
    // double *Chi  = box->v[vlu->index[1]];
    double *rh2  = box->v[vlrhs->index[1]];
    double *Chixx = box->v[vluDerivs->index[12]];
    double *Chiyy = box->v[vluDerivs->index[15]];
    double *Chizz = box->v[vluDerivs->index[17]];
    int i;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes") &&
       (box->ddX_dxdx[1][1][1]!=NULL || box->x_of_X[1]==NULL))
    {
      allDerivsOf_S(box, vlu->index[0], vluDerivs->index[0],
                    vluDerivs->index[3]);
      allDerivsOf_S(box, vlu->index[1], vluDerivs->index[9],
                    vluDerivs->index[12]);
    }
    else
    {
      FirstAndSecondDerivsOf_S(box, vlu->index[0], vluDerivs->index[0],
                               vluDerivs->index[3]);
      FirstAndSecondDerivsOf_S(box, vlu->index[1], vluDerivs->index[9],
                               vluDerivs->index[12]);
    }
    
    if(b==0 & 0)
    {
      double *Psix = box->v[vluDerivs->index[0]];
      double *Psiy = box->v[vluDerivs->index[1]];
      double *Psiz = box->v[vluDerivs->index[2]];
      double *Psixy = box->v[vluDerivs->index[4]];
      double *Psixz = box->v[vluDerivs->index[5]];
      double *Psiyz = box->v[vluDerivs->index[7]];
      // //FirstDerivsOf_S(box, vlu->index[0], vluDerivs->index[0]);
      filterSing(box, Psix); // ???
      filterSing(box, Psiy);
      filterSing(box, Psiz);
      FirstDerivsOf_S(box, vluDerivs->index[2],vluDerivs->index[6]);
      FirstDerivsOf_S(box, vluDerivs->index[1],vluDerivs->index[5]);
      FirstDerivsOf_S(box, vluDerivs->index[0],vluDerivs->index[3]);
      //filterSing(box, Psixx);
      //filterSing(box, Psiyy);
      //filterSing(box, Psizz);
      ////filterSing(box, Psiyz);
    }

//    forallpoints(box, i)  FPsi[i] = Psi[i]*Psi[i] - 9.0;
    forallpoints(box, i)
    {
      FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i] - rh1[i];
      FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i] - rh2[i];
    }

/*
double *Psi  = box->v[vlu->index[0]];
double *Psix = box->v[vluDerivs->index[0]];
double *Psiy = box->v[vluDerivs->index[1]];
double *Psiz = box->v[vluDerivs->index[2]];
double *Psixy = box->v[vluDerivs->index[4]];
double *Psixz = box->v[vluDerivs->index[5]];
double *Psiyz = box->v[vluDerivs->index[7]];

double *M = box->v[vluDerivs->index[4]];
double *c = box->v[vluDerivs->index[5]];
double *d = box->v[vluDerivs->index[7]];

spec_Coeffs(box, Psi, c);

double X,Y,Z;
double x,y,z;


if(b<4)
{
x=1.2; y=-0.5; z=0.5;
printf("nearestXYZ=%f ", nearestXYZ_of_xyz(box, &i, &X,&Y,&Z, x,y,z));
printf("(X,Y,Z)=(%.12f,%.12f,%.12f)\n", X,Y,Z);
//X=0.1; Y=0.3; 
Z=PI-PI/4;
XYZ_of_xyz(box, &X,&Y,&Z, x,y,z);
printf("(x,y,z)=(%f,%f,%f)   (X,Y,Z)=(%.12f,%.12f,%.12f)\n", x,y,z, X,Y,Z);
printf("##### Psi=%.12f\n", spec_interpolate(box, c, X,Y,Z));
}
*/
  }

  /* BCs */
  set_BCs(vlFu, vlu, vluDerivs, 1);
}


/* evaluate linearized Laplace */
void J_Poisson(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  tGrid *grid = vldu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *JdPsi = box->v[vlJdu->index[0]];
    // double *dPsi  = box->v[vldu->index[0]];
    double *dPsixx = box->v[vlduDerivs->index[3]];
    double *dPsiyy = box->v[vlduDerivs->index[6]];
    double *dPsizz = box->v[vlduDerivs->index[8]];
    double *JdChi = box->v[vlJdu->index[1]];
    // double *dChi  = box->v[vldu->index[1]];
    double *dChixx = box->v[vlduDerivs->index[12]];
    double *dChiyy = box->v[vlduDerivs->index[15]];
    double *dChizz = box->v[vlduDerivs->index[17]];
    int i;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes") && 
       (box->ddX_dxdx[1][1][1]!=NULL || box->x_of_X[1]==NULL))
    {
      allDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                    vlduDerivs->index[3]);
      allDerivsOf_S(box, vldu->index[1], vlduDerivs->index[9],
                    vlduDerivs->index[12]);
    }
    else
    {
      FirstAndSecondDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                               vlduDerivs->index[3]);
      FirstAndSecondDerivsOf_S(box, vldu->index[1], vlduDerivs->index[9],
                               vlduDerivs->index[12]);
    }
//    forallpoints(box, i)  JdPsi[i] = 2.0*Psi[i]*dPsi[i];
    forallpoints(box, i)
    {
      JdPsi[i] = dPsixx[i] + dPsiyy[i] + dPsizz[i];
      JdChi[i] = dChixx[i] + dChiyy[i] + dChizz[i];
    }
  }

  /* BCs */
  set_BCs(vlJdu, vldu, vlduDerivs, 0);
}

/* switch to fin. diff. and use UMFPACK as lin solver as precon */
void Precon_fd_Poisson_UMFPACK(tVarList *vlx, tVarList *vlb,
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

  /* when we call J_Poisson, we evaluate lin. eqns, but now with fd */
  INFO = UMFPACK_solve_wrapper(vlx, vlb, vlr, vlc1, vlc2, 
                               itmax, tol, &normres,
                               J_Poisson, Preconditioner_I);
  /* restore grid */
  copy_grid_withoutvars(grid_bak, grid, 0);
  free_grid(grid_bak);
  VLDisableFree(vlr);
}

/* set BCs for a varlist */
void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  int vind;

  for(vind=0; vind<vlu->n; vind++)
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[vlFu->index[vind]];
      double *Psi  = box->v[vlu->index[vind]];
      double *Psix = box->v[vluDerivs->index[9*vind]];
      double *Psiy = box->v[vluDerivs->index[9*vind+1]];
      double *Psiz = box->v[vluDerivs->index[9*vind+2]];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;

      /* BCs */
      if(Getv("Poisson_grid", "SphericalDF"))
      {
        forplane1(i,j,k, n1,n2,n3, 0)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*(vind+1)*nonlin;

        forplane1(i,j,k, n1,n2,n3, n1-1)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*(vind+1)*nonlin;
      }
      else if (Getv("Poisson_grid", "AnsorgNS"))
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
            double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("Poisson_temp2")];
            double *temp3 = box->v[Ind("Poisson_temp3")];
            double *temp4 = box->v[Ind("Poisson_temp4")];

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
          else if(Getv("Poisson_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("Poisson_temp2")];
            double *temp3 = box->v[Ind("Poisson_temp3")];
            double *temp4 = box->v[Ind("Poisson_temp4")];
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
          else if(Getv("Poisson_regularization", 
                       "regularity_on_axis_at_center") && (b==0 || b==3) )
          {
            double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("Poisson_temp2")];
            double *temp3 = box->v[Ind("Poisson_temp3")];
            double *temp4 = box->v[Ind("Poisson_temp4")];
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
          dP[1] = grid->box[2]->v[vluDerivs->index[vind*9]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vind*9]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vind*9+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vind*9+2]];
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
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          else // fix: B=0 is not on grid for ChebZeros!!!
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
              FPsi[Index(n1-1,0,k)] = U0;
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
          dP[1] = grid->box[3]->v[vluDerivs->index[vind*9]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vind*9+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vind*9+2]];
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
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          else // fix: B=0 is not on grid for ChebZeros!!!
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
              FPsi[Index(n1-1,0,k)] = U0;
            }
            free(line);
          }

          free(BM);
        }
        else errorexiti("b=%d should be impossible!", b);
      } /* end: else if (Getv("Poisson_grid", "AnsorgNS")) */
      else if (Getv("Poisson_grid", "4ABphi_2xyz"))
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
            double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("Poisson_temp2")];
            double *temp3 = box->v[Ind("Poisson_temp3")];
            double *temp4 = box->v[Ind("Poisson_temp4")];

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
          else if(Getv("Poisson_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("Poisson_temp2")];
            double *temp3 = box->v[Ind("Poisson_temp3")];
            double *temp4 = box->v[Ind("Poisson_temp4")];
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
          Pcoeffs = grid->box[5]->v[Ind("Poisson_temp1")];
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
          Pcoeffs = grid->box[4]->v[Ind("Poisson_temp1")];
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
          dP[1] = grid->box[2]->v[vluDerivs->index[vind*9]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vind*9]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vind*9+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vind*9+2]];
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
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          else // fix: B=0 is not on grid for ChebZeros!!!
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
              FPsi[Index(n1-1,0,k)] = U0;
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
          dP[1] = grid->box[3]->v[vluDerivs->index[vind*9]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vind*9+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vind*9+2]];
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
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          else // fix: B=0 is not on grid for ChebZeros!!!
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
              FPsi[Index(n1-1,0,k)] = U0;
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
          Pcoeffs = grid->box[0]->v[Ind("Poisson_temp1")];
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
              ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);

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
              ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
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
              ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
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
          Pcoeffs = grid->box[3]->v[Ind("Poisson_temp1")];
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
              ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
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
              ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
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
              ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                             
              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else errorexiti("b=%d should be impossible!", b);

      } /* end: else if (Getv("Poisson_grid", "4ABphi_2xyz")) */

    } /* end forallboxes */
}

/* compute A,B,phi from x,y,z */
void ABphi_of_xyz(tBox *box, double *A, double *B, double *phi,
                  double x, double y, double z)
{
  if(*A==0.0) *A+=1e-10;
  if(*A==1.0) *A-=1e-10;
  if(*B==0.0) *B+=1e-10;
  if(*B==1.0) *B-=1e-10;
  *phi = Arg(y,z);   if(*phi<0) *phi = 2.0*PI + *phi;
  XYZ_of_xyz(box, A,B,phi, x,y,z);
  if(*A<0.0) *A=0.0;
  if(*A>1.0) *A=1.0;
  if(*B<0.0) *B=0.0;
  if(*B>1.0) *B=1.0;
}

/* reset matrices and basis funcs to do finite differencing */
void convert_grid_to_fd_old(tGrid *grid)
{
  int b, bmax;

  convert_grid_to_fd(grid);
  /* Until June 7 2012, the above line was calling the code that is now in: 
     convert_grid_to_fd_matrices_onesidedBC(grid); */

  /* make use of the fact that the phi-direction is periodic */
  if(Getv("Poisson_grid", "4ABphi_2xyz")) bmax=4;
  else                                    bmax=grid->nboxes;
  for(b=0; b<bmax; b++)
  {
    tBox *box = grid->box[b];
    int n3 = box->n3;

    initdiffmatrix(box, 3, box->D3, box->DD3, n3,
                   fd2_coeffs, fd2_deriv_periodic, fd2_eval);
  }      
}


/**************************************************************/
/* old stuff: */

/* set BCs in the old way for Psi and Chi */
void set_BCs_old(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[vlFu->index[0]];
    double *Psi  = box->v[vlu->index[0]];
    double *Psix = box->v[vluDerivs->index[0]];
    double *Psiy = box->v[vluDerivs->index[1]];
    double *Psiz = box->v[vluDerivs->index[2]];
    double *FChi = box->v[vlFu->index[1]];
    double *Chi  = box->v[vlu->index[1]];
    double *Chix = box->v[vluDerivs->index[9]];
    double *Chiy = box->v[vluDerivs->index[10]];
    double *Chiz = box->v[vluDerivs->index[11]];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;

    /* BCs */
    if(Getv("Poisson_grid", "SphericalDF"))
    {
      forplane1(i,j,k, n1,n2,n3, 0)
      {
        FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*nonlin;
        FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 2.0*nonlin;
      }

      forplane1(i,j,k, n1,n2,n3, n1-1)
      {
        FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*nonlin;
        FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 1.0*nonlin;
      }
    }
    else if (Getv("Poisson_grid", "AnsorgNS"))
    {
      double *P;
      double *dP[4];
      double *C;
      double *dC[4];
      double *BM;

      /* special rho=0 case??? */
      if(b==0 || b==1 || b==2 || b==3)
      {
        int pl;
        char str[1000];
        snprintf(str, 999, "box%d_basis2", b);
        if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* loop over rho=0 boundary */
          for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
            forplane2(i,j,k, n1,n2,n3, pl)
            {
              if(k>0) /* phi>0: impose u_phi_phi=0 */
              {
                FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                FChi[Index(i,j,k)] = Chi_phi_phi[Index(i,j,k)];
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              {
                double Psi_rho = Psiy[Index(i,j,k)];
                double Chi_rho = Chiy[Index(i,j,k)];
                double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                double Chi_rho_phi_phi = Chi_y_phi_phi[Index(i,j,k)];
                FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                FChi[Index(i,j,k)] = Chi_rho + Chi_rho_phi_phi;
              }
            }
        }
        /* same as before, but also interpolate to rho=0 */
        else if(Getv("Poisson_regularization", "regularity_on_axis"))
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];
          double *line = (double *) calloc(n2, sizeof(double));
          double *BM[2];
          BM[0] = (double *) calloc(n2, sizeof(double));
          BM[1] = (double *) calloc(n2, sizeof(double));

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* obtain BM vectors for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[0], 0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[1], 1,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);

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

                /* find value Chi_phi_phi at B=0 or 1 */
                get_memline(Chi_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                FChi[Index(i,j,k)] = U0;
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              { /* Psi_rho = Psiy  
                   Chi_rho = Chiy
                   Psi_rho_phi_phi = Psi_y_phi_phi
                   Chi_rho_phi_phi = Chi_y_phi_phi */
                /* find value Psi_rho at B=0 or 1 */
                get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Psi_rho_phi_phi at B=0 or 1 */
                get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FPsi[Index(i,j,k)] = U0 + V0;

                /* find value Chi_rho at B=0 or 1 */
                get_memline(Chiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Chi_rho_phi_phi at B=0 or 1 */
                get_memline(Chi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FChi[Index(i,j,k)] = U0 + V0;
              }
            }
          }
          free(BM[0]);
          free(BM[1]);
          free(line);
        }
        /* same as before, but do it only in box0/3 at A,B=1,0 and A,B=1,1 */
        else if(Getv("Poisson_regularization", 
                     "regularity_on_axis_at_center") && (b==0 || b==3) )
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];
          double *line = (double *) calloc(n2, sizeof(double));
          double *BM[2];
          BM[0] = (double *) calloc(n2, sizeof(double));
          BM[1] = (double *) calloc(n2, sizeof(double));

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* obtain BM vectors for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[0], 0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[1], 1,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);

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

                /* find value Chi_phi_phi at B=0 or 1 */
                get_memline(Chi_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                FChi[Index(i,j,k)] = U0;
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              { /* Psi_rho = Psiy  
                   Chi_rho = Chiy
                   Psi_rho_phi_phi = Psi_y_phi_phi
                   Chi_rho_phi_phi = Chi_y_phi_phi */
                /* find value Psi_rho at B=0 or 1 */
                get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Psi_rho_phi_phi at B=0 or 1 */
                get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FPsi[Index(i,j,k)] = U0 + V0;

                /* find value Chi_rho at B=0 or 1 */
                get_memline(Chiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Chi_rho_phi_phi at B=0 or 1 */
                get_memline(Chi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FChi[Index(i,j,k)] = U0 + V0;
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
        P = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C = grid->box[1]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0.0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double *line = (double *) calloc(n1, sizeof(double));
          double U0;

          /* find values U0 in box0 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];

          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
          free(line);
        }
        free(BM);
      }
      else if(b==3)  /* in box3 */
      {
        BM = (double *) calloc(n1, sizeof(double));

        /* values at A=0 are equal in box3 and box2 */
        P = grid->box[2]->v[vlu->index[0]]; /* values in box2 */
        C = grid->box[2]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0.0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double *line = (double *) calloc(n1, sizeof(double));
          double U0;

          /* find values U0 in box3 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
//printf("U0=%f P[Index(i,j,k)]=%f\n", U0, P[Index(i,j,k)]);
          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
          free(line);
        }
        free(BM);
      }
      else if(b==1)
      {
        BM = (double *) calloc(max3(grid->box[0]->n1, box->n1, box->n2), 
                               sizeof(double));

        /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
        dP[1] = grid->box[2]->v[vluDerivs->index[0]];
        dC[1] = grid->box[2]->v[vluDerivs->index[9]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          FChi[Index(i,j,k)] = Chix[Index(i,j,k)] - dC[1][Index(i,j,k)];
        }

        /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[0]->v[vluDerivs->index[0]];
        dP[2] = grid->box[0]->v[vluDerivs->index[1]];
        dP[3] = grid->box[0]->v[vluDerivs->index[2]];
        dC[1] = grid->box[0]->v[vluDerivs->index[9]];
        dC[2] = grid->box[0]->v[vluDerivs->index[10]];
        dC[3] = grid->box[0]->v[vluDerivs->index[11]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, grid->box[0]->n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];
          double DP[4], DC[4];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box0 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[0]->n1;
            double *line = (double *) calloc(n1, sizeof(double));

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DC[m] += BM[l]*line[l];
            free(line);
          }
          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - DC[1]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - DC[2]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - DC[3]);
        }

        /* Psi=Chi=0 at infinity */
        if(Getv("box1_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
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
        P  = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C  = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[3]->v[vluDerivs->index[0]];
        dP[2] = grid->box[3]->v[vluDerivs->index[1]];
        dP[3] = grid->box[3]->v[vluDerivs->index[2]];
        dC[1] = grid->box[3]->v[vluDerivs->index[9]];
        dC[2] = grid->box[3]->v[vluDerivs->index[10]];
        dC[3] = grid->box[3]->v[vluDerivs->index[11]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, grid->box[3]->n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];
          double DP[4], DC[4];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box3 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[3]->n1;
            double *line = (double *) calloc(n1, sizeof(double));

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DC[m] += BM[l]*line[l];
            free(line);
          }
          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - DC[1]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - DC[2]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - DC[3]);
        }

        /* Psi=Chi=0 at infinity */
        if(Getv("box2_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
          }
          free(line);
        }

        free(BM);
      }
      else errorexiti("b=%d should be impossible!", b);
    } /* end: else if (Getv("Poisson_grid", "AnsorgNS")) */
    else if (Getv("Poisson_grid", "4ABphi_2xyz"))
    {
      double *P;
      double *dP[4];
      double *C;
      double *dC[4];
      double *X, *Y, *Z;
      double *Pcoeffs, *Ccoeffs;
      double Pinterp, Cinterp;
      double x,y,z;

      /* special rho=0 case??? */
      if(b==0 || b==1 || b==2 || b==3)
      {
        int pl;
        char str[1000];
        snprintf(str, 999, "box%d_basis2", b);
        if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* loop over rho=0 boundary */
          for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
            forplane2(i,j,k, n1,n2,n3, pl)
            {
              if(k>0) /* phi>0: impose u_phi_phi=0 */
              {
                FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                FChi[Index(i,j,k)] = Chi_phi_phi[Index(i,j,k)];
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              {
                double Psi_rho = Psiy[Index(i,j,k)];
                double Chi_rho = Chiy[Index(i,j,k)];
                double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                double Chi_rho_phi_phi = Chi_y_phi_phi[Index(i,j,k)];
                FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                FChi[Index(i,j,k)] = Chi_rho + Chi_rho_phi_phi;
              }
            }
        }
        /* same as before, but also interpolate to rho=0 */
        else if(Getv("Poisson_regularization", "regularity_on_axis"))
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];
          double *line = (double *) calloc(n2, sizeof(double));
          double *BM[2];
          BM[0] = (double *) calloc(n2, sizeof(double));
          BM[1] = (double *) calloc(n2, sizeof(double));

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* obtain BM vectors for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[0], 0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[1], 1,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);

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

                /* find value Chi_phi_phi at B=0 or 1 */
                get_memline(Chi_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                FChi[Index(i,j,k)] = U0;
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              { /* Psi_rho = Psiy  
                   Chi_rho = Chiy
                   Psi_rho_phi_phi = Psi_y_phi_phi
                   Chi_rho_phi_phi = Chi_y_phi_phi */
                /* find value Psi_rho at B=0 or 1 */
                get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Psi_rho_phi_phi at B=0 or 1 */
                get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FPsi[Index(i,j,k)] = U0 + V0;

                /* find value Chi_rho at B=0 or 1 */
                get_memline(Chiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Chi_rho_phi_phi at B=0 or 1 */
                get_memline(Chi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                FChi[Index(i,j,k)] = U0 + V0;
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
        P = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* values at A=Amin are interpolated from box5 */
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[5]->v[vlu->index[0]]; /* values in box5 */
        C = grid->box[5]->v[vlu->index[1]];
        Pcoeffs = grid->box[5]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[5]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[5], P, Pcoeffs);
        spec_Coeffs(grid->box[5], C, Ccoeffs);
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
        {
          int ind=Index(i,j,k);
          x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          Pinterp = spec_interpolate(grid->box[5], Pcoeffs, x,y,z);
          Cinterp = spec_interpolate(grid->box[5], Ccoeffs, x,y,z);
          FPsi[ind] = Psi[ind] - Pinterp;
          FChi[ind] = Chi[ind] - Cinterp;
        }
      }
      else if(b==3)  /* in box3 */
      {
        /* values at A=0 are equal in box3 and box2 */
        P = grid->box[2]->v[vlu->index[0]]; /* values in box2 */
        C = grid->box[2]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* values at A=Amin are interpolated from box4 */
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[4]->v[vlu->index[0]]; /* values in box4 */
        C = grid->box[4]->v[vlu->index[1]];
        Pcoeffs = grid->box[4]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[4]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[4], P, Pcoeffs);
        spec_Coeffs(grid->box[4], C, Ccoeffs);
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
        {
          int ind=Index(i,j,k);
          x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          Pinterp = spec_interpolate(grid->box[4], Pcoeffs, x,y,z);
          Cinterp = spec_interpolate(grid->box[4], Ccoeffs, x,y,z);
          FPsi[ind] = Psi[ind] - Pinterp;
          FChi[ind] = Chi[ind] - Cinterp;
        }
      }
      else if(b==1)  /* in box1 */
      {
        /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
        dP[1] = grid->box[2]->v[vluDerivs->index[0]];
        dC[1] = grid->box[2]->v[vluDerivs->index[9]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          FChi[Index(i,j,k)] = Chix[Index(i,j,k)] - dC[1][Index(i,j,k)];
        }

        /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[0]->v[vluDerivs->index[0]];
        dP[2] = grid->box[0]->v[vluDerivs->index[1]];
        dP[3] = grid->box[0]->v[vluDerivs->index[2]];
        dC[1] = grid->box[0]->v[vluDerivs->index[9]];
        dC[2] = grid->box[0]->v[vluDerivs->index[10]];
        dC[3] = grid->box[0]->v[vluDerivs->index[11]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];

          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - dC[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - dC[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - dC[3][Index(i,j,k)]);
        }
        
        /* Psi=Chi=0 at infinity */
        if(Getv("box1_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *BM = (double *) calloc(n2, sizeof(double));
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
          }
          free(line);
          free(BM);
        }
      }
      else if(b==2)  /* in box2 */
      {
        /* values at A=1 are equal in box1 and box2 */
        P  = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C  = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }

        /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[3]->v[vluDerivs->index[0]];
        dP[2] = grid->box[3]->v[vluDerivs->index[1]];
        dP[3] = grid->box[3]->v[vluDerivs->index[2]];
        dC[1] = grid->box[3]->v[vluDerivs->index[9]];
        dC[2] = grid->box[3]->v[vluDerivs->index[10]];
        dC[3] = grid->box[3]->v[vluDerivs->index[11]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];

          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - dC[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - dC[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - dC[3][Index(i,j,k)]);
        }

        /* Psi=Chi=0 at infinity */
        if(Getv("box2_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *BM = (double *) calloc(n2, sizeof(double));
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
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
        P = grid->box[0]->v[vlu->index[0]]; /* values in box0 */
        C = grid->box[0]->v[vlu->index[1]];
        Pcoeffs = grid->box[0]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[0]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[0], P, Pcoeffs);
        spec_Coeffs(grid->box[0], C, Ccoeffs);
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
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);

            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
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
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
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
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
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
        P = grid->box[3]->v[vlu->index[0]]; /* values in box3 */
        C = grid->box[3]->v[vlu->index[1]];
        Pcoeffs = grid->box[3]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[3]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[3], P, Pcoeffs);
        spec_Coeffs(grid->box[3], C, Ccoeffs);
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
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
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
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
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
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
      }
      else errorexiti("b=%d should be impossible!", b);

    } /* end: else if (Getv("Poisson_grid", "4ABphi_2xyz")) */

  }
}
