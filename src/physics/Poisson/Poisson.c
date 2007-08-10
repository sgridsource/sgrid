/* Poisson.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "Poisson.h"

#define Power pow


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
  enablevar(grid, Ind("Poisson_Psix"));
  enablevar(grid, Ind("Poisson_Psixx"));
  enablevar(grid, Ind("Poisson_Err_Psi"));
  enablevar(grid, Ind("Poisson_Chi"));
  enablevar(grid, Ind("Poisson_Chix"));
  enablevar(grid, Ind("Poisson_Chixx"));
  enablevar(grid, Ind("Poisson_Err_Chi"));
  
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
    double *Psi = box->v[Ind("Poisson_Psi")];
    double *Chi = box->v[Ind("Poisson_Chi")];

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
      Psi[i] = x*y*z; // 1.0/sqrt(x*x + y*y + z*z);
      Chi[i] = 0.0;
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
  double normresnonlin;

  tVarList *vlu, *vlFu, *vluDerivs;
  tVarList *vldu, *vlr, *vlduDerivs;

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vlFu = vlalloc(grid);
  vluDerivs= vlalloc(grid);
  
  /* add Poisson_Psi and Poisson_Chi to vlu */
  vlpush(vlu, Ind("Poisson_Psi"));
  vlpush(vlu, Ind("Poisson_Chi"));

  /* add derivs of Poisson_Psi and Poisson_Chi to vluDerivs */
  vlpush(vluDerivs, Ind("Poisson_Psix"));
  vlpush(vluDerivs, Ind("Poisson_Psixx"));
  vlpush(vluDerivs, Ind("Poisson_Chix"));
  vlpush(vluDerivs, Ind("Poisson_Chixx"));

  /* enable vlu, vluDerivs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 

  /* now duplicate vlu for result of F(u) */
  vlFu = AddDuplicateEnable(vlu,  "res");

  /* now duplicate vlu and vluDerivs for linarized Eqs. */
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlr        = AddDuplicateEnable(vlFu, "_l");
  vlduDerivs = vluDerivs; /* maybe: vlduDerivs=AddDuplicateEnable(vluDerivs, "_l"); */

  /* call Newton solver */
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, 0,
         itmax, tol, &normresnonlin, 1,
         bicgstab, Precon_I, vldu, vlr, vlduDerivs, vlu,
         linSolver_itmax, linSolver_tolFac);
/*
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, 0,
         itmax, tol, &normresnonlin, 1,
         linSol, Precon_I, vldu, vlr, vlduDerivs, vlu,
         linSolver_itmax, linSolver_tolFac);
*/
  /* free varlists */     
  VLDisableFree(vldu);
  /* VLDisableFree(vlduDerivs); */
  VLDisableFree(vlr);     
  vlfree(vlu);
  vlfree(vluDerivs);
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

      PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
      ChiErr[i] = Chi[i]-2.0/sqrt(x*x + y*y + z*z);
    }
  }
  return 0;
}


/* evaluate Poisson eqn for vlu */
void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2)
{
  tGrid *grid = vlu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[vlFu->index[0]];
    double *Psi  = box->v[vlu->index[0]];
    double *Psixx = box->v[vluDerivs->index[3]];
    double *Psiyy = box->v[vluDerivs->index[6]];
    double *Psizz = box->v[vluDerivs->index[8]];
    double *FChi = box->v[vlFu->index[1]];
    double *Chi  = box->v[vlu->index[1]];
    double *Chixx = box->v[vluDerivs->index[12]];
    double *Chiyy = box->v[vluDerivs->index[15]];
    double *Chizz = box->v[vluDerivs->index[17]];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes"))
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
//    forallpoints(box, i)  FPsi[i] = Psi[i]*Psi[i] - 9.0;
    forallpoints(box, i)
    {
      FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i];
      FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i];
    }

    /* BCs */
    forplane1(i,j,k, n1,n2,n3, 0)
    {
      FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0;
      FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 2.0;
    }

    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5;
      FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 1.0;
    }
  }
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
    double *dPsi  = box->v[vldu->index[0]];
    double *dPsixx = box->v[vlduDerivs->index[3]];
    double *dPsiyy = box->v[vlduDerivs->index[6]];
    double *dPsizz = box->v[vlduDerivs->index[8]];
    double *JdChi = box->v[vlJdu->index[1]];
    double *dChi  = box->v[vldu->index[1]];
    double *dChixx = box->v[vlduDerivs->index[12]];
    double *dChiyy = box->v[vlduDerivs->index[15]];
    double *dChizz = box->v[vlduDerivs->index[17]];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes"))
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

    /* BCs */
    forplane1(i,j,k, n1,n2,n3, 0)
    {
      JdPsi[Index(i,j,k)] = dPsi[Index(i,j,k)];
      JdChi[Index(i,j,k)] = dChi[Index(i,j,k)];
    }

    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      JdPsi[Index(i,j,k)] = dPsi[Index(i,j,k)];
      JdChi[Index(i,j,k)] = dChi[Index(i,j,k)];
    }
  }
}

/* evaluate linearized Laplace */
void J_Poisson_old(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  tGrid *grid = vldu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *Jdu = box->v[vlJdu->index[0]];
    double *du  = box->v[vldu->index[0]];
    double *duxx = box->v[vlduDerivs->index[3]];
    double *duyy = box->v[vlduDerivs->index[6]];
    double *duzz = box->v[vlduDerivs->index[8]];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes"))
      allDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                    vlduDerivs->index[3]);
    else
      FirstAndSecondDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                               vlduDerivs->index[3]);

//    forallpoints(box, i)  Jdu[i] = 2.0*u[i]*du[i];
    forallpoints(box, i)  Jdu[i] = duxx[i] + duyy[i] + duzz[i];

    /* BCs */
    forplane1(i,j,k, n1,n2,n3, 0)
      Jdu[Index(i,j,k)] = du[Index(i,j,k)];

    forplane1(i,j,k, n1,n2,n3, n1-1)
      Jdu[Index(i,j,k)] = du[Index(i,j,k)];
  }
}  

/* Einheitsmatrix als Precon */ 
void Precon_I(tVarList *vlJdu, tVarList *vldu,
              tVarList *vlduDerivs, tVarList *vlu)
{
  tGrid *grid = vldu->grid;
  int i,j,b;
  	
  for(j = 0; j < vldu->n; j++)
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *Jdu = box->v[vlJdu->index[j]];
      double *du  = box->v[vldu->index[j]];

      forallpoints(box, i)  Jdu[i] = du[i];
    }
}
