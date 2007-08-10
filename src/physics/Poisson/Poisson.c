/* Poisson.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "Poisson.h"

#define Power pow


void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);


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
      else if (Getv("Poisson_grid", "AnsorgNS"))
      {
        Psi[i] = 0.0;
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
        if(b==1)  rh1[i] = 4.0*PI;
        if(b==1)  rh2[i] = 8.0*PI;
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
  double normresnonlin;

  tVarList *vlu, *vlFu, *vluDerivs;
  tVarList *vldu, *vlr, *vlduDerivs;
  tVarList *vlrhs;

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
    
  /* enable vlu, vluDerivs, vlrhs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 
  enablevarlist(vlrhs);

  /* now duplicate vlu for result of F(u) */
  vlFu = AddDuplicateEnable(vlu,  "res");

  /* now duplicate vlu and vluDerivs for linarized Eqs. */
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlr        = AddDuplicateEnable(vlFu, "_l");
  vlduDerivs = vluDerivs; /* maybe: vlduDerivs=AddDuplicateEnable(vluDerivs, "_l"); */

  /* call Newton solver */
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, vlrhs,
         itmax, tol, &normresnonlin, 1,
         bicgstab, Precon_I, vldu, vlr, vlduDerivs, vlrhs,
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
        PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
        ChiErr[i] = Chi[i]-2.0/sqrt(x*x + y*y + z*z);
      }
      else if (Getv("Poisson_grid", "AnsorgNS"))
      {
        if(b==1||b==2)
        {
          PsiErr[i] = Psi[i]-1.0/sqrt((x-1)*(x-1) + y*y + z*z);
          ChiErr[i] = Chi[i]-2.0/sqrt((x+1)*(x+1) + y*y + z*z);
        }
      }
    }
  }
  return 0;
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
      FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i] - rh1[i];
      FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i] - rh2[i];
    }

    /* BCs */
    set_BCs(vlFu, vlu, vluDerivs, 1);
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
    set_BCs(vlJdu, vldu, vlduDerivs, 0);
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


void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
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

      if(b==0)  /* in box0 */
      {
        BM = (double *) calloc(n1, sizeof(double));

        /* values at A=0 are equal in box0 and box1 */
        P = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C = grid->box[1]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double line[n1]; // fixme!!!!!!!
          double U0;

          /* find values U0 in box0 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<=n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];

          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<=n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
        }
        free(BM);
      }
      else if(b==3)  /* in box3 */
      {
        BM = (double *) calloc(n1, sizeof(double));

        /* values at A=0 are equal in box3 and box2 */
        P = grid->box[2]->v[vlu->index[0]]; /* values in box2 */
        C = grid->box[2]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double line[n1]; // fixme!!!!!!!
          double U0;

          /* find values U0 in box3 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<=n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];

          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<=n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
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
          double DP[3], DC[3];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box0 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[0]->n1;
            double line[n1]; // fixme!!!!!!!

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<=n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<=n1; l++)  DC[m] += BM[l]*line[l];
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
        // fixme: B=0 is not on grid!!!
        for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
        {
          FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
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
          double DP[3], DC[3];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box3 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[3]->n1;
            double line[n1]; // fixme!!!!!!!

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<=n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<=n1; l++)  DC[m] += BM[l]*line[l];
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
        // fixme: B=0 is not on grid!!!
        for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
        {
          FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
        }

        free(BM);
      }
      else errorexiti("b=%d should be impossible!", b);
    } /* end: else if (Getv("Poisson_grid", "AnsorgNS")) */

  }
}
