/* ScalarOnKerr.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "ScalarOnKerr.h"



/* initialize ScalarOnKerr */
int ScalarOnKerr_startup(tGrid *grid)
{
  tVarList *ScalarOnKerrvars;
  int b;
  
  printf("Initializing ScalarOnKerr:\n");

  /* set boundary information for ScalarOnKerr evolution: 
     farlimit, falloff, propagation speed   */
  VarNameSetBoundaryInfo("ScalarOnKerr_psi",     0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsix",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsiy",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsiz",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsixx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsixy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsixz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsiyy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsiyz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_ddpsizz", 0, 1, 1.0);

  VarNameSetBoundaryInfo("ScalarOnKerr_Pi",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dPix", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dPiy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dPiz", 0, 1, 1.0);

  /* create a variable list for ScalarOnKerr evolutions  */
  ScalarOnKerrvars = vlalloc(grid);
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_psi"));
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_Pi"));
  if (0) prvarlist(ScalarOnKerrvars);
  enablevarlist(ScalarOnKerrvars);

  /* register evolved variables */
  evolve_vlregister(ScalarOnKerrvars);
  
  /* register evolution routine */
  evolve_rhsregister(ScalarOnKerr_evolve);

  /* register BC routine */
//  evolve_algebraicConditionsregister(set_boundary_ofPi);

  /* set initial data in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i,j,k;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];

    double *psi = box->v[Ind("ScalarOnKerr_psi")];
    double *Pi  = box->v[Ind("ScalarOnKerr_Pi")];

    forallpoints(box,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];

      psi[i] =  0.0;
      Pi[i]  =  0.0;
    }
  }

  //set_boundary_symmetry(level, ScalarOnKerrvars);

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(ScalarOnKerrvars);
  
  /* enable all derivative vars */
  enablevar(grid, Ind("ScalarOnKerr_dpsix"));
  enablevar(grid, Ind("ScalarOnKerr_ddpsixx"));
  enablevar(grid, Ind("ScalarOnKerr_dPix"));
  
  /* enable all metric vars */
  enablevar(grid, Ind("ScalarOnKerr_gtt"));
  enablevar(grid, Ind("ScalarOnKerr_guptt"));
  enablevar(grid, Ind("ScalarOnKerr_Gammattt"));

  enablevar(grid, Ind("ScalarOnKerr3d_gxx"));
  enablevar(grid, Ind("ScalarOnKerr3d_alpha"));
  enablevar(grid, Ind("ScalarOnKerr3d_betax"));
  enablevar(grid, Ind("ScalarOnKerr3d_Kxx"));
  enablevar(grid, Ind("ScalarOnKerr3d_TrK"));
  enablevar(grid, Ind("ScalarOnKerr3d_gupxx"));
  enablevar(grid, Ind("ScalarOnKerr3d_Gammaxxx"));
  enablevar(grid, Ind("ScalarOnKerr3d_dalphax"));

  /* enable temp vars */
  if(!Getv("physics", "ADMvars"))
  {
    enablevar(grid, Ind("temp1"));
    enablevar(grid, Ind("temp2"));
    enablevar(grid, Ind("temp3"));
  }


  /* set Kerr metric and Christoffels */
  Kerr(grid, Ind("x"), Ind("ScalarOnKerr_gtt"), Ind("ScalarOnKerr_guptt"),
       Ind("ScalarOnKerr_Gammattt"));

  /* set Kerr in 3+1 */
  Kerr3d(grid, Ind("x"), Ind("ScalarOnKerr3d_alpha"), Ind("ScalarOnKerr3d_betax"),
         Ind("ScalarOnKerr3d_gxx"), Ind("ScalarOnKerr3d_Kxx"),
         Ind("ScalarOnKerr3d_TrK"), Ind("ScalarOnKerr3d_gupxx"),
         Ind("ScalarOnKerr3d_Gammaxxx"), Ind("ScalarOnKerr3d_dalphax"));

  return 0;
}




/* evolve and set boundary points */
void ScalarOnKerr_evolve(tVarList *unew, tVarList *upre, double dt, 
                             tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    double *cpsi = vlldataptr(ucur, box, 0);
    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *pPi = vlldataptr(upre, box, 1);
    double *nPi = vlldataptr(unew, box, 1);
    int i;
    int ipsi 	= (ucur)->index[0];
    int iPi = (ucur)->index[1];
    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
    double *psixx = box->v[Ind("ScalarOnKerr_ddpsixx")];
    double *psixy = box->v[Ind("ScalarOnKerr_ddpsixx")+1];
    double *psixz = box->v[Ind("ScalarOnKerr_ddpsixx")+2];
    double *psiyy = box->v[Ind("ScalarOnKerr_ddpsixx")+3];
    double *psiyz = box->v[Ind("ScalarOnKerr_ddpsixx")+4];
    double *psizz = box->v[Ind("ScalarOnKerr_ddpsixx")+5];
    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    int i_gup = Ind("ScalarOnKerr_guptt");
    double *gtt = box->v[i_gup];
    double *gtx = box->v[i_gup+1];
    double *gty = box->v[i_gup+2];
    double *gtz = box->v[i_gup+3];
    double *gxx = box->v[i_gup+4];
    double *gxy = box->v[i_gup+5];
    double *gxz = box->v[i_gup+6];
    double *gyy = box->v[i_gup+7];
    double *gyz = box->v[i_gup+8];
    double *gzz = box->v[i_gup+9];
    int i_Gam = Ind("ScalarOnKerr_Gammattt");
    double *Gamttt = box->v[i_Gam];
    double *Gamttx = box->v[i_Gam+1];
    double *Gamtty = box->v[i_Gam+2];
    double *Gamttz = box->v[i_Gam+3];
    double *Gamtxx = box->v[i_Gam+4];
    double *Gamtxy = box->v[i_Gam+5];
    double *Gamtxz = box->v[i_Gam+6];
    double *Gamtyy = box->v[i_Gam+7];
    double *Gamtyz = box->v[i_Gam+8];
    double *Gamtzz = box->v[i_Gam+9];
    double *Gamxtt = box->v[i_Gam+10];
    double *Gamxtx = box->v[i_Gam+11];
    double *Gamxty = box->v[i_Gam+12];
    double *Gamxtz = box->v[i_Gam+13];
    double *Gamxxx = box->v[i_Gam+14];
    double *Gamxxy = box->v[i_Gam+15];
    double *Gamxxz = box->v[i_Gam+16];
    double *Gamxyy = box->v[i_Gam+17];
    double *Gamxyz = box->v[i_Gam+18];
    double *Gamxzz = box->v[i_Gam+19];
    double *Gamytt = box->v[i_Gam+20];
    double *Gamytx = box->v[i_Gam+21];
    double *Gamyty = box->v[i_Gam+22];
    double *Gamytz = box->v[i_Gam+23];
    double *Gamyxx = box->v[i_Gam+24];
    double *Gamyxy = box->v[i_Gam+25];
    double *Gamyxz = box->v[i_Gam+26];
    double *Gamyyy = box->v[i_Gam+27];
    double *Gamyyz = box->v[i_Gam+28];
    double *Gamyzz = box->v[i_Gam+29];
    double *Gamztt = box->v[i_Gam+30];
    double *Gamztx = box->v[i_Gam+31];
    double *Gamzty = box->v[i_Gam+32];
    double *Gamztz = box->v[i_Gam+33];
    double *Gamzxx = box->v[i_Gam+34];
    double *Gamzxy = box->v[i_Gam+35];
    double *Gamzxz = box->v[i_Gam+36];
    double *Gamzyy = box->v[i_Gam+37];
    double *Gamzyz = box->v[i_Gam+38];
    double *Gamzzz = box->v[i_Gam+39];
            
    /* compute the spatial derivs */
    allDerivsOf_S(box, ipsi,
                  Ind("ScalarOnKerr_dpsix"), Ind("ScalarOnKerr_ddpsixx"));
    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forallpoints(box, i)
    {
      double rPi, rpsi;
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double t = grid->time;
      double x0, y0;
      double g_ddpsi, gGt,gGx,gGy,gGz, gG_dpsi;
      /* g is upper metric */
      /* get all terms with less than 2 time derivs in g^ab d_a d_b psi */
      g_ddpsi = -2.0*(gtx[i]*Pix[i] +gty[i]*Piy[i] +gtz[i]*Piz[i]) + 
                gxx[i]*psixx[i] + gyy[i]*psiyy[i] + gzz[i]*psizz[i] +
                2.0*(gxy[i]*psixy[i] + gxz[i]*psixz[i] + gyz[i]*psiyz[i]);

      /* get gG[a] = g^bc Gamma^a_bc */
      gGt = gtt[i]*Gamttt[i] + 
            gxx[i]*Gamtxx[i] + gyy[i]*Gamtyy[i] + gzz[i]*Gamtzz[i] +
            2.0*(gtx[i]*Gamttx[i] + gty[i]*Gamtty[i] + gtz[i]*Gamttz[i] +
                 gxy[i]*Gamtxy[i] + gxz[i]*Gamtxz[i] + gyz[i]*Gamtyz[i]);
      gGx = gtt[i]*Gamxtt[i] + 
            gxx[i]*Gamxxx[i] + gyy[i]*Gamxyy[i] + gzz[i]*Gamxzz[i] +
            2.0*(gtx[i]*Gamxtx[i] + gty[i]*Gamxty[i] + gtz[i]*Gamxtz[i] +
                 gxy[i]*Gamxxy[i] + gxz[i]*Gamxxz[i] + gyz[i]*Gamxyz[i]);
      gGy = gtt[i]*Gamytt[i] + 
            gxx[i]*Gamyxx[i] + gyy[i]*Gamyyy[i] + gzz[i]*Gamyzz[i] +
            2.0*(gtx[i]*Gamytx[i] + gty[i]*Gamyty[i] + gtz[i]*Gamytz[i] +
                 gxy[i]*Gamyxy[i] + gxz[i]*Gamyxz[i] + gyz[i]*Gamyyz[i]);
      gGz = gtt[i]*Gamztt[i] + 
            gxx[i]*Gamzxx[i] + gyy[i]*Gamzyy[i] + gzz[i]*Gamzzz[i] +
            2.0*(gtx[i]*Gamztx[i] + gty[i]*Gamzty[i] + gtz[i]*Gamztz[i] +
                 gxy[i]*Gamzxy[i] + gxz[i]*Gamzxz[i] + gyz[i]*Gamzyz[i]);
      gG_dpsi = -gGt*cPi[i] + gGx*psix[i] + gGy*psiy[i] + gGz*psiz[i];

      /* source position */
      x0 = 10*cos(0.02*t);
      y0 = 10*sin(0.02*t);

      /* set RHS of psi and Pi */
      rPi  = (g_ddpsi - gG_dpsi)/gtt[i] + 
             -(1-Attenuation01( ((x-x0)*(x-x0)+(y-y0)*(y-y0)+z*z)/36,2,0.5));
//             -exp(-(x-x0)*(x-x0))*exp(-(y-y0)*(y-y0))*exp(-z*z); // source
      rpsi = -cPi[i];

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        nPi[i]  = pPi[i]  + dt * rPi;
        npsi[i] = ppsi[i] + dt * rpsi;
      }
      else
      {
        nPi[i]  = rPi;
        npsi[i] = rpsi;
      }
    }
  } /* end forallboxes */
      
  /* set BCs */
  set_boundary(unew, upre, dt, ucur);
//  set_psi_Pi_boundary(unew, upre, dt, ucur);

  /* filter near poles */
  coordinateDependentFilter(unew);

  /* filter high freq. angular modes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, jf,kf;
    double *nPi = vlldataptr(unew, box, 1);
    double *temp1 = box->v[Ind("temp1")];

    /* filter nPi */
    spec_Coeffs(box, nPi, temp1);
    kf=n3/3; kf*=2;
    jf=n2/3; jf*=2;
    forallijk(i,j,k)
      if(k>kf || j>jf) temp1[Index(i,j,k)]=0.0;
    spec_Eval(box, nPi, temp1);
  } /* end forallboxes */

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}


/* evolve and set boundary points */
void ScalarOnKerr_evolve_new(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    double *cpsi = vlldataptr(ucur, box, 0);
    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *pPi = vlldataptr(upre, box, 1);
    double *nPi = vlldataptr(unew, box, 1);
    int i;
    int ipsi 	= (ucur)->index[0];
    int iPi = (ucur)->index[1];
    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
    double *psixx = box->v[Ind("ScalarOnKerr_ddpsixx")];
    double *psixy = box->v[Ind("ScalarOnKerr_ddpsixx")+1];
    double *psixz = box->v[Ind("ScalarOnKerr_ddpsixx")+2];
    double *psiyy = box->v[Ind("ScalarOnKerr_ddpsixx")+3];
    double *psiyz = box->v[Ind("ScalarOnKerr_ddpsixx")+4];
    double *psizz = box->v[Ind("ScalarOnKerr_ddpsixx")+5];
    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    int i_alpha = Ind("ScalarOnKerr3d_alpha");
    double *alpha = box->v[i_alpha];
    int i_beta = Ind("ScalarOnKerr3d_betax");
    double *betax = box->v[i_beta];
    double *betay = box->v[i_beta+1];
    double *betaz = box->v[i_beta+2];
    int i_gup = Ind("ScalarOnKerr3d_gupxx");
    double *gxx = box->v[i_gup];
    double *gxy = box->v[i_gup+1];
    double *gxz = box->v[i_gup+2];
    double *gyy = box->v[i_gup+3];
    double *gyz = box->v[i_gup+4];
    double *gzz = box->v[i_gup+5];
    int i_TrK = Ind("ScalarOnKerr3d_TrK");
    double *TrK = box->v[i_TrK];
    int i_Gam = Ind("ScalarOnKerr3d_Gammaxxx");
    double *Gamxxx = box->v[i_Gam];
    double *Gamxxy = box->v[i_Gam+1];
    double *Gamxxz = box->v[i_Gam+2];
    double *Gamxyy = box->v[i_Gam+3];
    double *Gamxyz = box->v[i_Gam+4];
    double *Gamxzz = box->v[i_Gam+5];
    double *Gamyxx = box->v[i_Gam+6];
    double *Gamyxy = box->v[i_Gam+7];
    double *Gamyxz = box->v[i_Gam+8];
    double *Gamyyy = box->v[i_Gam+9];
    double *Gamyyz = box->v[i_Gam+10];
    double *Gamyzz = box->v[i_Gam+11];
    double *Gamzxx = box->v[i_Gam+12];
    double *Gamzxy = box->v[i_Gam+13];
    double *Gamzxz = box->v[i_Gam+14];
    double *Gamzyy = box->v[i_Gam+15];
    double *Gamzyz = box->v[i_Gam+16];
    double *Gamzzz = box->v[i_Gam+17];
    int i_dalpha = Ind("ScalarOnKerr3d_dalphax");
    double *dalphax = box->v[i_dalpha];
    double *dalphay = box->v[i_dalpha+1];
    double *dalphaz = box->v[i_dalpha+2];
            
    /* compute the spatial derivs */
    allDerivsOf_S(box, ipsi,
                  Ind("ScalarOnKerr_dpsix"), Ind("ScalarOnKerr_ddpsixx"));
    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forallpoints(box, i)
    {
      double rPi, rpsi;
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double t = grid->time;
      double x0, y0;
      double beta_dPi,ag_ddpsi, gGx,gGy,gGz,agG_dpsi, g_dpsi_da, aKPi, beta_dpsi;
      /* g is upper metric */
      /* terms on RHS of Pi eqn */ 
      beta_dPi = betax[i]*Pix[i] + betay[i]*Piy[i] + betaz[i]*Piz[i];
      ag_ddpsi = alpha[i]*
                  (gxx[i]*psixx[i] + gyy[i]*psiyy[i] + gzz[i]*psizz[i] +
                   2.0*(gxy[i]*psixy[i] + gxz[i]*psixz[i] + gyz[i]*psiyz[i]));
      
      /* get alpha gG[i] = alpha g^ik Gamma^i_jk */
      gGx = gxx[i]*Gamxxx[i] + gyy[i]*Gamxyy[i] + gzz[i]*Gamxzz[i] +
             2.0*(gxy[i]*Gamxxy[i] + gxz[i]*Gamxxz[i] + gyz[i]*Gamxyz[i]);
      gGy = gxx[i]*Gamyxx[i] + gyy[i]*Gamyyy[i] + gzz[i]*Gamyzz[i] +
             2.0*(gxy[i]*Gamyxy[i] + gxz[i]*Gamyxz[i] + gyz[i]*Gamyyz[i]);
      gGz = gxx[i]*Gamzxx[i] + gyy[i]*Gamzyy[i] + gzz[i]*Gamzzz[i] +
             2.0*(gxy[i]*Gamzxy[i] + gxz[i]*Gamzxz[i] + gyz[i]*Gamzyz[i]);
      agG_dpsi = alpha[i]*(gGx*psix[i] + gGy*psiy[i] + gGz*psiz[i]);

      /* get g_dpsi_da = g^ik dpsi_i dalpha_k */
      g_dpsi_da = gxx[i]*psix[i]*dalphax[i] + gyy[i]*psiy[i]*dalphay[i] +
                  gzz[i]*psiz[i]*dalphaz[i] +
                  gxy[i]*(psix[i]*dalphay[i] + psiy[i]*dalphax[i]) + 
                  gxz[i]*(psix[i]*dalphaz[i] + psiz[i]*dalphax[i]) +
                  gyz[i]*(psiy[i]*dalphaz[i] + psiz[i]*dalphay[i]);

      aKPi = alpha[i]*(TrK[i]*cPi[i]);

      /* term on RHS of psi eqn */
      beta_dpsi = betax[i]*psix[i] + betay[i]*psiy[i] + betaz[i]*psiz[i];

      /* source position */
      x0 = 10*cos(0.02*t);
      y0 = 10*sin(0.02*t);

      /* set RHS of psi and Pi */
      rPi  = beta_dPi - ag_ddpsi + agG_dpsi - g_dpsi_da + aKPi  +
              (-1.0/alpha[i])*
              (1-Attenuation01( ((x-x0)*(x-x0)+(y-y0)*(y-y0)+z*z)/36,2,0.5));
//              (-1/alpha[i])*exp(-(x-x0)*(x-x0))*exp(-(y-y0)*(y-y0))*exp(-z*z); // source
      rpsi = beta_dpsi - alpha[i]*cPi[i];

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        nPi[i]  = pPi[i]  + dt * rPi;
        npsi[i] = ppsi[i] + dt * rpsi;
      }
      else
      {
        nPi[i]  = rPi;
        npsi[i] = rpsi;
      }
    }
  }

  /* set BCs */
  set_boundary(unew, upre, dt, ucur);

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}

/* set BC for Pi */
void set_boundary_ofPi(tVarList *unew, tVarList *upre)
{
  int bi, pi, ijk;
  int ipsi = unew->index[0];
  int iPi  = unew->index[1];
  tGrid *grid = (unew)->grid;

  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *xp = box->v[Ind("x")];
    double *yp = box->v[Ind("y")];
    double *zp = box->v[Ind("z")];
    double *psinew = box->v[ipsi];
    double *Pinew = box->v[iPi];
    double *dpsi_dx = box->v[Ind("temp1")];
    double *dpsi_dy = box->v[Ind("temp2")];
    double *dpsi_dz = box->v[Ind("temp3")];
    double x,y,z, r, nx,ny,nz;

    cart_partials(box, psinew, dpsi_dx, dpsi_dy, dpsi_dz); 

    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      ijk = Index(i,j,k);
      x = xp[ijk];
      y = yp[ijk];
      z = zp[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r;
      ny = y/r;
      nz = z/r;

      Pinew[ijk] = (nx*dpsi_dx[ijk] + ny*dpsi_dy[ijk] + nz*dpsi_dz[ijk])
                    - psinew[ijk]/r;
    }
  }
}

/* set BC for old version */
void set_psi_Pi_boundary(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *cpsi = vlldataptr(ucur, box, 0);
    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *pPi = vlldataptr(upre, box, 1);
    double *nPi = vlldataptr(unew, box, 1);
    int ipsi = (ucur)->index[0];
    int iPi  = (ucur)->index[1];
    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    int i_gup = Ind("ScalarOnKerr_guptt");
    double *gtt = box->v[i_gup];
    double *gtx = box->v[i_gup+1];
    double *gty = box->v[i_gup+2];
    double *gtz = box->v[i_gup+3];
    double *gxx = box->v[i_gup+4];
    double *gxy = box->v[i_gup+5];
    double *gxz = box->v[i_gup+6];
    double *gyy = box->v[i_gup+7];
    double *gyz = box->v[i_gup+8];
    double *gzz = box->v[i_gup+9];
    int i_g3up = Ind("ScalarOnKerr3d_gupxx");
    double *g3xx = box->v[i_g3up];
    double *g3xy = box->v[i_g3up+1];
    double *g3xz = box->v[i_g3up+2];
    double *g3yy = box->v[i_g3up+3];
    double *g3yz = box->v[i_g3up+4];
    double *g3zz = box->v[i_g3up+5];
    int i_alpha = Ind("ScalarOnKerr3d_alpha");
    double *alpha = box->v[i_alpha];
    int i_beta = Ind("ScalarOnKerr3d_betax");
    double *betax = box->v[i_beta];
    double *betay = box->v[i_beta+1];
    double *betaz = box->v[i_beta+2];
            
    /* compute the spatial derivs */
    FirstDerivsOf_S(box, ipsi, Ind("ScalarOnKerr_dpsix"));
    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      double x,y,z;
      double r, nx,ny,nz, Nx,Ny,Nz;
      double rPi, rpsi, vx,vy,vz;

      ijk = Index(i,j,k);
      x = px[ijk];
      y = py[ijk];
      z = pz[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r;
      ny = y/r;
      nz = z/r;
      Nx = g3xx[ijk]*nx + g3xy[ijk]*ny + g3xz[ijk]*nz;
      Ny = g3xy[ijk]*nx + g3yy[ijk]*ny + g3yz[ijk]*nz;
      Nz = g3xz[ijk]*nx + g3yz[ijk]*ny + g3zz[ijk]*nz;
//Nx=nx; Ny=ny; Nz=nz;  // this lasts longer
//  Nx = gxx[ijk]*nx + gxy[ijk]*ny + gxz[ijk]*nz;
//  Ny = gxy[ijk]*nx + gyy[ijk]*ny + gyz[ijk]*nz;
//  Nz = gxz[ijk]*nx + gyz[ijk]*ny + gzz[ijk]*nz;

      vx = betax[ijk]-alpha[ijk]*Nx;
      vy = betay[ijk]-alpha[ijk]*Ny;
      vz = betaz[ijk]-alpha[ijk]*Nz;
//vx = -nx;
//vy = -ny;
//vz = -nz;

      /* set RHS of psi and Pi */
//      rPi  = vx* Pix[ijk] + vy* Piy[ijk] + vz* Piz[ijk] -0* cPi[ijk]/r;
//      rpsi = vx*psix[ijk] + vy*psiy[ijk] + vz*psiz[ijk] -0*cpsi[ijk]/r;
rPi  = vx* Pix[ijk] + vy* Piy[ijk] + vz* Piz[ijk] - 
       (cPi[ijk] + betax[ijk]*psix[ijk] + betay[ijk]*psiy[ijk] +
                   betaz[ijk]*psiz[ijk])/r;
//rpsi = vx*psix[ijk] + vy*psiy[ijk] + vz*psiz[ijk] -
//       0*(cpsi[ijk] + betax[ijk]*psix[ijk] + betay[ijk]*psiy[ijk] +
//                    betaz[ijk]*psiz[ijk])/r;

if(0)
{
printf("(%f,%f,%f), (%f,%f,%f), r=%f\n", nx,ny,nz, x,y,z, r);
printf("(%d,%d,%d), ijk=%d\n", i,j,k, ijk);
printf("%g %g ", psix[ijk], Pix[ijk]);
}

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        nPi[ijk]  = pPi[ijk]  + dt * rPi;
//        npsi[ijk] = ppsi[ijk] + dt * rpsi;
      }
      else
      {
        nPi[ijk]  = rPi;
//        npsi[ijk] = rpsi;
      }
    }
  }
}


/* compute and integrals of rho */
int ScalarOnKerr_analyze(tGrid *grid)
{
  int b;

// if( ! timeforoutput_index(grid, Ind("ScalarOnKerr_rho")) ) return 0;

// do something useful here 

  return 0;
}
