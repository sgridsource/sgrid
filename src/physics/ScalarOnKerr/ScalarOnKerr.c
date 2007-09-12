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

  VarNameSetBoundaryInfo("ScalarOnKerr_psidot",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsidotx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsidoty", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_dpsidotz", 0, 1, 1.0);

  /* create a variable list for ScalarOnKerr evolutions  */
  ScalarOnKerrvars = vlalloc(grid);
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_psi"));
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_psidot"));
  if (0) prvarlist(ScalarOnKerrvars);
  enablevarlist(ScalarOnKerrvars);

  /* register evolved variables */
  evolve_vlregister(ScalarOnKerrvars);
  
  /* register evolution routine */
  evolve_rhsregister(ScalarOnKerr_evolve);

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
    double *psidot= box->v[Ind("ScalarOnKerr_psidot")];

    forallpoints(box,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];

      psi[i]    =  0.0;
      psidot[i] =  0.0;
    }
  }

  //set_boundary_symmetry(level, ScalarOnKerrvars);

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(ScalarOnKerrvars);
  
  /* enable all derivative vars */
  enablevar(grid, Ind("ScalarOnKerr_dpsix"));
  enablevar(grid, Ind("ScalarOnKerr_ddpsixx"));
  enablevar(grid, Ind("ScalarOnKerr_dpsidotx"));
  
  /* enable all metric vars */
  enablevar(grid, Ind("ScalarOnKerr_gtt"));
  enablevar(grid, Ind("ScalarOnKerr_guptt"));
  enablevar(grid, Ind("ScalarOnKerr_Gammattt"));

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
    double *cpsidot = vlldataptr(ucur, box, 1);
    double *ppsidot = vlldataptr(upre, box, 1);
    double *npsidot = vlldataptr(unew, box, 1);
    int i;
    int ipsi 	= (ucur)->index[0];
    int ipsidot = (ucur)->index[1];
    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
    double *psixx = box->v[Ind("ScalarOnKerr_ddpsixx")];
    double *psixy = box->v[Ind("ScalarOnKerr_ddpsixx")+1];
    double *psixz = box->v[Ind("ScalarOnKerr_ddpsixx")+2];
    double *psiyy = box->v[Ind("ScalarOnKerr_ddpsixx")+3];
    double *psiyz = box->v[Ind("ScalarOnKerr_ddpsixx")+4];
    double *psizz = box->v[Ind("ScalarOnKerr_ddpsixx")+5];
    double *psidotx = box->v[Ind("ScalarOnKerr_dpsidotx")];
    double *psidoty = box->v[Ind("ScalarOnKerr_dpsidotx")+1];
    double *psidotz = box->v[Ind("ScalarOnKerr_dpsidotx")+2];
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
    FirstDerivsOf_S(box, ipsidot , Ind("ScalarOnKerr_dpsidotx"));

    /* loop over points and set RHS */
    forallpoints(box, i)
    {
      double rpsidot, rpsi;
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double t = grid->time;
      double x0, y0;
      double g_ddpsi, gGt,gGx,gGy,gGz, gG_dpsi;
      /* g is upper metric */
      /* get all terms with less than 2 time derivs in g^ab d_a d_b psi */
      g_ddpsi = 2.0*(gtx[i]*psidotx[i] +gty[i]*psidoty[i] +gtz[i]*psidotz[i]) + 
                gxx[i]*psixx[i] + gyy[i]*psiyy[i] + gzz[i]*psizz[i] +
                2.0*(gxy[i]*psixy[i] + gxz[i]*psixz[i] + gyz[i]*psiyz[i]);

      /* get gG[a] = g^bc Gamma^a_bc */
      gGt = gxx[i]*Gamtxx[i] + gyy[i]*Gamtyy[i] + gzz[i]*Gamtzz[i] +
            2.0*(gxy[i]*Gamtxy[i] + gxz[i]*Gamtxz[i] + gyz[i]*Gamtyz[i]);
      gGx = gxx[i]*Gamxxx[i] + gyy[i]*Gamxyy[i] + gzz[i]*Gamxzz[i] +
            2.0*(gxy[i]*Gamxxy[i] + gxz[i]*Gamxxz[i] + gyz[i]*Gamxyz[i]);
      gGy = gxx[i]*Gamyxx[i] + gyy[i]*Gamyyy[i] + gzz[i]*Gamyzz[i] +
            2.0*(gxy[i]*Gamyxy[i] + gxz[i]*Gamyxz[i] + gyz[i]*Gamyyz[i]);
      gGz = gxx[i]*Gamzxx[i] + gyy[i]*Gamzyy[i] + gzz[i]*Gamzzz[i] +
            2.0*(gxy[i]*Gamzxy[i] + gxz[i]*Gamzxz[i] + gyz[i]*Gamzyz[i]);
      gG_dpsi = gGt*cpsidot[i] + gGx*psix[i] + gGy*psiy[i] + gGz*psiz[i];

      /* source posistion */
      x0 = 10*cos(0.02*t);
      y0 = 10*sin(0.02*t);

      /* set RHS of psi and psidot */
      rpsidot = -(g_ddpsi - gG_dpsi)/gtt[i] + 
                exp(-(x-x0)*(x-x0))*exp(-(y-y0)*(y-y0))*exp(-z*z); // source
      rpsi    = cpsidot[i];

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        npsidot[i] = ppsidot[i] + dt * rpsidot;
        npsi[i]    = ppsi[i]    + dt * rpsi;
      }
      else
      {
        npsidot[i] = rpsidot;
        npsi[i]    = rpsi;
      }
    }
  }

  /* set BCs */
  set_boundary(unew, upre, dt, ucur);

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}




/* compute and integrals of rho */
int ScalarOnKerr_analyze(tGrid *grid)
{
  int b;

// if( ! timeforoutput_index(grid, Ind("ScalarOnKerr_rho")) ) return 0;

// do something useful here 

  return 0;
}
