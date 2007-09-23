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

//  /* register alternative BC routine */
//  evolve_algebraicConditionsregister(set_boundary_ofPi);

  /* filter all newly computed vars */
  evolve_algebraicConditionsregister(filter_unew);

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
  enablevar(grid, Ind("ScalarOnKerr_Gt"));

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
       Ind("ScalarOnKerr_Gammattt"), Ind("ScalarOnKerr_Gt"));

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
    int i_G = Ind("ScalarOnKerr_Gt");
    double *Gt = box->v[i_G];
    double *Gx = box->v[i_G+1];
    double *Gy = box->v[i_G+2];
    double *Gz = box->v[i_G+3];

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
      double g_ddpsi, gG_dpsi;
      /* g is upper metric */
      /* get all terms with less than 2 time derivs in g^ab d_a d_b psi */
      g_ddpsi = -2.0*(gtx[i]*Pix[i] +gty[i]*Piy[i] +gtz[i]*Piz[i]) + 
                gxx[i]*psixx[i] + gyy[i]*psiyy[i] + gzz[i]*psizz[i] +
                2.0*(gxy[i]*psixy[i] + gxz[i]*psixz[i] + gyz[i]*psiyz[i]);

      /* get G^a dpsi_a, where G[a] = g^bc Gamma^a_bc */
      gG_dpsi = -Gt[i]*cPi[i] + Gx[i]*psix[i] + Gy[i]*psiy[i] + Gz[i]*psiz[i];

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
//  set_boundary(unew, upre, dt, ucur);
  set_psi_Pi_boundary(unew, upre, dt, ucur);

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
//    double *cpsi = vlldataptr(ucur, box, 0);
//    double *ppsi = vlldataptr(upre, box, 0);
//    double *npsi = vlldataptr(unew, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *pPi = vlldataptr(upre, box, 1);
    double *nPi = vlldataptr(unew, box, 1);
//    int ipsi = (ucur)->index[0];
    int iPi  = (ucur)->index[1];
//    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
//    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
//    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
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
    int i_G = Ind("ScalarOnKerr_Gt");
    double *Gt = box->v[i_G];
    double *Gx = box->v[i_G+1];
    double *Gy = box->v[i_G+2];
    double *Gz = box->v[i_G+3];

    /* compute the spatial derivs */
//    FirstDerivsOf_S(box, ipsi, Ind("ScalarOnKerr_dpsix"));
    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      double x,y,z;
      double r, nx,ny,nz;
      double rPi;
      double betan, alpha, alpha2, gnn, Gn, gdn;
      double ap,bp,cp, lambdap, dnPi;

      ijk = Index(i,j,k);
      x = px[ijk];
      y = py[ijk];
      z = pz[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r;
      ny = y/r;
      nz = z/r;
      alpha2 = -1.0/gtt[ijk];
      alpha  = sqrt(alpha2);
      betan = alpha2*(gtx[ijk]*nx +gty[ijk]*ny +gtz[ijk]*nz);
      gnn = gxx[ijk]*nx*nx +gyy[ijk]*ny*ny +gzz[ijk]*nz*nz + 
            2.0*(gxy[ijk]*nx*ny +gxz[ijk]*nx*nz +gyz[ijk]*ny*nz);
      Gn = Gx[ijk]*nx + Gy[ijk]*ny + Gz[ijk]*nz;
      /* dn = d_i n_j = (delta_ij - n_i n_j)/r */
      /* gdn = g^ij d_i n_j */
      gdn = (gxx[ijk]+gyy[ijk]+gzz[ijk]-gnn)/r;
      lambdap = betan + sqrt( betan*betan + alpha2*gnn);
      ap = lambdap/(alpha2*gnn);
      bp = 1.0;
      cp = (gdn-Gn)/gnn;
      dnPi = nx*Pix[ijk] + ny*Piy[ijk] + nz*Piz[ijk];
      
      /* set RHS of Pi */
      rPi = -(bp*dnPi + cp*cPi[ijk])/ap;

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
        nPi[ijk]  = pPi[ijk]  + dt * rPi;
      else
        nPi[ijk]  = rPi;
    }
  }
}

/* filter all newly computed vars */
void filter_unew(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;

  coordinateDependentFilter(unew);

  /* filter high freq. angular modes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, jf,kf, vi;

    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *temp1 = box->v[Ind("temp1")];

      kf=n3/3; kf*=2;
      jf=n2/3; jf*=2;
      spec_Coeffs(box, var, temp1);
      forallijk(i,j,k)
        if(k>kf || j>jf) temp1[Index(i,j,k)]=0.0;
      spec_Eval(box, var, temp1);
    }
  } /* end forallboxes */
}


/* compute things */
int ScalarOnKerr_analyze(tGrid *grid)
{
  int b;

// if( ! timeforoutput_index(grid, Ind("ScalarOnKerr_rho")) ) return 0;

// do something useful here 

  return 0;
}


/******************************************************************/
/* another version of scalar evo: */
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
