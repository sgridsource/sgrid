/* ScalarOnKerr.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "ScalarOnKerr.h"
#include "DV_CircSchwSource/Constants.h"
#include "DV_CircSchwSource/Source.h"


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

  VarNameSetBoundaryInfo("ScalarOnKerr_Up",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarOnKerr_Um",   0, 1, 1.0);

  /* create a variable list for ScalarOnKerr evolutions  */
  ScalarOnKerrvars = vlalloc(grid);
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_psi"));
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_Pi"));
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_Up"));
  vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_Um"));
  if(Getv("ScalarOnKerr_1stOrder_inSpace", "yes"))
  {
    VarNameSetBoundaryInfo("ScalarOnKerr_phil", 0, 1, 1.0);
    VarNameSetBoundaryInfo("ScalarOnKerr_phim", 0, 1, 1.0);
    VarNameSetBoundaryInfo("ScalarOnKerr_phix", 0, 1, 1.0);
    VarNameSetBoundaryInfo("ScalarOnKerr_phiy", 0, 1, 1.0);
    VarNameSetBoundaryInfo("ScalarOnKerr_phiz", 0, 1, 1.0);
    vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_phil"));
    vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_phim"));
    vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_phix"));
    vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_phiy"));
    vlpush(ScalarOnKerrvars, Ind("ScalarOnKerr_phiz"));
  }
  if (0) prvarlist(ScalarOnKerrvars);
  enablevarlist(ScalarOnKerrvars);

  /* register evolved variables */
  evolve_vlregister(ScalarOnKerrvars);
  
  /* register evolution routine */
  evolve_rhsregister(ScalarOnKerr_evolve);

//  /* register alternative BC routine */
//  evolve_algebraicConditionsregister(set_boundary_ofPi);

  /* filter all newly computed vars */
  if(Getv("ScalarOnKerr_filter_unew", "simple"))
    evolve_algebraicConditionsregister(filter_unew);
  else if(Getv("ScalarOnKerr_filter_unew", "naive_Ylm"))
    evolve_algebraicConditionsregister(naive_Ylm_filter_unew);

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

  /* enable char. vars. */
  //enablevar(grid, Ind("ScalarOnKerr_Up"));
  //enablevar(grid, Ind("ScalarOnKerr_Um"));

  /* enable rho */
  enablevar(grid, Ind("ScalarOnKerr_rho"));

  /* enable all derivative vars */
  enablevar(grid, Ind("ScalarOnKerr_dpsix"));
  enablevar(grid, Ind("ScalarOnKerr_ddpsixx"));
  enablevar(grid, Ind("ScalarOnKerr_dPix"));
  if(Getv("ScalarOnKerr_1stOrder_inSpace", "yes"))
    enablevar(grid, Ind("ScalarOnKerr_dphixx"));

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
  double t = ucur->time;
  int firstorder = Getv("ScalarOnKerr_1stOrder_inSpace", "yes");
  // double x0, y0;
  double M = Getd("BHmass");
  double r0, Omega, Dr, q, q22;

  ///* old source position */
  //x0 = 10*cos(0.02*t);
  //y0 = 10*sin(0.02*t);

  /* source parameters */
  q  = 1.0;
  r0 = 10.0*M;
  Dr = 1.0*M;
  q22= 4.0*PI*q/(r0*sqrt(1.0-3*M/r0))*sqrt(5.0/(96.0*PI))*3.0;
  Omega = sqrt(M/(r0*r0*r0));

/* call Ian's func to set mass and radius */
set_mass_radius(M,r0);
  
  /* loop over all boxes */
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
    double *dpsix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *dpsiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *dpsiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
    double *ddpsixx = box->v[Ind("ScalarOnKerr_ddpsixx")];
    double *ddpsixy = box->v[Ind("ScalarOnKerr_ddpsixx")+1];
    double *ddpsixz = box->v[Ind("ScalarOnKerr_ddpsixx")+2];
    double *ddpsiyy = box->v[Ind("ScalarOnKerr_ddpsixx")+3];
    double *ddpsiyz = box->v[Ind("ScalarOnKerr_ddpsixx")+4];
    double *ddpsizz = box->v[Ind("ScalarOnKerr_ddpsixx")+5];
    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
    double *cphix, *cphiy, *cphiz;
    double *pphix, *pphiy, *pphiz;
    double *nphix, *nphiy, *nphiz;
    double *dphixx, *dphixy, *dphixz;
    double *dphiyx, *dphiyy, *dphiyz;
    double *dphizx, *dphizy, *dphizz;
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
    if(firstorder)
    {
      cphix = vlldataptr(ucur, box, 6);
      cphiy = vlldataptr(ucur, box, 7);
      cphiz = vlldataptr(ucur, box, 8);
      pphix = vlldataptr(upre, box, 6);
      pphiy = vlldataptr(upre, box, 7);
      pphiz = vlldataptr(upre, box, 8);
      nphix = vlldataptr(unew, box, 6);
      nphiy = vlldataptr(unew, box, 7);
      nphiz = vlldataptr(unew, box, 8);
      dphixx = box->v[Ind("ScalarOnKerr_dphixx")];
      dphixy = box->v[Ind("ScalarOnKerr_dphixx")+1];
      dphixz = box->v[Ind("ScalarOnKerr_dphixx")+2];
      dphiyx = box->v[Ind("ScalarOnKerr_dphixx")+3];
      dphiyy = box->v[Ind("ScalarOnKerr_dphixx")+4];
      dphiyz = box->v[Ind("ScalarOnKerr_dphixx")+5];
      dphizx = box->v[Ind("ScalarOnKerr_dphixx")+6];
      dphizy = box->v[Ind("ScalarOnKerr_dphixx")+7];
      dphizz = box->v[Ind("ScalarOnKerr_dphixx")+8];
      cart_partials(box, cphix, dphixx,dphixy,dphixz);
      cart_partials(box, cphiy, dphiyx,dphiyy,dphiyz);
      cart_partials(box, cphiz, dphizx,dphizy,dphizz);
    }
    else
    {
      allDerivsOf_S(box, ipsi,
                    Ind("ScalarOnKerr_dpsix"), Ind("ScalarOnKerr_ddpsixx"));
    }
    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forallpoints(box, i)
    {
      double rPi, rpsi;
      double rphix, rphiy, rphiz;
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double rho, r, phi, theta, Y22;
      double g_ddpsi, gG_dpsi;  /* g is upper metric */
      double psix,psiy,psiz, psixx,psixy,psixz,psiyy,psiyz,psizz;
      double psidot, psidotx,psidoty,psidotz, psidotdot;

      /* set derivs*/
      if(firstorder)
      {
        psix = cphix[i];
        psiy = cphiy[i];
        psiz = cphix[i];
        psixx = dphixx[i];
        psixy = 0.5*(dphixy[i]+dphiyx[i]);
        psixz = 0.5*(dphixz[i]+dphizx[i]);
        psiyy = dphiyy[i];
        psiyz = 0.5*(dphiyz[i]+dphizy[i]);
        psizz = dphizz[i];
      }
      else
      {
        psix = dpsix[i];
        psiy = dpsiy[i];
        psiz = dpsiz[i];
        psixx = ddpsixx[i];
        psixy = ddpsixy[i];
        psixz = ddpsixz[i];
        psiyy = ddpsiyy[i];
        psiyz = ddpsiyz[i];
        psizz = ddpsizz[i];
      }
      /* set psidot and derivs from Pi */
      psidot  = cPi[i];
      psidotx = Pix[i];
      psidoty = Piy[i];
      psidotz = Piz[i];

      /* g is upper metric */
      /* get all terms with less than 2 time derivs in g^ab d_a d_b psi */
      g_ddpsi = 2.0*(gtx[i]*psidotx +gty[i]*psidoty +gtz[i]*psidotz) + 
                gxx[i]*psixx + gyy[i]*psiyy + gzz[i]*psizz +
                2.0*(gxy[i]*psixy + gxz[i]*psixz + gyz[i]*psiyz);

      /* get G^a dpsi_a, where G[a] = g^bc Gamma^a_bc */
      gG_dpsi = Gt[i]*psidot + Gx[i]*psix + Gy[i]*psiy + Gz[i]*psiz;

      /* source rho */
      r     = sqrt(x*x + y*y + z*z);
      theta = 0.5*PI - asin(z/r);
      phi   = Arg(x,y); // returns value in (-PI,PI]
      Y22 = sqrt(5.0/(96.0*PI))*1.5*(1.0 - cos(2.0*theta));
      rho = (q22/(4.0*PI*r0))*(exp( -(r-r0)*(r-r0)/(Dr*Dr) )/(sqrt(PI)*Dr))*
            cos(2.0*(Omega*t - phi)) * Y22;

/* use Ian's source */
//rho = SourceInKerrSchild(1.204119982655925 + t, x, y, z);

      /* compute psidotdot */
      psidotdot = -(g_ddpsi - gG_dpsi + 4.0*PI*rho)/gtt[i];
                   //-(1-Attenuation01( ((x-x0)*(x-x0)+(y-y0)*(y-y0)+z*z)/36,2,0.5)); // old source
                   //  -exp(-(x-x0)*(x-x0))*exp(-(y-y0)*(y-y0))*exp(-z*z); // oldest source
                
      /* set RHS of psi and Pi */
      rPi  = psidotdot;
      rpsi = psidot;
      /* set RHS of phix, ... if needed */
      if(firstorder) { rphix = psidotx;  rphiy = psidoty;  rphiz = psidotz; }

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        nPi[i]  = pPi[i]  + dt * rPi;
        npsi[i] = ppsi[i] + dt * rpsi;
        if(firstorder)
        { 
          nphix[i] = pphix[i] + dt * rphix;
          nphiy[i] = pphiy[i] + dt * rphiy;
          nphiz[i] = pphiz[i] + dt * rphiz;
        }
      }
      else /* if(dt==0.0) */
      {
        nPi[i]  = rPi;
        npsi[i] = rpsi;
        if(firstorder)
        { 
          nphix[i] = rphix;
          nphiy[i] = rphiy;
          nphiz[i] = rphiz;
        }
      } /* end else */
    }
  } /* end forallboxes */


  /* special nPi filter */
  if(!Getv("ScalarOnKerr_special_nPi_filter", "no"))
  {
    tVarList *vl_Pi = vlalloc(grid);
    vlpush(vl_Pi, unew->index[1]);
    if(Getv("ScalarOnKerr_special_nPi_filter", "simple"))
      filter_unew(vl_Pi, 0);
    else if(Getv("ScalarOnKerr_special_nPi_filter", "naive_Ylm"))
      naive_Ylm_filter_unew(vl_Pi, 0);
    vlfree(vl_Pi);
  }
//filter_unew_radially(unew, NULL);

  /* set char. vars. and use them to compute unew */
  set_Up_Um_onBoundary(unew, upre, dt, ucur);
  compute_unew_from_Up_Um_onBoundary(unew, upre, dt, ucur);
//interpolate_between_boxes(unew, NULL);

  /* set BCs */
////  set_boundary(unew, upre, dt, ucur);
//  set_psi_Pi_phi_boundary(unew, upre, dt, ucur);
  set_psi_Pi_boundary(unew, upre, dt, ucur);

  if(Getv("ScalarOnKerr_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}


/* set BC for old version: this is the 2007 version which sets only rPi */
//void set_psi_Pi_boundary_2007_PiVersion(tVarList *unew, tVarList *upre, double dt, 
void set_psi_Pi_boundary(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  
  /* not: forallboxes(grid,b) , use only outermost box */
  b = grid->nboxes - 1;
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
      double betan, alpha2, gnn, Gn, gdn;
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
      cp = 0.5*(gdn-Gn)/gnn;
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


/* similar to set_psi_Pi_boundary, but set BC on d_t psi instead of d_t Pi */
void set_psi_Pi_boundary_New(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  
  /* not: forallboxes(grid,b) , use only outermost box */
  b = grid->nboxes - 1;
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *cpsi = vlldataptr(ucur, box, 0);
    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
//    double *cPi = vlldataptr(ucur, box, 1);
//    double *pPi = vlldataptr(upre, box, 1);
//    double *nPi = vlldataptr(unew, box, 1);
    int ipsi = (ucur)->index[0];
//    int iPi  = (ucur)->index[1];
    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
//    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
//    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
//    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
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
    FirstDerivsOf_S(box, ipsi, Ind("ScalarOnKerr_dpsix"));
//    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      double x,y,z;
      double r, nx,ny,nz;
      double rpsi;
      double betan, alpha2, gnn, Gn, gdn;
      double ap,bp,cp, lambdap, dnpsi;

      ijk = Index(i,j,k);
      x = px[ijk];
      y = py[ijk];
      z = pz[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r;
      ny = y/r;
      nz = z/r;
      alpha2 = -1.0/gtt[ijk];
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
      cp = 0.5*(gdn-Gn)/gnn;
      dnpsi = nx*psix[ijk] + ny*psiy[ijk] + nz*psiz[ijk];
      
      /* set RHS of psi */
      rpsi = -(bp*dnpsi + cp*cpsi[ijk])/ap;

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
        npsi[ijk]  = ppsi[ijk]  + dt * rpsi;
      else
        npsi[ijk]  = rpsi;
    }
  }
}


/* set outer BC of psi, Pi, phi for first order version */
void set_psi_Pi_phi_boundary(tVarList *unew, tVarList *upre, double dt, 
                             tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;

  /* not: forallboxes(grid,b) , use only outermost box */
  b = grid->nboxes - 1;
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *cpsi = vlldataptr(ucur, box, 0);
//    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *pPi = vlldataptr(upre, box, 1);
    double *nPi = vlldataptr(unew, box, 1);
//    int ipsi = (ucur)->index[0];
//    int iPi  = (ucur)->index[1];
//    double *psix = box->v[Ind("ScalarOnKerr_dpsix")];
//    double *psiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
//    double *psiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
//    double *Pix = box->v[Ind("ScalarOnKerr_dPix")];
//    double *Piy = box->v[Ind("ScalarOnKerr_dPix")+1];
//    double *Piz = box->v[Ind("ScalarOnKerr_dPix")+2];
    double *cphix = vlldataptr(ucur, box, 6);
    double *cphiy = vlldataptr(ucur, box, 7);
    double *cphiz = vlldataptr(ucur, box, 8);
    double *pphix = vlldataptr(upre, box, 6);
    double *pphiy = vlldataptr(upre, box, 7);
    double *pphiz = vlldataptr(upre, box, 8);
    double *nphix = vlldataptr(unew, box, 6);
    double *nphiy = vlldataptr(unew, box, 7);
    double *nphiz = vlldataptr(unew, box, 8);
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
//    FirstDerivsOf_S(box, iPi , Ind("ScalarOnKerr_dPix"));

    /* loop over points and set RHS */
    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      double x,y,z;
      double r, nx,ny,nz, cos_th,sin_th, sin_ph,cos_ph, mx,my,mz, lx,ly,lz;
      double rPi, rphix,rphiy,rphiz;
      double betan, alpha2, gnn, Gn, gdn;
      double ap,bp,cp, lambdap, am,bm,cm, lambdam, nphin, cphin;
      double nphil,nphim, nUp,nUm, cUp,cUm;
     
      ijk = Index(i,j,k);
      x = px[ijk];
      y = py[ijk];
      z = pz[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r; /* nx=sin_th cos_ph, ny=sin_th sin_ph, nz=cos_th */
      ny = y/r;
      nz = z/r;
      cos_th = nz;
      sin_th = sqrt(1.0-cos_th*cos_th);
      if(sin_th>dequaleps) { cos_ph=nx/sin_th; sin_ph=ny/sin_th; }
      else { cos_ph=1.0; sin_ph=0.0; } 
      lx = cos_th*cos_ph;
      ly = cos_th*sin_ph;
      lz = -sin_th;
      mx = -sin_ph;
      my = cos_ph;
      mz = 0.0;
      alpha2 = -1.0/gtt[ijk];
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
      cp = 0.5*(gdn-Gn)/gnn;
      lambdam = betan - sqrt( betan*betan + alpha2*gnn );
      am = lambdam/(alpha2*gnn);
      bm = bp;
      cm = cp;

      nphin = nx*nphix[ijk] + ny*nphiy[ijk] + nz*nphiz[ijk];
      nphil = lx*nphix[ijk] + ly*nphiy[ijk] + lz*nphiz[ijk];
      nphim = mx*nphix[ijk] + my*nphiy[ijk] + mz*nphiz[ijk];

      /* nUp = ap*nPi[ijk] + bp*nphin + cp*npsi[ijk];
         nUm = am*nPi[ijk] + bm*nphin + cm*npsi[ijk]; */
      nUp = 0.0;  /* <-set d/dt Up = 0 */
      nUm = am*nPi[ijk] + bm*nphin + cm*npsi[ijk];

nUp = ap*nPi[ijk] + bp*nphin + cp*npsi[ijk];

cphin = nx*cphix[ijk] + ny*cphiy[ijk] + nz*cphiz[ijk];
cUp = ap*cPi[ijk] + bp*cphin + cp*cpsi[ijk];
cUm = am*cPi[ijk] + bm*cphin + cm*cpsi[ijk];
nUp += -160*(cUp-0);

      /* set RHS of Pi and phi */
      rPi = (nUp - nUm)/(ap-am);
//rPi = nUp - bp*nphin - cp*npsi[ijk];
//rPi=0.1;
      nphin = nUm - am*rPi -cm*npsi[ijk];
      rphix = nphin*nx + nphil*lx + nphim*mx;
      rphiy = nphin*ny + nphil*ly + nphim*my;
      rphiz = nphin*nz + nphil*lz + nphim*mz;

      /* set new vars or RHS, depending in which integrator is used */
      if(dt!=0.0)
      {
        nPi[ijk]   = pPi[ijk]   + dt * rPi;
        nphix[ijk] = pphix[ijk] + dt * rphix;
        nphiy[ijk] = pphiy[ijk] + dt * rphiy;
        nphiz[ijk] = pphiz[ijk] + dt * rphiz;
      }
      else
      {
        nPi[ijk]   = rPi;
        nphix[ijk] = rphix;
        nphiy[ijk] = rphiy;
        nphiz[ijk] = rphiz;
      }  /* end if(dt==0.0) */
    }
  }
}


/* compute and transfer Up and Um on boundaries */
void set_Up_Um_onBoundary(tVarList *unew, tVarList *upre, double dt, 
                          tVarList *ucur)
{
  tGrid *grid = unew->grid;
  int b;

  /* add my own "penalty" terms to npsi */
  for(b=1; b<grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    tBox *lbox = grid->box[b-1]; /* shell inside shell(=box) */
    int ln1=lbox->n1;
    int ln2=lbox->n2;
    int ijk, l_ijk, i,j,k;
    double *npsi =  vlldataptr(unew,  box, 0);
    double *lnpsi = vlldataptr(unew, lbox, 0);
    double *cpsi =  vlldataptr(ucur,  box, 0);
    double *lcpsi = vlldataptr(ucur, lbox, 0);
    double tau2 = 0.05/grid->dt; // 4.0;
//double npsival;
////double tau2 = 0.0025*ln1*(ln1+1)+ 0.0025*n1*(n1+1);
//double tau1=0.0;
//double tau2=0.1;
//double *cpsiX = box->v[Ind("temp1")];
//double *lcpsiX=lbox->v[Ind("temp1")];
//double cpsir, lcpsir;
//
//spec_Deriv1(box, 1, cpsi, cpsiX);
//if(b==1) spec_Deriv1(lbox, 1, lcpsi, lcpsiX);
    
    /* loop over boundary points */
    forplane1(i,j,k, n1,n2,n3, 0) /* assume that all boxes have same n2,n3 */
    {
      ijk  = Index(i,j,k);
      l_ijk= Ind_n1n2(ln1-1, j, k, ln1,ln2);
//npsival = 0.5*(npsi[ijk] + lnpsi[l_ijk]);
//lnpsi[l_ijk] = npsi[ijk] = npsival;
//cpsir = cpsiX[ijk];
//lcpsir=lcpsiX[l_ijk];
//lnpsi[l_ijk] -= tau1*(lcpsi[l_ijk] - cpsi[ijk]) + 
//                tau2*(lcpsir - cpsir);
//npsi[ijk]    -= tau1*(cpsi[ijk] - lcpsi[l_ijk]) +
//                tau2*(cpsir - lcpsir);
//      /* ??? assume ??? that psi on left should be equal to psi on right */
      lnpsi[l_ijk] -= tau2*(lcpsi[l_ijk] - cpsi[ijk]);
      npsi[ijk]    -= tau2*(cpsi[ijk] - lcpsi[l_ijk]);
    }
  } /* end for b */

  /* compute nUp and nUm in each box */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int ijk, pi;
    double *npsi = vlldataptr(unew, box, 0);
    double *nPi = vlldataptr(unew, box, 1);
    double *nUp = vlldataptr(unew, box, 2);
    double *nUm = vlldataptr(unew, box, 3);
    double *npsix = box->v[Ind("temp1")];
    double *npsiy = box->v[Ind("temp2")];
    double *npsiz = box->v[Ind("temp3")];
    double *cpsi = vlldataptr(ucur, box, 0);
    double *cPi = vlldataptr(ucur, box, 1);
    double *cUp = vlldataptr(ucur, box, 2);
    double *cUm = vlldataptr(ucur, box, 3);
    double *cpsix = box->v[Ind("ScalarOnKerr_dpsix")];
    double *cpsiy = box->v[Ind("ScalarOnKerr_dpsix")+1];
    double *cpsiz = box->v[Ind("ScalarOnKerr_dpsix")+2];
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

    /* compute spatial derivs of npsi */
    cart_partials(box, npsi, npsix, npsiy, npsiz);
    cart_partials(box, cpsi, cpsix, cpsiy, cpsiz);

    /* loop over points and set RHS */
    forPointList_inbox(boxBoundaryPointList, box, pi , ijk)
    {
      double x,y,z;
      double r, nx,ny,nz;
      double rPi;
      double betan, alpha2, gnn, Gn, gdn;
      double ap,bp,cp, lambdap, am,bm,cm, lambdam, npsir, cpsir;

      x = px[ijk];
      y = py[ijk];
      z = pz[ijk];
      r = sqrt(x*x + y*y + z*z);
      nx = x/r;
      ny = y/r;
      nz = z/r;
      alpha2 = -1.0/gtt[ijk];
      betan = alpha2*(gtx[ijk]*nx +gty[ijk]*ny +gtz[ijk]*nz);
      gnn = gxx[ijk]*nx*nx +gyy[ijk]*ny*ny +gzz[ijk]*nz*nz + 
            2.0*(gxy[ijk]*nx*ny +gxz[ijk]*nx*nz +gyz[ijk]*ny*nz);
      Gn = Gx[ijk]*nx + Gy[ijk]*ny + Gz[ijk]*nz;
      /* dn = d_i n_j = (delta_ij - n_i n_j)/r */
      /* gdn = g^ij d_i n_j */
      gdn = (gxx[ijk]+gyy[ijk]+gzz[ijk]-gnn)/r;
      lambdap = betan + sqrt( betan*betan + alpha2*gnn );
      ap = lambdap/(alpha2*gnn);
      bp = 1.0;
      cp = 0.5*(gdn-Gn)/gnn;
      lambdam = betan - sqrt( betan*betan + alpha2*gnn );
      am = lambdam/(alpha2*gnn);
      bm = bp;
      cm = cp;

      npsir = ( nx*npsix[ijk] + ny*npsiy[ijk] + nz*npsiz[ijk] );
      cpsir = ( nx*cpsix[ijk] + ny*cpsiy[ijk] + nz*cpsiz[ijk] );

      /* set char vars both unew and ucur */
      nUp[ijk] = ap*nPi[ijk] + bp*npsir + cp*npsi[ijk];
      nUm[ijk] = am*nPi[ijk] + bm*npsir + cm*npsi[ijk];
      cUp[ijk] = ap*cPi[ijk] + bp*cpsir + cp*cpsi[ijk];
      cUm[ijk] = am*cPi[ijk] + bm*cpsir + cm*cpsi[ijk];
//if(ijk==20 || ijk==39)
//{
//printf("b=%d ijk=%d (x,y,z)=(%g,%g,%g): nUp[ijk]=%g nUm[ijk]=%g\n",
//b,ijk,x,y,z, nUp[ijk], nUm[ijk]);
//}
    }
  }

  /* add penalty terms to nUp and nUm on inner boundaries */
  for(b=1; b<grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    tBox *lbox = grid->box[b-1]; /* shell inside shell(=box) */
    int ln1=lbox->n1;
    int ln2=lbox->n2;
    int ijk, l_ijk, i,j,k;
    double *nUp =  vlldataptr(unew, box, 2);
    double *nUm =  vlldataptr(unew, box, 3);
    double *lnUp = vlldataptr(unew, lbox, 2);
    double *lnUm = vlldataptr(unew, lbox, 3);
    double *cUp =  vlldataptr(ucur, box, 2);
    double *cUm =  vlldataptr(ucur, box, 3);
    double *lcUp = vlldataptr(ucur, lbox, 2);
    double *lcUm = vlldataptr(ucur, lbox, 3);
    double *npsi =  vlldataptr(unew,  box, 0);
    double *lnpsi = vlldataptr(unew, lbox, 0);
    double *cpsi =  vlldataptr(ucur,  box, 0);
    double *lcpsi = vlldataptr(ucur, lbox, 0);
    // double ltau = 0.2*ln1*(ln1+1);
    // double tau  = 0.2*n1*(n1+1);
    double tau = 0.5/grid->dt;
//double tau = 800;

    /* loop over boundary points */
    forplane1(i,j,k, n1,n2,n3, 0) /* assume that all boxes have same n2,n3 */
    {
      ijk  = Index(i,j,k);
      l_ijk= Ind_n1n2(ln1-1, j, k, ln1,ln2);
      lnUp[l_ijk] -= tau*(lcUp[l_ijk] - cUp[ijk]); 
      nUm[ijk]    -= tau*(cUm[ijk] - lcUm[l_ijk]);
//      lnpsi[l_ijk] -= tau2*(lcpsi[l_ijk] - cpsi[ijk]);
//      npsi[ijk]    -= tau2*(cpsi[ijk] - lcpsi[l_ijk]);
//if(ijk==20)
//{
//printf("lb=%d l_ijk=%d: lnUp[l_ijk]=%g lnUm[l_ijk]=%g\n"
//       " b=%d   ijk=%d:  nUp[ijk]  =%g  nUm[ijk]=  %g\n",
//b-1,l_ijk, lnUp[l_ijk],lnUm[l_ijk], b,ijk, nUp[ijk],nUm[ijk]);
//}
    }
  } /* end for b */

  /* add penalty term to nUp on outer boundary */
  b=grid->nboxes-1;
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *nUp =  vlldataptr(unew, box, 2);
    double *cUp =  vlldataptr(ucur, box, 2);
    double tau = 0.5/grid->dt;

    /* loop over boundary points */
    forplane1(i,j,k, n1,n2,n3, n1-1) /* assume that all boxes have same n2,n3 */
    {
      ijk  = Index(i,j,k);
      nUp[ijk] -= tau*(cUp[ijk] - 0.0);
    }
  } /* end outer boundary penalty terms */
}

/* compute unew from Up and Um on boundaries */
void compute_unew_from_Up_Um_onBoundary(tVarList *unew, tVarList *upre,
                                        double dt, tVarList *ucur)
{
  tGrid *grid = unew->grid;
  int b;
  
  /* compute nUp and nUm in each box */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int ijk, i,j,k;
    double *npsi = vlldataptr(unew, box, 0);
    double *nPi = vlldataptr(unew, box, 1);
    double *nUp = vlldataptr(unew, box, 2);
    double *nUm = vlldataptr(unew, box, 3);
    double *npsiX = box->v[Ind("temp1")];
//    double *npsix = box->v[Ind("temp1")];
//    double *npsiy = box->v[Ind("temp2")];
//    double *npsiz = box->v[Ind("temp3")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *dXdx = box->v[Ind("dXdx")];
    double *dXdy = box->v[Ind("dXdx")+1];
    double *dXdz = box->v[Ind("dXdx")+2];
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

    /* compute spatial derivs of npsi */
    spec_Deriv1(box, 1, npsi, npsiX);
    /* move this-^ inside the for(i=0; i<n1; i+=n1-1) loop if npsi is 
       changed for each i */

    /* loop over lower and upper boundary points */
    for(i=0; i<n1; i+=n1-1)
    {
      if(b==0 && i==0) continue; /* do nothing on inner bound of box0 */
      if(b==grid->nboxes-1 && i==n1-1) break; /*do nothing on outer bound of last box */
      for(k=0; k<n3; k++)
      for(j=0; j<n2; j++)
      {
        double x,y,z;
        double r, nx,ny,nz;
        double rPi;
        double betan, alpha2, gnn, Gn, gdn;
        double ap,bp,cp, lambdap, am,bm,cm, lambdam;
        double dXdr, npsir, D00, npsir_plus_cnpsi; //, D00, npsir_plus_cnpsi, bnpsir;

        ijk = Index(i,j,k);
        x = px[ijk];
        y = py[ijk];
        z = pz[ijk];
        r = sqrt(x*x + y*y + z*z);
        nx = x/r;
        ny = y/r;
        nz = z/r;
        alpha2 = -1.0/gtt[ijk];
        betan = alpha2*(gtx[ijk]*nx +gty[ijk]*ny +gtz[ijk]*nz);
        gnn = gxx[ijk]*nx*nx +gyy[ijk]*ny*ny +gzz[ijk]*nz*nz + 
              2.0*(gxy[ijk]*nx*ny +gxz[ijk]*nx*nz +gyz[ijk]*ny*nz);
        Gn = Gx[ijk]*nx + Gy[ijk]*ny + Gz[ijk]*nz;
        /* dn = d_i n_j = (delta_ij - n_i n_j)/r */
        /* gdn = g^ij d_i n_j */
        gdn = (gxx[ijk]+gyy[ijk]+gzz[ijk]-gnn)/r;
        lambdap = betan + sqrt( betan*betan + alpha2*gnn );
        ap = lambdap/(alpha2*gnn);
        bp = 1.0;
        cp = 0.5*(gdn-Gn)/gnn;
        lambdam = betan - sqrt( betan*betan + alpha2*gnn );
        am = lambdam/(alpha2*gnn);
        bm = bp;
        cm = cp;

        /* compute dXdr = dXdx dxdr +... */
        dXdr = dXdx[ijk]*nx + dXdy[ijk]*ny + dXdz[ijk]*nz;
        /* npsir := d/dr psi = dXdr d/dX psi */
        npsir = dXdr * npsiX[ijk];
        /* D00 = dXdr * box->D1[(i+1)*(i+1)-1]; */

        /* compute nPi and from nUp or nUm on boundary */
        /* nUp[ijk] = ap*nPi[ijk] + bp*npsir + cp*npsi[ijk];
           nUm[ijk] = am*nPi[ijk] + bm*npsir + cm*npsi[ijk]; */
        nPi[ijk] = (nUp[ijk] - nUm[ijk])/(ap-am);
//        /* Note: if only nUp/m was correct at upper/lower bounday: */
//        if(i==0)  nPi[ijk] = (nUm[ijk] - bm*npsir - cm*npsi[ijk])/am;
//        else      nPi[ijk] = (nUp[ijk] - bp*npsir - cp*npsi[ijk])/ap;

//        /* Note: U0 = psi is a zero speed mode, so we might want to 
//           leave psi untouched!!! 
//           I.E.: we do not use the following : */
        /* U0 = psi is a zero speed mode, but really we need to set npsir */
        /* deriv = D00 u[0] + sum_{j=1...n-1} D0j u[j] */
        /* npsir_plus_cnpsi = D00 psi[0] + sum_{j=1...n-1} D0j psi[j] +
                              c*psi[0]
                            = c*psi[0] + D00 psi[0] + (D0j psi[j]-D00 psi[0]) 
         ==> psi[0] = (npsir_plus_cnpsi - (D0j psi[j]-D00 psi[0]))/(c+D00);
           npsi[ijk] = (npsir_plus_cnpsi - (npsir-D00*npsi[ijk]))/(cm+D00); */
        npsir_plus_cnpsi = nUm[ijk] - am*nPi[ijk];
        D00 = dXdr * box->D1[(i+1)*(i+1)-1];
        npsi[ijk] = (npsir_plus_cnpsi - (npsir-D00*npsi[ijk]))/(cm+D00);
// change npsi[ijk +/- 1] instead of npsi[ijk]
//double D01;
//int ijk_in, i_in;
//npsir_plus_cnpsi = nUm[ijk] - am*nPi[ijk];
//if(i==0) { i_in = 1;         ijk_in = ijk+1; }
//else     { i_in = n1*i+n1-2; ijk_in = ijk-1; }
//D01 = box->D1[i_in];// <-test: correct this with the dXdr at the right place
//npsi[ijk_in] = (npsir_plus_cnpsi - (npsir-D01*npsi[ijk_in]))/(cm+D01);
      }
    }
  }

  /* average npsi and nPi between boxes */
  for(b=1; b<grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    tBox *lbox = grid->box[b-1]; /* shell inside shell(=box) */
    int ln1=lbox->n1;
    int ln2=lbox->n2;
    int ijk, l_ijk, i,j,k;
    double *npsi  = vlldataptr(unew, box, 0);
    double *lnpsi = vlldataptr(unew, lbox, 0);
    double *nPi  = vlldataptr(unew, box, 1);
    double *lnPi = vlldataptr(unew, lbox, 1);
    double npsival, nPival;

    /* loop over boundary points */
    forplane1(i,j,k, n1,n2,n3, 0) /* assume that all boxes have same n1,n2,n3 */
    {
      ijk  = Index(i,j,k);
      l_ijk= Ind_n1n2(ln1-1, j, k, ln1,ln2);
      npsival = 0.5*(npsi[ijk] + lnpsi[l_ijk]);
      nPival  = 0.5*(nPi[ijk] + lnPi[l_ijk]);
//npsival=lnpsi[l_ijk];
//npsival=npsi[ijk];
//nPival= lnPi[l_ijk];
//      npsi[ijk] = lnpsi[l_ijk] = npsival;
//      nPi[ijk]  = lnPi[l_ijk] = nPival;
//if(j==1 && k==0 )
//{
//printf("#lb=%d l_ijk=%d: lnPi[l_ijk]=%g lnpsi[l_ijk]=%g\n"
//       "# b=%d   ijk=%d:  nPi[ijk]  =%g  npsi[ijk]=  %g\n",
//       b-1,l_ijk, lnPi[l_ijk], lnpsi[l_ijk], b,ijk, nPi[ijk], npsi[ijk]);
//}
    }
  } /* end for b */
}

/* interpolate between boxes */
void interpolate_between_boxes(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;
  
  /* loop over boxes */
  for(b=1; b<grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    tBox *lbox = grid->box[b-1]; /* shell inside shell(=box) */
    int ln1=lbox->n1;
    int ln2=lbox->n2;
    int ijk, l_ijk, i,j,k;
    double *npsi = vlldataptr(unew,  box, 0);
    double *lnpsi= vlldataptr(unew, lbox, 0);
    double *nPi  = vlldataptr(unew,  box, 1);
    double *lnPi = vlldataptr(unew, lbox, 1);
    double *X  =  box->v[Ind("X")];
    double *lX = lbox->v[Ind("X")];
    double *BM;
    double *lBM;

    /* set BM and lBM if necessary */
    BM = lBM = NULL;
    if( !dequal(X[0],lX[ln1 - 2]) )
    {
      lBM = (double *) calloc(ln1, sizeof(double));
      spec_Basis_times_CoeffMatrix_direc(lbox, 1, lBM, X[0]);
    }
    if( !dequal(lX[ln1-1],X[1]) )
    {
      BM = (double *) calloc(n1, sizeof(double));
      spec_Basis_times_CoeffMatrix_direc(box, 1, BM, lX[ln1-1]);
    }

    /* loop over boundary points */
    forplane1(i,j,k, n1,n2,n3, 0) /* assume that all boxes have same n2,n3 */
    {
      ijk  = Index(i,j,k);
      l_ijk= Ind_n1n2(ln1-2, j, k, ln1,ln2); 
      /* set values at i=0 in box from values in lbox */
      if( dequal(X[ijk],lX[l_ijk]) )
      {
        npsi[ijk]= lnpsi[l_ijk];
        nPi[ijk] = lnPi[l_ijk];
      }
      else
      {
        npsi[ijk]=scalarproduct_vectors(lBM, lnpsi +Ind_n1n2(0,j,k, ln1,ln2), ln1);
        nPi[ijk]=scalarproduct_vectors(lBM, lnPi +Ind_n1n2(0,j,k, ln1,ln2), ln1);
      }
      /* set values at i=ln1-1 in lbox from values in box */
      if( dequal(lX[l_ijk+1],X[ijk+1]) )
      {
        lnpsi[l_ijk+1]= npsi[ijk+1];
        lnPi[l_ijk+1] = nPi[ijk+1];
      }
      else
      {
        lnpsi[l_ijk+1]=scalarproduct_vectors(BM, npsi +Index(0,j,k), n1);
        lnPi[l_ijk+1] =scalarproduct_vectors(BM, nPi + Index(0,j,k), n1);
      }
    }
    free(lBM);
    free(BM);
  } /* end for b */
}


/* filter all newly computed vars */
void filter_unew(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;
  double n2frac = Getd("ScalarOnKerr_filter_n2frac");
  double n3frac = Getd("ScalarOnKerr_filter_n3frac");
  int shift2 = Geti("ScalarOnKerr_filter_shift2");
  int shift3 = Geti("ScalarOnKerr_filter_shift3");
  int ellipse= Getv("ScalarOnKerr_filter_YZregion", "ellipse");

  coordinateDependentFilter(unew);

  /* filter high freq. angular modes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, jf,kf, vi;
    double jmax, kmax;

    jf=0.5*n2frac*n2; jf*=2; jf+=shift2-2;
    jmax=jf;
        
    kf=0.5*n3frac*n3; kf*=2; kf+=shift3-2;
    kmax=kf;

    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *temp1 = box->v[Ind("temp1")];

      spec_Coeffs(box, var, temp1);
      forallijk(i,j,k)
      {
        int tj = 2*((j+1)/2);
        int tk = 2*((k+1)/2);

        if(ellipse) { if((tj/jmax)*(tj/jmax) + (tk/kmax)*(tk/kmax) > 1.0)
                        temp1[Index(i,j,k)]=0.0; }
        else        if(k>kf || j>jf) temp1[Index(i,j,k)]=0.0;
          
      }
      spec_Eval(box, var, temp1);
    }
  } /* end forallboxes */
}

/* filter all newly computed vars */
void naive_Ylm_filter_unew(tVarList *unew, tVarList *upre)
{
  Naive_YlmFilter(unew);
}

/* filter all newly computed vars */
void filter_unew_radially(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;

  /* filter high freq. angular modes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, f, vi;

    f=3*n1/2;

    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *temp1 = box->v[Ind("temp1")];

      spec_analysis1(box, 1, box->Mcoeffs1, var, temp1);
      forallijk(i,j,k) if(i>f) temp1[Index(i,j,k)]=0.0;
      spec_synthesis1(box, 1, box->Meval1, var, temp1);
    }
  } /* end forallboxes */
}


/* compute things */
int ScalarOnKerr_analyze(tGrid *grid)
{
  int b;

  if( ! timeforoutput_index(grid, Ind("ScalarOnKerr_rho")) ) return 0;

  /* set rho on grid */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *rho = box->v[Ind("ScalarOnKerr_rho")];

    forallpoints(box,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double t = grid->time;

      rho[i] = SourceInKerrSchild(1.204119982655925 + t, x, y, z);
    }
  }
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
