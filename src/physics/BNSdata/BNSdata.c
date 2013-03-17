/* BNSdata.c */
/* Wolfgang Tichy 2007 */


#include "sgrid.h"
#include "BNSdata.h"

#define Power pow

/* struct types used in root finder newton_linesrch_itsP */
typedef struct T_grid_x1_2_struct {
  tGrid *grid;   /* grid */
  double x1;  /* edge of NS1 */
  double x2;  /* edge of NS2 */
  int    it;  /* iteration number */
  double tol; /* tolerance */
} t_grid_x1_2_struct;

typedef struct T_grid_bXYZ1_bXYZ2_struct {
  tGrid *grid; /* grid */
  int b1;      /* box1 */
  double X1 ;  /* X pos 1 */
  double Y1 ;  /* Y pos 1 */
  double Z1 ;  /* Z pos 1 */
  int b2;      /* box2 */
  double X2 ;  /* X pos 2 */
  double Y2 ;  /* Y pos 2 */
  double Z2 ;  /* Z pos 2 */
} t_grid_bXYZ1_bXYZ2_struct;

typedef struct T_grid_grid0_m01_m02_struct {
  tGrid *grid;  /* grid where we operate */
  tGrid *grid0; /* grid from which we interpolate vars when domains change */
  double m01;   /* desired m01 */
  double m02;   /* desired m02 */
} t_grid_grid0_m01_m02_struct;


/* global vars */
extern double rf_surf1; /* radius of star1 */
extern double rf_surf2; /* radius of star2 */
extern double P_core1;  /* core pressure of star1 */
extern double P_core2;  /* core pressure of star2 */
extern tParameter *pdb;
extern int npdb;
extern int npdbmax;
tGrid *central_q_errors_VectorFunc__grid; /* grid var for central_q_errors_VectorFunc */
tGrid *xouts_error_VectorFunc__grid; /* grid var for xouts_error_VectorFunc */
int    xouts_error_VectorFunc__it;   /* it for xouts_error_VectorFunc */
double xouts_error_VectorFunc__tol;  /* tol for xouts_error_VectorFunc */
double xouts_error_VectorFunc__xout1;  /* xout1 we currently try to achieve */
double xouts_error_VectorFunc__xout2;  /* xout2 we currently try to achieve */
tGrid *xmaxs_error_VectorFunc__grid; /* grid var for xmaxs_error_VectorFunc */
double xmaxs_error_VectorFunc__xmax1;  /* xmax1 we currently try to achieve */
double xmaxs_error_VectorFunc__xmax2;  /* xmax2 we currently try to achieve */
tGrid *dqdx_at_Xmax1_2_VectorFunc__grid;  /* grid var for dqdx_at_Xmax1_2_VectorFunc */
int dqdx_at_Xmax1_2_VectorFunc__bi1; /* boxind of max1 */
int dqdx_at_Xmax1_2_VectorFunc__bi2; /* boxind of max2 */
double dqdx_at_Xmax1_2_VectorFunc__Xmax1; /* pos. of max1 */
double dqdx_at_Xmax1_2_VectorFunc__Ymax1; /* pos. of max1 */
double dqdx_at_Xmax1_2_VectorFunc__Xmax2; /* pos. of max2 */
double dqdx_at_Xmax1_2_VectorFunc__Ymax2; /* pos. of max2 */
              


/* global var lists */
tVarList *vlu, *vlFu, *vluDerivs;
tVarList *vldu, *vlJdu, *vlduDerivs;


/* functions in this file */
void save_surfacepos_in_BNSdata_surface_sigma_pm(int star, tGrid *grid,
                                                 tGrid *grid_bak, double t);
void compute_ABphi_from_xyz(tBox *box, double *A, double *B, double *phi,
                            double x, double y, double z);
void make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(tGrid *grid,
     tVarList **vlw,  tVarList **vlwDerivs,  tVarList **vlFw, 
     tVarList **vldw, tVarList **vldwDerivs, tVarList **vlJdw, char *Name);
void free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(
     tVarList *vlw,  tVarList *vlwDerivs,  tVarList *vlFw,
     tVarList *vldw, tVarList *vldwDerivs, tVarList *vlJdw);
double GridL2Norm_of_vars_in_string(tGrid *grid, char *str);
double GridL2Norm_of_vars_in_string_withZeroErr_inbox12(tGrid *grid, char *str);
int BNS_Eqn_Iterator_for_vars_in_string(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr, char *str);
int BNS_ordered_Eqn_Iterator(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
int BNS_Eqn_Iterator(tGrid *grid, int itmax, double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
int BNS_Eqn_sequence1(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
int BNS_Eqn_sequence2(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
int BNS_Eqn_sequence3(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
double GetInnerRestMass(tGrid *grid, int bi);
void m01_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m02_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m0_errors_VectorFuncP(int n, double *vec, double *fvec, void *p);
void compute_new_q_and_adjust_domainshapes_InterpFromGrid0(tGrid *grid, 
                                                           tGrid *grid0, 
                                                           int innerdom);
void compute_new_q_and_adjust_domainshapes(tGrid *grid, int innerdom);
void m01_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m02_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
void find_Varmax_along_x_axis_usingBNSdata_temp123(tGrid *grid, int varind, 
                                                int *bi, double *X, double *Y);
void find_qmax1_along_x_axis(tGrid *grid, int *bi, double *X, double *Y);
void central_q_errors_VectorFunc(int n, double *vec, double *fvec);
void estimate_q_errors_VectorFunc(int n, double *vec, double *fvec);
double BNSdata_find_position_of_qmax(tGrid *grid, int *bi, 
                                     double *X, double *Y, double *Z);
double BNSdata_find_xyz_of_qmax(tGrid *grid, int *bi, 
                                double *x, double *y, double *z);
void set_BNSdata_actual_xyzmax_pars(tGrid *grid);





/* set conformal rotational velocity for one star */
void BNS_set_wB(tGrid *grid, int star, double xc,double yc,double zc)
{
  int corot1 = Getv("BNSdata_rotationstate1","corotation");
  int corot2 = Getv("BNSdata_rotationstate2","corotation");
  int wBfac_Psi6    = Getv("BNSdata_wB_factor","Psi6");
  int wBfac_h       = Getv("BNSdata_wB_factor","h");
  int wBfac_ooalpha = Getv("BNSdata_wB_factor","1/alpha");
  int b;
  double BNSdata_n = Getd("BNSdata_n");
  double omegax1   = Getd("BNSdata_omegax1");
  double omegay1   = Getd("BNSdata_omegay1");
  double omegaz1   = Getd("BNSdata_omegaz1");
  double omegax2   = Getd("BNSdata_omegax2");
  double omegay2   = Getd("BNSdata_omegay2");
  double omegaz2   = Getd("BNSdata_omegaz2");
  int corot;
  double omegax, omegay, omegaz;

  if(star==1)
  {
    corot=corot1;
    omegax=omegax1;  omegay=omegay1;  omegaz=omegaz1;
  }
  else if(star==2)
  {
    corot=corot2;
    omegax=omegax2;  omegay=omegay2;  omegaz=omegaz2;
  }
  else
    errorexit("BNS_set_wB: star must be 1 or 2.");

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
    double *BNSdata_wBx = box->v[Ind("BNSdata_wBx")];
    double *BNSdata_wBy = box->v[Ind("BNSdata_wBx")+1];
    double *BNSdata_wBz = box->v[Ind("BNSdata_wBx")+2];
    double *BNSdata_Psi = box->v[Ind("BNSdata_Psi")];
    double *BNSdata_alphaP = box->v[Ind("BNSdata_alphaP")];
    double *BNSdata_q = box->v[Ind("BNSdata_q")];

    if( star==2 && (b==0 || b==1 || b==5) ) continue;
    if( star==1 && (b==3 || b==2 || b==4) ) continue;
    
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

      /* set wB */
      if(!corot)
      {
        double Psi1 = BNSdata_Psi[i];
        double Psi4 = Psi1*Psi1*Psi1*Psi1;
        double Psi6 = Psi1*Psi1*Psi4;
        double ooalpha = Psi1/BNSdata_alphaP[i]; /* 1/alpha */
        double h = (BNSdata_n+1.0)*BNSdata_q[i] + 1.0; /* h = (n+1) q + 1 */
        // double u0;
        double vx,vy,vz;
        double Att1, wBfac;
        double A = pX[i];

        if(b==1 || b==2)  Att1 = 1.0-Attenuation01((A-0.1)/0.8, 2.0, 0.5);
        else              Att1 = 1.0;

        /* omega cross r-rc */
        vx = ( omegay1* (z-zc) - omegaz1* (y-yc) )*Att1;
        vy = ( omegaz1* (x-xc) - omegax1* (z-zc) )*Att1;
        vz = ( omegax1* (y-yc) - omegay1* (x-xc) )*Att1;
        /* 1/u0^2 ~ alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),*/
        /* u0 = 1.0/sqrt(fabs(exp(2*Phi1) - Psi4*(vx*vx + vy*vy + vz*vz))); */

        /* set wBfac */
        /* wBfac = h*u0 * Psi^6  ???  */
        wBfac = 1.0;
        if(wBfac_Psi6)    wBfac = wBfac * Psi6;
        if(wBfac_h)       wBfac = wBfac * h;
        if(wBfac_ooalpha) wBfac = wBfac * ooalpha; /* u0 ~ 1/alpha */

        BNSdata_wBx[i] = vx * wBfac;
        BNSdata_wBy[i] = vy * wBfac;
        BNSdata_wBz[i] = vz * wBfac;
      }
      else
      {
        BNSdata_wBx[i] = 0.0; 
        BNSdata_wBy[i] = 0.0;
        BNSdata_wBz[i] = 0.0;
      }
    } /* end forallpoints */
  }
}


/* initial shift as in gr-qc/0007028v2*/
void BNSdata_initial_shift(int star, double fac, 
                           double m1, double m2, double Omega, double r12,
                           double rs, double x, double y, double z,
                           double *Bx, double *By, double *Bz)
{
  double F, epsW, epsXi, epsxy, Wy, dWydx,dWydy,dWydz, dXidx,dXidy,dXidz;
  double r = sqrt(x*x + y*y + z*z);
  double nx = x;
  double ny = y;
  double nz = z;
  double xz_on = 0.0;

  /* vector n = x/r */
  if(r>0) { nx=x/r; ny=y/r; nz=z/r; }

  /* function F */
  if(star==1) { F=fac*m1*Omega*r12/(1.0+m1/m2); epsW=+1; epsXi=+1; epsxy=-1; }
  else        { F=fac*m2*Omega*r12/(1.0+m2/m1); epsW=+1; epsXi=+1; epsxy=+1; }

  /* vector W , Wx = Wz = 0.0. And derivs of W and derivs of Xi */
  if(r<rs)
  {
    Wy = (6.0*F/rs)*( 1.0 - r*r/(3*rs*rs) )*epsW;
    dWydx = ( -(4.0*F/(rs*rs))*(r/rs)*nx )*epsW;
    dWydy = ( -(4.0*F/(rs*rs))*(r/rs)*ny )*epsW;
    dWydz = ( -(4.0*F/(rs*rs))*(r/rs)*nz )*epsW;
    dXidx = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * nx )*epsXi;
    dXidy = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * ny
              +((2.0*F)/rs) * (1.0 - 3.0*r*r/(5.0*rs*rs)) )*epsXi;
    dXidz = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * nz )*epsXi;
  }
  else
  {
    Wy = (4.0*F/r)*epsW;
    dWydx = ( (-4.0*F/(r*r))*nx )*epsW;    
    dWydy = ( (-4.0*F/(r*r))*ny )*epsW;    
    dWydz = ( (-4.0*F/(r*r))*nz )*epsW;    
    dXidx = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*nx )*epsXi;
    dXidy = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*ny
              +((4.0*F)/5.0) * ((rs*rs)/(r*r*r)) )*epsXi;
    dXidz = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*nz )*epsXi;
  }

  /* shift */
  *Bx = (          - 0.125*(dXidx + y*dWydx) )*epsxy *xz_on;
  *By = ( 0.875*Wy - 0.125*(dXidy + y*dWydy) )*epsxy;
  *Bz =            - 0.125*(dXidz + y*dWydz) *xz_on;
}

/* initialize BNSdata */
int BNSdata_startup(tGrid *grid)
{
  int corot1 = Getv("BNSdata_rotationstate1","corotation");
  int corot2 = Getv("BNSdata_rotationstate2","corotation");
  int b;
  int TOVav        = Getv("BNSdata_guess", "TOVaverage");
  int TOVprod      = Getv("BNSdata_guess", "TOVproduct");
  int initShift    = Getv("BNSdata_guess", "TaniguchiShift");
  int initFromChkp = Getv("BNSdata_guess", "initialize_from_checkpoint");
  double Omega     = Getd("BNSdata_Omega");
  double omegax1   = Getd("BNSdata_omegax1");
  double omegay1   = Getd("BNSdata_omegay1");
  double omegaz1   = Getd("BNSdata_omegaz1");
  double omegax2   = Getd("BNSdata_omegax2");
  double omegay2   = Getd("BNSdata_omegay2");
  double omegaz2   = Getd("BNSdata_omegaz2");
  double xCM       = Getd("BNSdata_x_CM");
  double kappa     = Getd("BNSdata_kappa");
  double BNSdata_b = Getd("BNSdata_b");
  double BNSdata_n = Getd("BNSdata_n");
  double Gamma     = 1.0 + 1.0/BNSdata_n;
  double rs1, m1, Phic1, Psic1, m01;
  double rs2, m2, Phic2, Psic2, m02;
  double xout1 = grid->box[0]->x_of_X[1](
                     (void *) grid->box[0], 0, 0.0,0.0,0.0);
  double xin1  = grid->box[0]->x_of_X[1](
                     (void *) grid->box[0], 0, 0.0,1.0,0.0);
  double xin2  = grid->box[3]->x_of_X[1](
                     (void *) grid->box[3], 0, 0.0,1.0,0.0);
  double xout2 = grid->box[3]->x_of_X[1](
                     (void *) grid->box[3], 0, 0.0,0.0,0.0);
  double xc1 = 0.5*(xout1+xin1);
  double xc2 = 0.5*(xin2+xout2);
  double ysh1;

  printf("Initializing BNSdata:\n");
/*
{
double A,B, x,y,z, Xp,Rp;
int dom;
for(A=0; A<=1; A++)
for(B=0; B<=1; B++)
for(dom=0; dom<=3; dom++)
{
xyz_of_AnsorgNS(grid->box[dom], -1, dom, A,B,1, &x,&y,&z, &Xp,&Rp);
printf("dom=%d (A,B)=(%g,%g): (x,y,z)=(%g,%.9e,%.9e)\n"
       "                      (Xp,Rp)=(%.9e,%.9e)\n\n", dom,A,B, x,y,z, Xp,Rp);
}
printf("================================\n");
printf("\n");
for(dom=0; dom<=3; dom++)
{
for(A=0; A<=1; A+=0.5)
for(B=0; B<=1; B+=0.5)
{
xyz_of_AnsorgNS(grid->box[dom], -1, dom, A,B,1, &x,&y,&z, &Xp,&Rp);
printf("dom=%d (A,B)=(%g,%g): (x,y,z)=(%g,%.9e,%.9e)\n"
       "                      (Xp,Rp)=(%.9e,%.9e)\n\n", dom,A,B, x,y,z, Xp,Rp);
}
printf("\n");
}
grid->time  = -777; 
write_grid(grid);
exit(77);
}
*/
  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("BNSdata_Psi",   1, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_Bx",    0, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_By",    0, 1, 1.0);
  VarNameSetBoundaryInfo("BNSdata_Bz",    0, 1, 1.0);
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
  enablevar(grid, Ind("BNSdata_SigmaX"));
  enablevar(grid, Ind("BNSdata_SigmaXX"));
  enablevar(grid, Ind("BNSdata_SigmaXXX"));
  enablevar(grid, Ind("BNSdata_lSigmaX"));
  enablevar(grid, Ind("BNSdata_lSigmaXX"));
  enablevar(grid, Ind("BNSdata_lSigmaXXX"));
  enablevar(grid, Ind("BNSdata_wBx"));
  enablevar(grid, Ind("BNSdata_q"));
  enablevar(grid, Ind("BNSdata_wBxx"));
  enablevar(grid, Ind("BNSdata_qx"));
  enablevar(grid, Ind("BNSdata_VRx"));
  enablevar(grid, Ind("BNSdata_temp1"));
  enablevar(grid, Ind("BNSdata_temp2"));
  enablevar(grid, Ind("BNSdata_temp3"));
  enablevar(grid, Ind("BNSdata_temp4"));
  enablevar(grid, Ind("BNSdata_qg"));
  enablevar(grid, Ind("BNSdata_Psiold"));
  enablevar(grid, Ind("BNSdata_Boldx"));
  enablevar(grid, Ind("BNSdata_alphaPold"));
  enablevar(grid, Ind("BNSdata_Sigmaold"));
  enablevar(grid, Ind("BNSdata_qgold"));
  enablevar(grid, Ind("BNSdata_qcorot"));
  enablevar(grid, Ind("BNSdata_qnocent"));
  enablevar(grid, Ind("BNSdata_surface_sigma_pm"));

  /* enable some lapse and shift of ADMvars */
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));
  
  /* set values of A,B,phi in box4/5 */
  set_BNSdata_ABphi(grid);

  /* set rs, m, Phic, Psic, m0 for both stars */
  TOV_init(P_core1, kappa, Gamma, 1, &rs1, &m1, &Phic1, &Psic1, &m01);
  TOV_init(P_core2, kappa, Gamma, 1, &rs2, &m2, &Phic2, &Psic2, &m02);

  /* set qmax1/2 */
  Setd("BNSdata_qmax1", pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
                        pow(P_core1, 1.0/(1.0 + BNSdata_n)));
  Setd("BNSdata_qmax2", pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
                        pow(P_core2, 1.0/(1.0 + BNSdata_n)));
  /* set cart positions of qmax1/2 */
  if(Getd("BNSdata_xmax1")<=0.0) Setd("BNSdata_xmax1", xc1);
  if(Getd("BNSdata_xmax2")>=0.0) Setd("BNSdata_xmax2", xc2);

  /* set cart positions of inner star edges */
  if(Getd("BNSdata_xin1")<=0.0) Setd("BNSdata_xin1", xin1);
  if(Getd("BNSdata_xin2")>=0.0) Setd("BNSdata_xin2", xin2);

  /* load data from some old checkpoint file */
  if(initFromChkp && GetsLax("outdir_previous_iteration")!=NULL)
  {
    char filename[10000];
    snprintf(filename, 9999, "%s/%s", 
             Gets("outdir_previous_iteration"), "checkpoint.0");
    prdivider(1);
    printf("loading initial guess from %s\n", filename);
    BNSgrid_load_initial_guess_from_checkpoint(grid, filename);
  }
  else if(initFromChkp && strlen(Gets("BNSdata_initfile"))>0)
  {
    char filename[10000];
    snprintf(filename, 9999, "%s", Gets("BNSdata_initfile"));
    prdivider(1);
    printf("loading initial guess from %s\n", filename);
    BNSgrid_load_initial_guess_from_checkpoint(grid, filename);
  }
  else /* use some TOV data */
  {
    /* get yshift1 (for testing) */
    if(m02==0.0) ysh1 = Getd("BNSdata_yshift1");
    else         ysh1 = 0.0;

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
      double *BNSdata_Psi    = box->v[Ind("BNSdata_Psi")];
      double *BNSdata_alphaP = box->v[Ind("BNSdata_alphaP")];
      double *BNSdata_q      = box->v[Ind("BNSdata_q")];
      double *BNSdata_Bx     = box->v[Ind("BNSdata_Bx")];
      double *BNSdata_By     = box->v[Ind("BNSdata_Bx")+1];
      double *BNSdata_Bz     = box->v[Ind("BNSdata_Bx")+2];
      double *BNSdata_Sigma  = box->v[Ind("BNSdata_Sigma")];
      double *BNSdata_wBx    = box->v[Ind("BNSdata_wBx")];
      double *BNSdata_wBy    = box->v[Ind("BNSdata_wBx")+1];
      double *BNSdata_wBz    = box->v[Ind("BNSdata_wBx")+2];
      double r1, m1_r, P1, Phi1, Psi1, m01_r, q1;
      double r2, m2_r, P2, Phi2, Psi2, m02_r, q2;
      double Bx1,By1,Bz1;
      double Bx2,By2,Bz2;

//double x,y,z;
//printf("m1=%g m2=%g Omega=%g  rs1=%g rs2=%g xc1=%g xc2=%g  xc1-xc2=%g\n",
//m1,m2, Omega, rs1,rs2, xc1,xc2, xc1-xc2);
//y=z=0;
//y=0.1; z=0.2;
//for(y=-5; y<5; y+=0.1)
//{
//z=0;
//x=0.1*y;
//BNSdata_initial_shift(1, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs1,
//                      -(x-xc1), -(y-ysh1), z,  &Bx1,&By1,&Bz1);
//BNSdata_initial_shift(2, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs2,
//                       x-xc2, y, z,  &Bx2,&By2,&Bz2);
//printf("%g  %g  %g\n", y, Bx1,Bx2);
//}
//exit(11);

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
        /* set Psi, alphaP, q */
        if(TOVav || TOVprod || (b==0 || b==1 || b==5) )
        {
          r1 = sqrt((x-xc1)*(x-xc1) + (y-ysh1)*(y-ysh1) + z*z);
          TOV_m_P_Phi_Psi_m0_OF_rf(r1, rs1, kappa, Gamma,
                                   P_core1, Phic1, Psic1,
                                   &m1_r, &P1, &Phi1, &Psi1, &m01_r);
          q1 = pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
               pow(P1, 1.0/(1.0 + BNSdata_n));
          BNSdata_initial_shift(1, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs1, 
                                -(x-xc1), -(y-ysh1), z,  &Bx1,&By1,&Bz1);
        }
        else { m1_r = P1 = Phi1 = m01_r = q1 = Bx1=By1=Bz1 = 0.0;  Psi1 = 1.0; }
        if(TOVav || TOVprod || (b==2 || b==3 || b==4) )
        {
          r2 = sqrt((x-xc2)*(x-xc2) + y*y + z*z);
          TOV_m_P_Phi_Psi_m0_OF_rf(r2, rs2, kappa, Gamma,
                                   P_core2, Phic2, Psic2,
                                   &m2_r, &P2, &Phi2, &Psi2, &m02_r);
          q2 = pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
               pow(P2, 1.0/(1.0 + BNSdata_n));
          BNSdata_initial_shift(2, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs2, 
                                x-xc2, y, z,  &Bx2,&By2,&Bz2);
        }
        else { m2_r = P2 = Phi2 = m02_r = q2 = Bx2=By2=Bz2 = 0.0;  Psi2 = 1.0; }
        
        /* set the data */
        if(TOVprod)
        {
          BNSdata_Psi[i]   = Psi1*Psi2;
          BNSdata_alphaP[i]= exp(Phi1+Phi2)*Psi1*Psi2;
          BNSdata_q[i]     = q1 + q2;
        }
        else /* for TOVav or TOV */
        {
          BNSdata_Psi[i]   = Psi1 + Psi2 - 1.0;
          BNSdata_alphaP[i]= exp(Phi1)*Psi1 + exp(Phi2)*Psi2 - 1.0;
          BNSdata_q[i]     = q1 + q2;
        }
        /* set inertial shift B^i if wanted */
        if(initShift)
        {
          BNSdata_Bx[i] = Bx1 + Bx2;
          BNSdata_By[i] = By1 + By2;
          BNSdata_Bz[i] = Bz1 + Bz2;
        }
        /* set Sigma and wB if needed */
        if( (b==0 || b==1 || b==5) && (!corot1) )
        {
          double Att;
          double A = pX[i];
          if(b==1) Att = 1.0-Attenuation01((A-0.1)/0.8, 2.0, 0.5);
          else     Att = 1.0;

          BNSdata_Sigma[i] = Omega*(xc1-xCM) * y * Att;
        }
        if( (b==3 || b==2 || b==4) && (!corot2) )
        {
          double Att;
          double A = pX[i];
          if(b==2) Att = 1.0-Attenuation01((A-0.1)/0.8, 2.0, 0.5);
          else     Att = 1.0;

          BNSdata_Sigma[i] = Omega*(xc2-xCM) * y * Att;
        }
      } /* end forallpoints */
    }

    /* set wB in both stars */
    BNS_set_wB(grid, 1, xc1,0.0,0.0);
    BNS_set_wB(grid, 2, xc2,0.0,0.0);
    
    /* if NS1 is shifted in y-direc. (for testing) adjust grid on right side */
    if(ysh1 != 0.0 && Getv("BNSdata_adjustdomain01","yes"))
    {
      /* adjust grid so that new q=0 is at A=0 */
      compute_new_q_and_adjust_domainshapes(grid, 0);

      /* set q to zero if q<0, and also in region 1 & 2 */
      set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);
    }
  } /* end intialization using TOV data */

  /* print out maxima */
  printf("BNSdata_startup: BNSdata_qmax1 = %g  BNSdata_qmax2 = %g\n"
         "                 BNSdata_xmax1 = %g  BNSdata_xmax2 = %g\n"
         "                  BNSdata_xin1 = %g   BNSdata_xin2 = %g\n",
         Getd("BNSdata_qmax1"), Getd("BNSdata_qmax2"),
         Getd("BNSdata_xmax1"), Getd("BNSdata_xmax2"),
         Getd("BNSdata_xin1"), Getd("BNSdata_xin2"));

  /* recompute q from the other fields */
  if(Getv("BNSdata_init_q_fromfields", "yes"))
  {
    /* adjust grid so that new q=0 is at A=0 */
    compute_new_q_and_adjust_domainshapes(grid, 0);
    compute_new_q_and_adjust_domainshapes(grid, 3);
  
    /* set q to zero if q<0, and also in region 1 & 2 */
    set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);
  }

  /* set qg=q, qgold=qg */
  varcopy(grid, Ind("BNSdata_qg"), Ind("BNSdata_q"));  
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));  

  /* save initial surface positions in BNSdata_surface_sigma_pm */
  save_surfacepos_in_BNSdata_surface_sigma_pm(1, grid, grid, 99999.0);
  save_surfacepos_in_BNSdata_surface_sigma_pm(2, grid, grid, 99999.0);

  /* set BNSdata_actual_xyzmax pars */
  set_BNSdata_actual_xyzmax_pars(grid);
  return 0;
}


/* functions to shift/translate a field/var */
void BNS_translate_var_in_box(tBox *box, int iv, int ivx,
                              double dx, double dy, double dz)
{
  int b = box->b;
  double *var = box->v[iv];
  double *dvdx = box->v[ivx];
  double *dvdy = box->v[ivx+1];
  double *dvdz = box->v[ivx+2];
  int i;

  /* get gradient of var */
  FirstDerivsOf_S(box, iv, ivx);
  
  /* v_translated(x) = v(x-dx) ~ v(x) - [vq(x)/dx] dx */
  //printf("b=%d: dx=%g dy=%g dz=%g\n", b, dx,dy,dz);
  forallpoints(box, i)
    var[i] = var[i] - (dvdx[i]*dx + dvdy[i]*dy + dvdz[i]*dz);
}
/* translate var on one side of grid specified by star */
void BNS_translate_var_aroundStar(tGrid *grid, int star, int iv, int ivx,
                                  double dx, double dy, double dz)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    /* continue with next box if on wrong side */
    if(star==1 && (b>=2 && b<=4)) continue;
    if(star==2 && (b<=1 || b==5)) continue;     
    BNS_translate_var_in_box(box, iv, ivx, dx,dy,dz);
  }
}   

/* center Psi, beta^i, alphaP, Sigma around each star */
int BNSdata_center_fields_if_desired(tGrid *grid, int it)
{
  if(    !Getv("BNSdata_center_fields", "no")    &&
     it>=Geti("BNSdata_center_fields_first_at")  && 
         Geti("BNSdata_center_fields_first_at")>=0 )
  {
    int b,i;
    double m01, m02;
    int iPsi = Ind("BNSdata_Psi");
    int iBx  = Ind("BNSdata_Bx");
    int iBy  = Ind("BNSdata_By");
    int iBz  = Ind("BNSdata_Bz");
    int ialP = Ind("BNSdata_alphaP");
    int iSig = Ind("BNSdata_Sigma");
    int itx  = Ind("BNSdata_temp1");
    double dx,dy,dz, x1,y1,z1, x2,y2,z2;
    double fac = Getd("BNSdata_center_fields_fac");
    double xmax1 = Getd("BNSdata_xmax1");
    double xmax2 = Getd("BNSdata_xmax2");
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("BNSdata_center_fields", "center_yz"))  xon=0;

    /* get global max of q in NS1/2 */
    x1 = Getd("BNSdata_actual_xmax1");
    y1 = Getd("BNSdata_actual_ymax1"); 
    z1 = Getd("BNSdata_actual_zmax1");
    x2 = Getd("BNSdata_actual_xmax2");
    y2 = Getd("BNSdata_actual_ymax2");
    z2 = Getd("BNSdata_actual_zmax2");

    printf("Centering fields:  BNSdata_center_fields_fac = %s\n"
           " BNSdata_center_fields = %s\n",
           Gets("BNSdata_center_fields_fac"),
           Gets("BNSdata_center_fields"));

    /* get dx,dy,dz and center around star1 */
    dx = -fac*(x1-xmax1);
    dy = -fac*y1;
    dz = -fac*z1;
    dx = dx*xon;  /* <-- switch centering on/off */
    dy = dy*yon;
    dz = dz*zon;
    BNS_translate_var_aroundStar(grid,1, iPsi,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,1, iBx ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,1, iBy ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,1, iBz ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,1, ialP,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,1, iSig,itx, dx,dy,dz);

    /* get dx,dy,dz and center around star2 */
    dx = -fac*(x2-xmax2);
    dy = -fac*y2;
    dz = -fac*z2;
    dx = dx*xon;  /* <-- switch centering on/off */
    dy = dy*yon;
    dz = dz*zon;
    BNS_translate_var_aroundStar(grid,2, iPsi,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,2, iBx ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,2, iBy ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,2, iBz ,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,2, ialP,itx, dx,dy,dz);
    BNS_translate_var_aroundStar(grid,2, iSig,itx, dx,dy,dz);

    /* recalc q as well? */
    if(Getv("BNSdata_center_fields", "reset_q"))
    {
      if(Getv("BNSdata_center_fields", "adjust_domainshapes"))
      {
        compute_new_q_and_adjust_domainshapes(grid, 0);
        compute_new_q_and_adjust_domainshapes(grid, 3);
      }
      else
        BNS_compute_new_centered_q(grid);

      /* set q to zero if q<0 or in region 1 and 2 */
      set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);
    }
   
    /* enforce uniqueness on axis:
       set vars equal to val at phi=0 for all phi>0 */
    if(Getv("BNSdata_uniqueness_on_axis", "yes"))
      BNS_enforce_uniqueness_on_axis(vlu);

    /* set actual positions of maxima again */
    set_BNSdata_actual_xyzmax_pars(grid);

    /* print new masses */
    m01 = GetInnerRestMass(grid, 0);
    m02 = GetInnerRestMass(grid, 3);
    printf("     => m01=%.19g m02=%.19g\n", m01, m02);
  }
  return 0;
}

/* keep inner edges of stars in place */
int BNSdata_keep_xin_if_desired(tGrid *grid, int it)
{
  if(    !Getv("BNSdata_keep_xin", "no")    &&
     it>=Geti("BNSdata_keep_xin_first_at")  && 
         Geti("BNSdata_keep_xin_first_at")>=0 )
  {
    int b,i;
    double m01, m02;
    int iPsi = Ind("BNSdata_Psi");
    int iBx  = Ind("BNSdata_Bx");
    int iBy  = Ind("BNSdata_By");
    int iBz  = Ind("BNSdata_Bz");
    int ialP = Ind("BNSdata_alphaP");
    int iSig = Ind("BNSdata_Sigma");
    int itx  = Ind("BNSdata_temp1");
    double dx, xin1, xin2;
    double fac = Getd("BNSdata_keep_xin_fac");
    double BNSdata_xin1 = Getd("BNSdata_xin1");
    double BNSdata_xin2 = Getd("BNSdata_xin2");

    /* get pos of actual xin1/2 */
    xin1 = grid->box[0]->x_of_X[1]((void *) grid->box[0], -1, 0.0,1.0,0.0);
    xin2 = grid->box[3]->x_of_X[1]((void *) grid->box[3], -1, 0.0,1.0,0.0);
    printf("Actual:          xin1=%.12g         xin2=%.12g\n", xin1, xin2);
    printf("Desired: BNSdata_xin1=%.12g BNSdata_xin2=%.12g\n",
           BNSdata_xin1, BNSdata_xin2);
    printf("Shifting fields:  BNSdata_keep_xin_fac = %s\n"
           " BNSdata_keep_xin = %s\n",
           Gets("BNSdata_keep_xin_fac"),
           Gets("BNSdata_keep_xin"));

    /* get dx and shift around star1 */
    dx = -fac*(xin1-BNSdata_xin1);
    BNS_translate_var_aroundStar(grid,1, iPsi,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,1, iBx ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,1, iBy ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,1, iBz ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,1, ialP,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,1, iSig,itx, dx,0.0,0.0);

    /* get dx and shift around star2 */       
    dx = -fac*(xin2-BNSdata_xin2);
    BNS_translate_var_aroundStar(grid,2, iPsi,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,2, iBx ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,2, iBy ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,2, iBz ,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,2, ialP,itx, dx,0.0,0.0);
    BNS_translate_var_aroundStar(grid,2, iSig,itx, dx,0.0,0.0);

    /* recalc q as well? */
    if(Getv("BNSdata_keep_xin", "reset_q"))
    {
      if(Getv("BNSdata_keep_xin", "adjust_domainshapes"))
      {
        compute_new_q_and_adjust_domainshapes(grid, 0);
        compute_new_q_and_adjust_domainshapes(grid, 3);
      }
      else
        BNS_compute_new_centered_q(grid);

      /* set q to zero if q<0 or in region 1 and 2 */
      set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);
    }
   
    /* enforce uniqueness on axis:
       set vars equal to val at phi=0 for all phi>0 */
    if(Getv("BNSdata_uniqueness_on_axis", "yes"))
      BNS_enforce_uniqueness_on_axis(vlu);

    /* set actual positions of maxima again */
    set_BNSdata_actual_xyzmax_pars(grid);

    /* print new masses */
    m01 = GetInnerRestMass(grid, 0);
    m02 = GetInnerRestMass(grid, 3);
    printf("     => m01=%.19g m02=%.19g\n", m01, m02);
  }
  return 0;
}

/* functions to compute new q from fields, and to also shift them so
   that the max in q is centered where we want it on the x-axis */
void BNS_compute_new_centered_q(tGrid *grid)
{
  int iq = Ind("BNSdata_q");
  int iqg= Ind("BNSdata_qg");

  BNS_compute_new_q(grid, iq);
  if(Getv("BNSdata_center_new_q_flag", "yes"))
  {
    int iqx= Ind("BNSdata_qx");
    int b, i;
    int bi1, bi2;
    double dx,dy,dz, x1,y1,z1, x2,y2,z2;
    double fac = Getd("BNSdata_center_new_q_fac");
    double xmax1 = Getd("BNSdata_xmax1");
    double xmax2 = Getd("BNSdata_xmax2");
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("BNSdata_center_new_q", "center_yz"))  xon=0;

    /* get global max of q in NS1/2 */
    x1 = Getd("BNSdata_actual_xmax1");
    y1 = Getd("BNSdata_actual_ymax1"); 
    z1 = Getd("BNSdata_actual_zmax1");
    x2 = Getd("BNSdata_actual_xmax2");
    y2 = Getd("BNSdata_actual_ymax2");
    z2 = Getd("BNSdata_actual_zmax2");
           
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *BNSdata_q = box->v[iq];
      double *dqdx = box->v[iqx];
      double *dqdy = box->v[iqx+1];
      double *dqdz = box->v[iqx+2];
      
      /* get gradient of q */
      FirstDerivsOf_S(box, iq, iqx);
      
      /* q_centered(x) = q(x+dx) ~ q(x) + [dq(x)/dx] dx */
      if(b<=1 || b==5) { dx = x1-xmax1; dy = y1; dz = z1; }
      else             { dx = x2-xmax2; dy = y2; dz = z2; }
      dx = dx*xon;  /* <-- switch centering on/off */
      dy = dy*yon;
      dz = dz*zon;
      //printf("b=%d: dx=%g dy=%g dz=%g\n", b, dx,dy,dz);
      forallpoints(box, i)
        BNSdata_q[i] = BNSdata_q[i] + 
                       fac*(dqdx[i]*dx + dqdy[i]*dy + dqdz[i]*dz);
    }
  }
  varcopy(grid, iqg, iq); /* set qg=q */
}
double BNS_compute_new_centered_q_atXYZ(tGrid *grid, int bi,
                                        double X, double Y, double Z)
{
  double q;
  q = BNS_compute_new_q_atXYZ(grid,bi, X,Y,Z);
  if(Getv("BNSdata_center_new_q_flag", "yes"))
  {
    int iq = Ind("BNSdata_temp4");
    int iqx= Ind("BNSdata_qx");
    tBox *box = grid->box[bi];
    double *cx= box->v[Ind("BNSdata_temp1")];
    double *cy= box->v[Ind("BNSdata_temp2")];
    double *cz= box->v[Ind("BNSdata_temp3")];
    double dx,dy,dz, xm,ym,zm;
    double dqdx, dqdy, dqdz;
    double fac = Getd("BNSdata_center_new_q_fac");
    double xmax;
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("BNSdata_center_new_q", "center_yz"))  xon=0;

    if(bi<=1 || bi==5)
    {
      /* desired xmax1 */
      xmax = Getd("BNSdata_xmax1");
      /* get global max of q in NS1/2 */
      xm = Getd("BNSdata_actual_xmax1");
      ym = Getd("BNSdata_actual_ymax1"); 
      zm = Getd("BNSdata_actual_zmax1");
    }
    else
    {
      /* desired xmax1 */
      xmax = Getd("BNSdata_xmax2");
      /* get global max of q in NS1/2 */
      xm = Getd("BNSdata_actual_xmax2");
      ym = Getd("BNSdata_actual_ymax2");
      zm = Getd("BNSdata_actual_zmax2");
    }

    /* get gradient of q and its coeffs, set coeffs of dq in 
       BNSdata_temp1/2/3 */
    BNS_compute_new_q(grid, iq);
    FirstDerivsOf_S(box, iq, iqx);
    spec_Coeffs(box, box->v[iqx], cx);
    spec_Coeffs(box, box->v[iqx+1], cy);
    spec_Coeffs(box, box->v[iqx+2], cz);

    /* find dqdx, dqdy, dqdz by interpolation */
    dqdx = spec_interpolate(box, cx, X,Y,Z);
    dqdy = spec_interpolate(box, cy, X,Y,Z);
    dqdz = spec_interpolate(box, cz, X,Y,Z);

    /* q_centered(x) = q(x+dx) ~ q(x) + [dq(x)/dx] dx */
    dx = xm-xmax;
    dy = ym;
    dz = zm;
    dx = dx*xon;  /* <-- switch centering on/off */
    dy = dy*yon;
    dz = dz*zon;
//printf("bi=%d: dx=%g dy=%g dz=%g\n", bi, dx,dy,dz);
    q = q + fac*(dqdx*dx + dqdy*dy + dqdz*dz);
  }
  return q;
}

/* center q if needed, depending on how pars are set */
int BNSdata_center_q_if_desired(tGrid *grid, int it)
{
  if(    !Getv("BNSdata_center_new_q", "no")    &&
     it>=Geti("BNSdata_center_new_q_first_at")  && 
         Geti("BNSdata_center_new_q_first_at")>=0 )
  {
    int b,i;
    double m01, m02;

    printf("Centering q:  BNSdata_center_new_q_fac = %s\n"
           " BNSdata_center_new_q = %s\n",
           Gets("BNSdata_center_new_q_fac"), Gets("BNSdata_center_new_q"));
    Sets("BNSdata_center_new_q_flag", "yes");  /* activate centering of q */
    if(Getv("BNSdata_center_new_q", "adjust_domainshapes"))
    {
      compute_new_q_and_adjust_domainshapes(grid, 0);
      compute_new_q_and_adjust_domainshapes(grid, 3);
    }
    else 
    {
      BNS_compute_new_centered_q(grid);
    }
    Sets("BNSdata_center_new_q_flag", "no");  /* deactivate centering of q */
    /* set q to zero if q<0 or in region 1 and 2 */
    set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

    /* enforce uniqueness on axis:
       set vars equal to val at phi=0 for all phi>0 */
    if(Getv("BNSdata_uniqueness_on_axis", "yes"))
      BNS_enforce_uniqueness_on_axis(vlu);

    /* set actual positions of maxima again */
    set_BNSdata_actual_xyzmax_pars(grid);

    /* print new masses */
    m01 = GetInnerRestMass(grid, 0);
    m02 = GetInnerRestMass(grid, 3);
    printf("     => m01=%.19g m02=%.19g\n", m01, m02);
  }
  return 0;
}


/* find both qmax and reset BNSdata_qmax1/2, BNSdata_xmax1/2 accordingly */
void find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars(tGrid *grid)
{
  int bi1, bi2;
  double xmax1, qmax1, xmax2, qmax2;
  double Xmax1,Ymax1, Xmax2,Ymax2;

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,0.0);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,0.0);
  else
    xmax2 = Xmax2;
  /* set qmax1/2 */
  Setd("BNSdata_qmax1", qmax1);
  Setd("BNSdata_qmax2", qmax2);
  /* set cart positions of qmax1/2 */
  Setd("BNSdata_xmax1", xmax1);
  Setd("BNSdata_xmax2", xmax2);
  printf("find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars:\n");
  printf("  resetting: BNSdata_qmax1 = %g  BNSdata_qmax2=%g\n"
         "             BNSdata_xmax1 = %g  BNSdata_xmax2=%g\n",
         Getd("BNSdata_qmax1"), Getd("BNSdata_qmax2"),
         Getd("BNSdata_xmax1"), Getd("BNSdata_xmax2"));
}


/* adjust C1/2 and Omega as in Pedro's NS ini dat paper.
   BUT with this algorithm the iterations fail... */
int adjust_C1_C2_Omega_q_Pedro(tGrid *grid, int it, double tol)
{
  double Cvec[3];
  double m0errorvec[3];
  double m01, m02, m0_error, dm01, dm02;
  double Omega, dOmega=1e-3; 
  double L2norm1,L2norm2,L2norm3;
  int check, stat, bi, i;
  t_grid_grid0_m01_m02_struct pars[1];

  /* rest masses before adjusting q */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("BNSdata_solve: rest mass in inner domains before computing new q:\n"
         "   m01=%g  m02=%g\n", m01, m02);

//    /* set new q on grid */
//    BNS_compute_new_centered_q(grid);
//
//    /* compute masses after computing new q */
//    m01 = GetInnerRestMass(grid, 0);
//    m02 = GetInnerRestMass(grid, 3);
//    printf("rest mass in inner domains after computing new q: "
//           "m01=%g m02=%g\n", m01, m02);
//
//    /* compute masses after adjusting grid */
//    compute_new_q_and_adjust_domainshapes(grid, 0);
//    compute_new_q_and_adjust_domainshapes(grid, 3);
//    m01 = GetInnerRestMass(grid, 0);
//    m02 = GetInnerRestMass(grid, 3);
//    printf("rest mass after adjusting grid to new q: "
//           "m01=%g m02=%g\n", m01, m02);
//
  /* compute error in masses */
  dm01 = m01 - Getd("BNSdata_m01");
  dm02 = m02 - Getd("BNSdata_m02");
  m0_error = dm01*dm01 + dm02*dm02;
  m0_error = sqrt(m0_error)/(Getd("BNSdata_m01")+Getd("BNSdata_m02"));
  printf("BNSdata_solve step %d: rest mass error = %.4e "
         "(before adjusting C1/2)\n", it, m0_error);

  /* save Omega */
  Omega = Getd("BNSdata_Omega");
  //dOmega= Omega * m0_error;
  printf("BNSdata_solve step %d: old Omega = %.4e  dOmega = %.4e\n",
         it, Getd("BNSdata_Omega"), dOmega);
         
  /* save old qg in BNSdata_qgold */
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* compute L2-diff between new q and qold for Omega - dOmega */
  Setd("BNSdata_Omega", Omega - dOmega);
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  L2norm3 = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega + dOmega */
  Setd("BNSdata_Omega", Omega + dOmega);
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  L2norm2 = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega */
  Setd("BNSdata_Omega", Omega);
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  L2norm1 = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1"));

  printf("BNSdata_solve step %d: L2norm1=%g\n", it, L2norm1);
  printf("BNSdata_solve step %d: L2norm2=%g\n", it, L2norm2);
  printf("BNSdata_solve step %d: L2norm3=%g\n", it, L2norm3);
  if(L2norm1<=L2norm2 && L2norm1<=L2norm3)
  { Setd("BNSdata_Omega", Omega);  dOmega=dOmega*0.5; }
  if(L2norm2<L2norm1  && L2norm2<=L2norm3) 
    Setd("BNSdata_Omega", Omega+dOmega);
  if(L2norm3<L2norm1  && L2norm3<L2norm2)
    Setd("BNSdata_Omega", Omega-dOmega);
  printf("BNSdata_solve step %d: new Omega = %.4e  dOmega = %.4e\n",
         it, Getd("BNSdata_Omega"), dOmega);
  BNS_compute_new_centered_q(grid);
  
/*
{
int b;
  forallboxes(grid, b)
  {
    double *t = grid->box[b]->v[Ind("BNSdata_temp1")];;
    double *q = grid->box[b]->v[Ind("BNSdata_q")];;
    double *qo = grid->box[b]->v[Ind("BNSdata_qold")];;
    forallpoints(grid->box[b], i)
    { t[i] = -(b+100);  q[i] = -(b+20);  qo[i] = -(b+30);}
  }
}
*/

  /* set desired masses for this iteration */
  pars->m01 = Getd("BNSdata_m01"); // + dm01*0.9;
  pars->m02 = Getd("BNSdata_m02"); // + dm02*0.9;
  printf("BNSdata_solve step %d: "
         "adjusting q,C1,C2 to achieve: m01=%g  m02=%g\n",
         it, pars->m01, pars->m02);

  /* print C1/2 we used before */
  printf("old: BNSdata_C1=%g BNSdata_C2=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_C2"));
  printf("     => m01=%.19g m02=%.19g\n", m01, m02);

  /* choose C1/2 such that rest masses are not too big or too small */
  for(i=0; i<1000; i++)
  {
    double *q_b1 = grid->box[1]->v[Ind("BNSdata_q")];
    double *q_b2 = grid->box[2]->v[Ind("BNSdata_q")];

    BNS_compute_new_centered_q(grid);
    m01 = GetInnerRestMass(grid, 0);
    m02 = GetInnerRestMass(grid, 3);

    check = 0;

    if(m01 > 1.1*Getd("BNSdata_m01"))
      { Setd("BNSdata_C1", 0.999*Getd("BNSdata_C1"));  check=1; }
    else if(m01 < 0.9*Getd("BNSdata_m01"))
      { Setd("BNSdata_C1", 1.002*Getd("BNSdata_C1"));  check=1; }

    if(m02 > 1.1*Getd("BNSdata_m02") && Getd("BNSdata_m02")>0)
      { Setd("BNSdata_C2", 0.999*Getd("BNSdata_C2"));  check=1; }
    else if(m02 < 0.9*Getd("BNSdata_m02") && Getd("BNSdata_m02")>0)
      { Setd("BNSdata_C2", 1.002*Getd("BNSdata_C2"));  check=1; }

    if(check==0) break;
  }

  /* refine guess for C1/2 */
  pars->grid = grid;
  Cvec[1] = Getd("BNSdata_C1");
  stat = newton_linesrch_itsP(Cvec, 1, &check, m01_guesserror_VectorFuncP,
                              (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  Setd("BNSdata_C1", Cvec[1]);

  Cvec[1] = Getd("BNSdata_C2");
  if(Getd("BNSdata_m02")>0)
    stat = newton_linesrch_itsP(Cvec, 1, &check, m02_guesserror_VectorFuncP,
                                (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  Setd("BNSdata_C2", Cvec[1]);

  /* print guess for C1/2 */                                        
  printf("guess: BNSdata_C1=%g BNSdata_C2=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_C2"));
//Yo(3);
//CheckIfFinite(grid,  "BNSdata_q");

  /***********************************************************************/
  /* do newton_linesrch_itsP iterations of Cvec until m0errorvec is zero */
  /***********************************************************************/
  pars->grid = grid;
  pars->grid0= grid;
  /* adjust C1 and thus m01 */
  Cvec[1] = Getd("BNSdata_C1");
  stat = newton_linesrch_itsP(Cvec, 1, &check, m01_error_VectorFuncP,
                              (void *) pars, 1000, tol*0.01);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
  Setd("BNSdata_C1", Cvec[1]);

  /* adjust C2 and thus m02 */
  Cvec[1] = Getd("BNSdata_C2");
  if(Getd("BNSdata_m02")>0)
    stat = newton_linesrch_itsP(Cvec, 1, &check, m02_error_VectorFuncP,
                                (void *) pars, 1000, tol*0.01);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
  Setd("BNSdata_C2", Cvec[1]);
  printf("new: BNSdata_C1=%g BNSdata_C2=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_C2"));

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* print new masses */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("     => m01=%.19g m02=%.19g\n", m01, m02);

  return 0;
}

/* adjust C1/2 and Omega by keeping central q constant.
   Then recompute q [as in Bonazzolla, Gourgoulkon, Marck (BGM)]. */
int adjust_C1_C2_Omega_q_BGM(tGrid *grid, int it, double tol)
{
  double C1OC2xCMvec[5];  /* contains C1, Omega, C2, xCM */
  double m01, m02;
  int check, stat, bi, i;

//Setd("BNSdata_Omega", 0.01);

  /* save old q in BNSdata_qold */
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* rest masses before adjusting q */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("BNSdata_solve: rest mass in inner domains before computing new q:\n"
         "   m01=%g  m02=%g\n", m01, m02);

  /* print C1/2 we used before */
  printf("old: BNSdata_C1=%g BNSdata_Omega=%g BNSdata_C2=%g BNSdata_x_CM=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_Omega"), 
         Getd("BNSdata_C2"), Getd("BNSdata_x_CM"));


  /**********************************************************************/
  /* do newton_linesrch_its iterations of Cvec until m0errorvec is zero */
  /**********************************************************************/
  central_q_errors_VectorFunc__grid = grid;

  /* guess C1, C2 */
  C1OC2xCMvec[1] = Getd("BNSdata_C1");
  C1OC2xCMvec[2] = Getd("BNSdata_C2");
  stat = newton_linesrch_its(C1OC2xCMvec, 2, &check,
                             estimate_q_errors_VectorFunc, 1000, tol*0.01);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  Setd("BNSdata_C1", C1OC2xCMvec[1]);
  Setd("BNSdata_C2", C1OC2xCMvec[2]);

//  /* adjust C1, C2, Omega and xCM */
//  C1OC2xCMvec[1] = Getd("BNSdata_C1");
//  C1OC2xCMvec[2] = Getd("BNSdata_Omega");
//  C1OC2xCMvec[3] = Getd("BNSdata_C2");
//  C1OC2xCMvec[4] = Getd("BNSdata_x_CM");
//  stat = newton_linesrch_its(C1OC2xCMvec, 3, &check,
//                             central_q_errors_VectorFunc, 1000, tol*0.01);
//  // or maybe:
//  //stat = newton_linesrch_its(C1OC2xCMvec, 4, &check,
//  //                           central_q_errors_VectorFunc, 1000, tol*0.01);
//  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
//  Setd("BNSdata_C1", C1OC2xCMvec[1]);
//  Setd("BNSdata_Omega", C1OC2xCMvec[2]);
//  Setd("BNSdata_C2", C1OC2xCMvec[3]);
//  Setd("BNSdata_x_CM", C1OC2xCMvec[4]);

//printf("HACK: skip adjust C2, and set C2=C1\n");
//Setd("BNSdata_C2", C1OC2xCMvec[1]);

  printf("new: BNSdata_C1=%g BNSdata_Omega=%g BNSdata_C2=%g BNSdata_x_CM=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_Omega"), 
         Getd("BNSdata_C2"), Getd("BNSdata_x_CM"));

  /* adjust grid so that new q=0 is at A=0 */
  compute_new_q_and_adjust_domainshapes(grid, 3);
  compute_new_q_and_adjust_domainshapes(grid, 0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* print new masses */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("     => m01=%.19g m02=%.19g\n", m01, m02);

  return 0;
}


/* set the pars BNSdata_desired_VolAvSigma12 to values that min BC violation */
void set_BNSdata_desired_VolAvSigma12_toMinBCerr(tGrid *grid, int index_Sigma)
{
  int ijk;
  int AddInnerVolIntToBC=Getv("BNSdata_Sigma_surface_BCs","AddInnerVolIntToBC");
  int InnerVolIntZero = Getv("BNSdata_Sigma_surface_BCs", "InnerVolIntZero");
  double *Sigma;
  double VolAvSigma1, VolAvSigma2;

  /* do nothing? */
  if(Getv("BNSdata_set_desired_VolAvSigmas", "no")) return;

  /* set VolAvSigma1/2 to sum of Sigmas */
  VolAvSigma1 = 0.0;
  Sigma      = grid->box[0]->v[index_Sigma];
  forallpoints(grid->box[0], ijk) VolAvSigma1 += Sigma[ijk]; 
  VolAvSigma2 = 0.0;
  Sigma      = grid->box[3]->v[index_Sigma];
  forallpoints(grid->box[3], ijk) VolAvSigma2 += Sigma[ijk]; 

  /* use VolInt instead in some other cases */
  if (AddInnerVolIntToBC || InnerVolIntZero)
  {
    VolAvSigma1 = VolumeIntegral_inBNSgridBox(grid, 0, index_Sigma); 
    VolAvSigma2 = VolumeIntegral_inBNSgridBox(grid, 3, index_Sigma); 
  }
  Setd("BNSdata_desired_VolAvSigma1", VolAvSigma1);
  Setd("BNSdata_desired_VolAvSigma2", VolAvSigma2);
  printf(" setting: BNSdata_desired_VolAvSigma1 / 2 = %g / %g\n",
  VolAvSigma1, VolAvSigma2);
}


/* backup grid,pdb to grid_bak,pdb_bak.
   But do it only if BNSdata_domainshape_diff_tol<1e30. 
   Call as:    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak); */
void backup_grid_pdb(tGrid *grid, tParameter *pdb,
                     tGrid *grid_bak, tParameter *pdb_bak)
{
  /* do nothing if we tolerate large differences */ 
  if(Getd("BNSdata_domainshape_diff_tol")>=1e30) return;

  /* make exact copies of grid and pdb */
  copy_grid(grid, grid_bak, 0);
  copy_pdb(pdb, npdb, pdb_bak);
}

/* save star surfaces (given by sigma_pm) in BNSdata_surface_sigma_pm on
   both grids. */
void save_surfacepos_in_BNSdata_surface_sigma_pm(int star, tGrid *grid,
                                                 tGrid *grid_bak, double t)
{
  static double latest_saved_t1 = 1e30;
  static double latest_saved_t2 = 1e30;
  double latest_saved_t;
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int surfi  = Ind("BNSdata_surface_sigma_pm");
  int b, i,j,k, n1,n2,n3, ijk,ijk1;
  double *surf_sig;
  double *surf;
  double *surf_bak;

  /* save surface for one star */
  if(star==1)      { b=0; latest_saved_t=latest_saved_t1; }
  else if(star==2) { b=3; latest_saved_t=latest_saved_t2; }
  else 
    errorexit("save_surfacepos_in_BNSdata_surface_sigma_pm: star must be 1 or 2.");
  n1=(grid)->box[b]->n1;
  n2=(grid)->box[b]->n2;
  n3=(grid)->box[b]->n3;
  surf_sig = (grid)->box[b]->v[sigpmi];
  surf     = (grid)->box[b]->v[surfi];
  surf_bak = (grid_bak)->box[b]->v[surfi];

  /* save previous surfaces only if t!=latest_saved_t */
  if(!dequal(t,latest_saved_t))
  {
    /* save previous surfaces on both grids */
    for (k = 0; k < n3; k++)
    for (j = 0; j < n2; j++)
    for (i = n1-1; i > 0; i--)
    {
      ijk = Index(i,j,k);
      ijk1= Index(i-1,j,k);
      surf[ijk]    = surf[ijk1];
      if(surf_bak) surf_bak[ijk]= surf_bak[ijk1];
    }
    /* store t at which we saved surfaces for the last time */
    if(star==1)      latest_saved_t1 = t;
    else if(star==2) latest_saved_t2 = t;
  }

  /* set current surface on both grids */
  forplane1(i,j,k, n1,n2,n3, 0)
  {
    ijk = Index(i,j,k);
    surf[ijk] = surf_sig[ijk];
    if(surf_bak) surf_bak[ijk] = surf_sig[ijk];
  }
}

/* compare grid,pdb to grid_bak,pdb_bak for one star and keep grid,pdb
   if the surface difference is small, otherwise restore grid,pdb
   to grid_bak,pdb_bak
   Call as:
    restore_grid_pdb_if_change_in_star_is_large(1/2, grid,pdb, grid_bak,pdb_bak);
*/
void restore_grid_pdb_if_change_in_star_is_large(int star,
                                      tGrid *grid, tParameter *pdb,
                                      tGrid *grid_bak, tParameter *pdb_bak)
{
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int surfi  = Ind("BNSdata_surface_sigma_pm");
  int b, i,j,k, n1,n2,n3, ijk,ijk1, n;
  double *surf_sig_bak;
  double *surf;
  double diff_new_old, diff_new_sig;

  /* save current and previous surfaces in BNSdata_surface_sigma_pm  */
  save_surfacepos_in_BNSdata_surface_sigma_pm(star, grid, grid_bak, grid->time);

  /* do nothing if we tolerate large differences */  
  if(Getd("BNSdata_domainshape_diff_tol")>=1e30) return;

  /* compute diff between grid_bak and grid */
  if(star==1) b=0;
  else if(star==2) b=3;
  else 
    errorexit("restore_grid_pdb_if_change_in_star_is_large: star must be 1 or 2.");
  n1=(grid)->box[b]->n1;
  n2=(grid)->box[b]->n2;
  n3=(grid)->box[b]->n3;
  surf_sig_bak = grid_bak->box[b]->v[sigpmi];
  surf         = grid->box[b]->v[surfi];

  /* L2-norm diffs: surf^n - surf_sig_bak, surf^n - surf^{n-1} */
  diff_new_old = diff_new_sig = 0.0; n=0;
  forplane1(i,j,k, n1,n2,n3, 0)
  {
    ijk = Index(i,j,k);
    diff_new_sig += (surf[ijk]-surf_sig_bak[ijk])*(surf[ijk]-surf_sig_bak[ijk]);
    ijk1 = Index(i+1,j,k); /* index of previous surface */
    diff_new_old += (surf[ijk]-surf[ijk1])*(surf[ijk]-surf[ijk1]);
    n++;
  }
  diff_new_sig = sqrt(diff_new_sig/n);
  diff_new_old = sqrt(diff_new_old/n);
  printf("restore_grid_pdb_if_change_in_star_is_large:\n"
         " diff_new_old=%g diff_new_sig=%g\n", diff_new_old, diff_new_sig);

  /* Check if diff_new_old is small. 
     Compare arXiv:0804.3787, III C. Algorithm point 3. */
  /* if diff_new_old is small, keep new grid */
  if(diff_new_old < Getd("BNSdata_domainshape_diff_tol")*diff_new_sig)
  { 
    printf(" adjusted domain shape around star%d.\n", star);
    return;
  }
  else /* otherwise restore old grid */
  {
    copy_grid(grid_bak, grid, 0);
    copy_pdb(pdb_bak, npdb, pdb);
    printf(" kept domain shape around star%d:", star);
    printf(" BNSdata_C1=%g BNSdata_C2=%g\n",
           Getd("BNSdata_C1"), Getd("BNSdata_C2"));
  }
}

/* average surface shapes for one star */
void average_current_and_old_surfaceshape(double weight, int innerdom,
                                          tGrid *grid, tParameter *pdb,
                                          tGrid *grid_bak, tParameter *pdb_bak)
{
  tGrid *grid2;
  int interp_qgold = !Getv("BNSdata_new_q", "FromFields");
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int b, i, outerdom;
  void (*Interp_From_Grid1_To_Grid2)(tGrid *grid1, tGrid *grid2, int vind,
                                     int innerdom);
  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("average_current_and_old_surfaceshape: "
                 "innerdom is not 0 or 3");

  /* do nothing if weight=1 */  
  if(weight==1.0) return;

  printf("average_current_and_old_surfaceshape:  weight=%f  innerdom=%d\n", 
         weight, innerdom);

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* loop over AnsorgNS boxes */
  for(b=0; b<=3; b++)
  {
    double *surf_bak;
    double *surf2;

    /* do nothing if we are in the wrong box */
    if(b!=innerdom && b!=outerdom) continue;

    surf_bak = grid_bak->box[b]->v[sigpmi];
    surf2    = grid2->box[b]->v[sigpmi];

    /* set val on grid2 to weighted average of grid and grid_bak */      
    forallpoints(grid->box[b], i)
      surf2[i] = weight*surf2[i] + (1.0-weight)*surf_bak[i];
  }

  /* initialize coords of grid2 on side of innerdom */
  BNSgrid_init_Coords_pm(grid2, innerdom);
  /* use interpolator that does only side of innerdom */
  Interp_From_Grid1_To_Grid2 = Interp_Var_From_Grid1_To_Grid2_pm;

  /* interpolate some vars from grid onto new grid2 */
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"),innerdom);
  if( (innerdom==0 && !Getv("BNSdata_rotationstate1","corotation")) ||
      (innerdom==3 && !Getv("BNSdata_rotationstate2","corotation"))   )
  {
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"),innerdom);
  }
  //  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"),innerdom);
  if(interp_qgold) 
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);

  /* q on grid2 can just be computed */
  BNS_compute_new_centered_q(grid2);
//  /* set q to zero if q<0 or in region 1 and 2 */
//  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);
}

/* add a const to sigma_pm to change surface shape for one star */
void add_const_to_sigma_pm(tGrid *grid, double dsig, int innerdom)
{
  tGrid *grid2;
  int interp_qgold = !Getv("BNSdata_new_q", "FromFields");
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int b, i, outerdom;
  void (*Interp_From_Grid1_To_Grid2)(tGrid *grid1, tGrid *grid2, int vind,
                                     int innerdom);
  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("add_const_to_sigma_pm: innerdom is not 0 or 3");

  /* do nothing if dsig=0 */  
  if(dsig==0.0) return;

  printf(" add_const_to_sigma_pm:  dsig=%f  innerdom=%d\n", 
         dsig, innerdom);

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* loop over AnsorgNS boxes */
  for(b=0; b<=3; b++)
  {
    double *surf;
    double *surf2;

    /* do nothing if we are in the wrong box */
    if(b!=innerdom && b!=outerdom) continue;

    surf  = grid->box[b]->v[sigpmi];
    surf2 = grid2->box[b]->v[sigpmi];

    /* add dsig on grid2 to */      
    forallpoints(grid->box[b], i)  surf2[i] = surf[i] + dsig;
  }

  /* initialize coords of grid2 on side of innerdom */
  BNSgrid_init_Coords_pm(grid2, innerdom);
  /* use interpolator that does only side of innerdom */
  Interp_From_Grid1_To_Grid2 = Interp_Var_From_Grid1_To_Grid2_pm;

  /* interpolate some vars from grid onto new grid2 */
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_q"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"),innerdom);
  if( (innerdom==0 && !Getv("BNSdata_rotationstate1","corotation")) ||
      (innerdom==3 && !Getv("BNSdata_rotationstate2","corotation"))   )
  {
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"),innerdom);
  }
  if(interp_qgold)
  {
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);
  }
  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);
}

void smooth_BNSdata_Sigma_NearBoundary(tGrid *grid, int itmax, double tol,
     int (*linear_solver)(tVarList *x, tVarList *b,
              tVarList *r, tVarList *c1,tVarList *c2,
              int itmax, double tol, double *normres,
              void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *),
              void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)))
{
  tGrid *grid_bak;
  tParameter *pdb_bak;
  double dsig_ch = Getd("BNSdata_SmoothSigmaRegion");
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int xi     = Ind("x");
  double dsigp = grid->box[0]->v[sigpmi][0] * dsig_ch;
  double dsigm = grid->box[3]->v[sigpmi][0] * dsig_ch;
  double normresnonlin;
  double xout1, xin1, xout2, xin2;
  int b0n1 = grid->box[0]->n1;
  int b0n2 = grid->box[0]->n2;
  int b3n1 = grid->box[3]->n1;
  int b3n2 = grid->box[3]->n2;
  double sigp_00, sigp_10, sigm_00, sigm_10;

  /* do nothing if no change */
  if(dsig_ch==0.0) return;
  if(Getv("BNSdata_SmoothSigma", "no")) return;

  prdivider(1);
  printf("smooth_BNSdata_Sigma_NearBoundary:  dsigp=%g  dsigm=%g\n",
         dsigp, dsigm);

  /* make exact copies of grid and pdb */
  grid_bak = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid_bak, 0);
  pdb_bak  = make_empty_pdb(npdbmax);
  copy_pdb(pdb, npdb, pdb_bak);

  /* now change Coordinates_AnsorgNS_sigma_pm on grid, so that boundaries
     move inwards */
  add_const_to_sigma_pm(grid, dsigp, 0);
  add_const_to_sigma_pm(grid, dsigm, 3);
  /* get inner and outer edges of both stars */
  xout1 = grid->box[0]->v[xi][0]; /* sigma_p at B=phi=0 */
  xin1  = grid->box[0]->v[xi][Ind_n1n2(0,b0n2-1,0, b0n1,b0n2)]; /* sigma_p at B=1, phi=0 */
  xout2 = grid->box[3]->v[xi][0]; /* sigma_m at B=phi=0 */
  xin2  = grid->box[3]->v[xi][Ind_n1n2(0,b3n2-1,0, b3n1,b3n2)]; /* sigma_m at B=1, phi=0 */
  printf(" -> xin1=%g xout1=%g  xin2=%g xout2=%g\n", xin1,xout1, xin2, xout2);
  sigp_00 = grid->box[0]->v[sigpmi][0]; /* sigma_p at B=phi=0 */
  sigp_10 = grid->box[0]->v[sigpmi][Ind_n1n2(0,b0n2-1,0, b0n1,b0n2)]; /* sigma_p at B=1, phi=0 */
  sigm_00 = grid->box[3]->v[sigpmi][0]; /* sigma_m at B=phi=0 */
  sigm_10 = grid->box[3]->v[sigpmi][Ind_n1n2(0,b3n2-1,0, b3n1,b3n2)]; /* sigma_m at B=1, phi=0 */
  printf(" -> sigp_10=%g sigp_00=%g  sigm_10=%g sigm_00=%g\n",
         sigp_10, sigp_00, sigm_10, sigm_00);

  /* do not touch BNSdata_Sigma in inner boxes, but solve in outer */
  printf("set BNSdata_Sigma by interpolation outside newly shrunk boxes 0 and 3...\n");
  Sets("BNSdata_KeepInnerSigma", "yes");
  BNS_Eqn_Iterator_for_vars_in_string(grid, itmax, tol, 
               &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
  Sets("BNSdata_KeepInnerSigma", "no");

  /* interpolate new BNSdata_Sigma to grid_bak */
  Interp_Var_From_Grid1_To_Grid2_pm(grid, grid_bak, Ind("BNSdata_Sigma"), 0);
  Interp_Var_From_Grid1_To_Grid2_pm(grid, grid_bak, Ind("BNSdata_Sigma"), 3);
  printf("smooth_BNSdata_Sigma_NearBoundary: modified BNSdata_Sigma.\n");

  /* copy grid_bak and pdb_bak back into grid and pdb */
  copy_grid(grid_bak, grid, 0);
  copy_pdb(pdb_bak, npdb, pdb);
  /* free grid_bak pdb_bak */
  free_grid(grid_bak);
  free_pdb(pdb_bak, npdb); 
}

/* add a const to sigma_pm to change surface shape for one star */
void set_sigma_pm_const(tGrid *grid, double sig, int innerdom)
{
  tGrid *grid2;
  int interp_qgold = !Getv("BNSdata_new_q", "FromFields");
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int b, i, outerdom;
  void (*Interp_From_Grid1_To_Grid2)(tGrid *grid1, tGrid *grid2, int vind,
                                     int innerdom);
  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("add_const_to_sigma_pm: innerdom is not 0 or 3");

  printf(" set_sigma_pm_const:  sig=%f  innerdom=%d\n", sig, innerdom);

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* loop over AnsorgNS boxes */
  for(b=0; b<=3; b++)
  {
    double *surf2;

    /* do nothing if we are in the wrong box */
    if(b!=innerdom && b!=outerdom) continue;

    surf2 = grid2->box[b]->v[sigpmi];

    /* set sig on grid2 to */      
    forallpoints(grid->box[b], i)  surf2[i] = sig;
  }

  /* initialize coords of grid2 on side of innerdom */
  BNSgrid_init_Coords_pm(grid2, innerdom);
  /* use interpolator that does only side of innerdom */
  Interp_From_Grid1_To_Grid2 = Interp_Var_From_Grid1_To_Grid2_pm;

  /* interpolate some vars from grid onto new grid2 */
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_q"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"),innerdom);
  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"),innerdom);
  if( (innerdom==0 && !Getv("BNSdata_rotationstate1","corotation")) ||
      (innerdom==3 && !Getv("BNSdata_rotationstate2","corotation"))   )
  {
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"),innerdom);
  }
  if(interp_qgold)
  {
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);
  }
  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);
}

/* find a sigma such that we get a sphere that fits inside the star */
double FindMaxSphere_insideStar(tGrid *grid, int innerdom)
{
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int b, i, outerdom;
  double MaxSphere_sig;

  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("FindMaxSphere_insideStar: innerdom is not 0 or 3");

  /* loop over outer AnsorgNS boxes to find min |sigma_pm| */
  for(b=1; b<=2; b++)
  {
    double *surf;
    /* do nothing if we are in the wrong box */
    if(b!=outerdom) continue;

    surf = grid->box[b]->v[sigpmi];
    MaxSphere_sig = surf[0];

    /* set min |sigma| on grid */      
    forallpoints(grid->box[b], i)
      if(fabs(surf[i])>fabs(MaxSphere_sig)) MaxSphere_sig=surf[i];;
  }
  printf(" FindMaxSphere_insideStar:  MaxSphere_sig=%g  innerdom=%d\n", 
         MaxSphere_sig, innerdom);
  return MaxSphere_sig;
}


/* change boundary temporarily to a sphere, then solve for BNSdata_Sigma
   and interpolate back onto real grid */
void solve_BNSdata_Sigma_WithSphereBoundary(tGrid *grid, int itmax, double tol,
                                            double *normres,
     int (*linear_solver)(tVarList *x, tVarList *b,
              tVarList *r, tVarList *c1,tVarList *c2,
              int itmax, double tol, double *normres,
              void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *),
              void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)))
{
  tGrid *grid_bak;
  tParameter *pdb_bak;
  double dsig_ch = Getd("BNSdata_SmoothSigmaRegion");
  double Sphere_sigma_p = Getd("BNSdata_Sphere_sigma_p");
  double Sphere_sigma_m = Getd("BNSdata_Sphere_sigma_m");
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  int xi     = Ind("x");
  double dsigp = grid->box[0]->v[sigpmi][0] * dsig_ch;
  double dsigm = grid->box[3]->v[sigpmi][0] * dsig_ch;
  double xout1, xin1, xout2, xin2;
  int b0n1 = grid->box[0]->n1;
  int b0n2 = grid->box[0]->n2;
  int b3n1 = grid->box[3]->n1;
  int b3n2 = grid->box[3]->n2;
  double sigp_00, sigp_10, sigm_00, sigm_10;
  double MaxSphere_sigp, MaxSphere_sigm;
  double ds;
  char *BCs_sav;
  char *str;

  prdivider(1);
  printf("solve_BNSdata_Sigma_WithSphereBoundary:\n");
  
  /* get max possible sphere sizes */
  MaxSphere_sigp = FindMaxSphere_insideStar(grid, 0);
  MaxSphere_sigm = FindMaxSphere_insideStar(grid, 3);

  /* if Sphere_sigma_p/m - MaxSphere_sigp/m is about dsigp/m do nothing */
  ds = Sphere_sigma_p - MaxSphere_sigp;
  if(ds<=0.9*dsigp || ds>=1.1*dsigp) Sphere_sigma_p = MaxSphere_sigp + dsigp;
  ds = Sphere_sigma_m - MaxSphere_sigm;
  if(ds>=0.9*dsigm || ds<=1.1*dsigm) Sphere_sigma_m = MaxSphere_sigm + dsigm;
  Setd("BNSdata_Sphere_sigma_p", Sphere_sigma_p);
  Setd("BNSdata_Sphere_sigma_m", Sphere_sigma_m);
  printf(" Sphere_sigma_p=%g  Sphere_sigma_m=%g\n",
         Sphere_sigma_p, Sphere_sigma_m);

  /* make exact copies of grid and pdb */
  grid_bak = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid_bak, 0);
  pdb_bak  = make_empty_pdb(npdbmax);
  copy_pdb(pdb, npdb, pdb_bak);

  /* now change Coordinates_AnsorgNS_sigma_pm on grid, so that boundaries
     move inwards */
  set_sigma_pm_const(grid, Sphere_sigma_p, 0);
  set_sigma_pm_const(grid, Sphere_sigma_m, 3);
  /* get inner and outer edges of both spheres */
  xout1 = grid->box[0]->v[xi][0]; /* sigma_p at B=phi=0 */
  xin1  = grid->box[0]->v[xi][Ind_n1n2(0,b0n2-1,0, b0n1,b0n2)]; /* sigma_p at B=1, phi=0 */
  xout2 = grid->box[3]->v[xi][0]; /* sigma_m at B=phi=0 */
  xin2  = grid->box[3]->v[xi][Ind_n1n2(0,b3n2-1,0, b3n1,b3n2)]; /* sigma_m at B=1, phi=0 */
  printf(" -> xin1=%g xout1=%g  xin2=%g xout2=%g\n", xin1,xout1, xin2, xout2);
  sigp_00 = grid->box[0]->v[sigpmi][0]; /* sigma_p at B=phi=0 */
  sigp_10 = grid->box[0]->v[sigpmi][Ind_n1n2(0,b0n2-1,0, b0n1,b0n2)]; /* sigma_p at B=1, phi=0 */
  sigm_00 = grid->box[3]->v[sigpmi][0]; /* sigma_m at B=phi=0 */
  sigm_10 = grid->box[3]->v[sigpmi][Ind_n1n2(0,b3n2-1,0, b3n1,b3n2)]; /* sigma_m at B=1, phi=0 */
  printf(" -> sigp_10=%g sigp_00=%g  sigm_10=%g sigm_00=%g\n",
         sigp_10, sigp_00, sigm_10, sigm_00);

  /* do not touch BNSdata_Sigma in inner boxes, but solve in outer */
  printf("solve for BNSdata_Sigma with new shrunk boundary...\n");
  BCs_sav = cmalloc( strlen(Gets("BNSdata_Sigma_surface_BCs"))+10 );
  str     = cmalloc( strlen(BCs_sav)+100 );
  strcpy(BCs_sav, Gets("BNSdata_Sigma_surface_BCs"));
  sprintf(str, "%s %s", BCs_sav, "EllEqn");
  Sets("BNSdata_Sigma_surface_BCs", str); 
  BNS_Eqn_Iterator_for_vars_in_string(grid, itmax, tol, normres, 
                                      linear_solver, 1, "BNSdata_Sigma");
  Sets("BNSdata_Sigma_surface_BCs", BCs_sav);
  free(str);
  free(BCs_sav);

  /* interpolate new BNSdata_Sigma to grid_bak */
  Interp_Var_From_Grid1_To_Grid2_pm(grid, grid_bak, Ind("BNSdata_Sigma"), 0);
  Interp_Var_From_Grid1_To_Grid2_pm(grid, grid_bak, Ind("BNSdata_Sigma"), 3);
  printf("solve_BNSdata_Sigma_WithSphereBoundary: done.\n");
  prdivider(1);

  /* copy grid_bak and pdb_bak back into grid and pdb */
  copy_grid(grid_bak, grid, 0);
  copy_pdb(pdb_bak, npdb, pdb);
  /* free grid_bak pdb_bak */
  free_grid(grid_bak);
  free_pdb(pdb_bak, npdb);
}

/* filter Coordinates_AnsorgNS_sigma_pm and adjust the shape of the 
   boundary between domain0/1 or domain3/2 accordingly */
void filter_Coordinates_AnsorgNS_sigma_pm(tGrid *grid, int innerdom)
{
  tGrid *grid2;
  int interp_qgold = !Getv("BNSdata_new_q", "FromFields");
  int outerdom;
  int sigpmi = Ind("Coordinates_AnsorgNS_sigma_pm");
  void (*Interp_From_Grid1_To_Grid2)(tGrid *grid1, tGrid *grid2, int vind,
                                     int innerdom);
  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("filter_Coordinates_AnsorgNS_sigma_pm: "
                 "innerdom is not 0 or 3");

  /* do nothing if filter is off */
  if(Getv("BNSdata_domainshape_filter", "no")) return;
  printf("filter_Coordinates_AnsorgNS_sigma_pm: innerdom=%d outerdom=%d\n", 
         innerdom, outerdom);

  if(Getv("BNSdata_domainshape_filter", "LowPassInB"))
    printf(" LowPassInB with nfB=%d inside"
           " reset_Coordinates_AnsorgNS_sigma_pm\n",
           Geti("BNSdata_domainshape_filter_nfB"));
  if(Getv("BNSdata_domainshape_filter", "LowPassInB_dsigma_pm_dB_01_EQ_0"))
    printf(" dsigma_pm_dB01=0 on grid with n2=%d in"
           " reset_Coordinates_AnsorgNS_sigma_pm\n",
           Geti("BNSdata_domainshape_filter_nfB"));
  if(Getv("BNSdata_domainshape_filter", "keep_sigma_pm_B1"))
    printf(" kept previous sigma_pm at B=1 inside"
           " reset_Coordinates_AnsorgNS_sigma_pm\n");
  if(Getv("BNSdata_domainshape_filter", "min_dsigma_pm_dB_1"))
    printf(" minimized dsigma_pm_dB at B=1 inside"
           " reset_Coordinates_AnsorgNS_sigma_pm\n");
  if(Getv("BNSdata_domainshape_filter", "dsigma_pm_dB_01_EQ_0"))
    printf(" set dsigma_pm_dB=0 at B=0,1 inside"
           " reset_Coordinates_AnsorgNS_sigma_pm\n");

  /* filter Coordinates_AnsorgNS_sigma_pm on grid2 */
  if(Getv("BNSdata_domainshape_filter", "Bphi2/3"))
  {
    /* make new grid2, which is an exact copy of grid */
    grid2 = make_empty_grid(grid->nvariables, 0);
    copy_grid(grid, grid2, 0);

    /* filter with 2/3 rule */
    BNSdata_filter_with2o3rule_inBphi(grid2, sigpmi, innerdom);

    /* initialize coords of grid2 on side of innerdom */
    BNSgrid_init_Coords_pm(grid2, innerdom);
    /* use interpolator that does only side of innerdom */
    Interp_From_Grid1_To_Grid2 = Interp_Var_From_Grid1_To_Grid2_pm;

    /* interpolate some vars from grid onto new grid2 */
    //  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"),innerdom);
    //  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"),innerdom);
    Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"),innerdom);
    if( (innerdom==0 && !Getv("BNSdata_rotationstate1","corotation")) ||
        (innerdom==3 && !Getv("BNSdata_rotationstate2","corotation"))   )
    {
      Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"),innerdom);
      Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"),innerdom);
      Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"),innerdom);
      Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"),innerdom);
    }
    if(interp_qgold)
      Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);

    BNS_compute_new_centered_q(grid2);

//  /* set q to zero if q<0 or in region 1 and 2 */
//  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

    /* copy grid2 back into grid, and free grid2 */
    copy_grid(grid2, grid, 0);
    free_grid(grid2);
  }
}

/* Adjust C1/2 and thus q by demanding that m01 and m02 stay the same. */
int adjust_C1_C2_q_keep_restmasses(tGrid *grid, int it, double tol)
{
  double Cvec[3];
  double m0errorvec[3];
  double m01, m02, m0_error, dm01, dm02;
  int check, stat, bi, i;
  tGrid *grid_bak;
  tParameter *pdb_bak;
  t_grid_grid0_m01_m02_struct pars[1];

  /* grid and pdb for backups */
  grid_bak = make_empty_grid(grid->nvariables, 0);
  pdb_bak  = make_empty_pdb(npdbmax);

  /* rest masses before adjusting q */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("adjust_C1_C2_q_keep_restmasses: in BNSdata_solve step %d:\n", it);
  printf(" rest mass in inner domains before computing new q:"
         " m01=%g m02=%g\n", m01, m02);

  /* compute error in masses */
  dm01 = m01 - Getd("BNSdata_m01");
  dm02 = m02 - Getd("BNSdata_m02");
  m0_error = dm01*dm01 + dm02*dm02;
  m0_error = sqrt(m0_error)/(Getd("BNSdata_m01")+Getd("BNSdata_m02"));
  printf(" rest mass error = %.4e "
         "(before adjusting C1/2)\n", m0_error);

  /* set desired masses for this iteration */
  pars->m01 = Getd("BNSdata_m01"); // + dm01*0.9;
  pars->m02 = Getd("BNSdata_m02"); // + dm02*0.9;
  printf(" adjusting q,C1,C2 to achieve: m01=%g  m02=%g  tol*0.01=%g\n",
         pars->m01, pars->m02, tol*0.01);

  /* print C1/2 we used before */
  printf(" old: BNSdata_C1=%g BNSdata_C2=%g\n",
         Getd("BNSdata_C1"), Getd("BNSdata_C2"));
  /* printf("     => m01=%g m02=%g\n", m01, m02); */

  /* adjust C1/2 */
  if(!Getv("BNSdata_adjust_C1C2", "no"))
  {
    /* refine guesses for C1,C2? */
    if(Getv("BNSdata_adjust_C1C2", "refineguess"))
    {
      /* choose C1/2 such that rest masses are not too big or too small */
      for(i=0; i<1000; i++)
      {
        double *q_b1 = grid->box[1]->v[Ind("BNSdata_q")];
        double *q_b2 = grid->box[2]->v[Ind("BNSdata_q")];

        BNS_compute_new_centered_q(grid);
        m01 = GetInnerRestMass(grid, 0);
        m02 = GetInnerRestMass(grid, 3);

        check = 0;

        if(m01 > 1.1*Getd("BNSdata_m01"))
          { Setd("BNSdata_C1", 0.999*Getd("BNSdata_C1"));  check=1; }
        else if(m01 < 0.9*Getd("BNSdata_m01"))
          { Setd("BNSdata_C1", 1.002*Getd("BNSdata_C1"));  check=1; }
    
        if(m02 > 1.1*Getd("BNSdata_m02") && Getd("BNSdata_m02")>0)
          { Setd("BNSdata_C2", 0.999*Getd("BNSdata_C2"));  check=1; }
        else if(m02 < 0.9*Getd("BNSdata_m02") && Getd("BNSdata_m02")>0)
          { Setd("BNSdata_C2", 1.002*Getd("BNSdata_C2"));  check=1; }

        if(check==0) break;
      }

      /* refine guess for C1/2 */
      pars->grid = grid;
      Cvec[1] = Getd("BNSdata_C1");
      stat = newton_linesrch_itsP(Cvec, 1, &check, m01_guesserror_VectorFuncP,
                                  (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
      if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
      Setd("BNSdata_C1", Cvec[1]);

      Cvec[1] = Getd("BNSdata_C2");
      if(Getd("BNSdata_m02")>0)
        stat = newton_linesrch_itsP(Cvec, 1, &check, m02_guesserror_VectorFuncP,
                                    (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
      if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
      Setd("BNSdata_C2", Cvec[1]);

      /* print guess for C1/2 */                                        
      printf(" guess: BNSdata_C1=%g BNSdata_C2=%g\n",
             Getd("BNSdata_C1"), Getd("BNSdata_C2"));
    }
    /***********************************************************************/
    /* do newton_linesrch_itsP iterations of Cvec until m0errorvec is zero */
    /***********************************************************************/
    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("BNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* adjust C1 and thus m01 */
    Cvec[1] = Getd("BNSdata_C1");
    stat = newton_linesrch_itsP(Cvec, 1, &check, m01_error_VectorFuncP,
                                (void *) pars, 1000, tol*0.01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
    Setd("BNSdata_C1", Cvec[1]);

    /* filter domain shape to keep star1's surface smooth during iterations */
    filter_Coordinates_AnsorgNS_sigma_pm(grid, 0);

    /* see if we keep the new domain shape and C1 */
    restore_grid_pdb_if_change_in_star_is_large(1,grid,pdb, grid_bak,pdb_bak);

    /* average domain shape of star1 on grid and grid_bak */
    average_current_and_old_surfaceshape(Getd("BNSdata_domainshape_weight"),
                                         0, grid,pdb, grid_bak,pdb_bak);

    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("BNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* adjust C2 and thus m02 */
    Cvec[1] = Getd("BNSdata_C2");
    if(Getd("BNSdata_m02")>0)
      stat = newton_linesrch_itsP(Cvec, 1, &check, m02_error_VectorFuncP,
                                  (void *) pars, 1000, tol*0.01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
    Setd("BNSdata_C2", Cvec[1]);

    /* filter domain shape to keep star2's surface smooth during iterations */
    filter_Coordinates_AnsorgNS_sigma_pm(grid, 3);

    /* see if we keep the new domain shape and C2 */
    restore_grid_pdb_if_change_in_star_is_large(2, grid,pdb, grid_bak,pdb_bak);

    /* average domain shape of star2 on grid and grid_bak */
    average_current_and_old_surfaceshape(Getd("BNSdata_domainshape_weight"),
                                         3, grid,pdb, grid_bak,pdb_bak);

    printf("adjust_C1_C2_q_keep_restmasses:\n");
    printf(" new: BNSdata_C1=%g BNSdata_C2=%g\n",
           Getd("BNSdata_C1"), Getd("BNSdata_C2"));
  }
  else /* keep C1/2, but adjust q */
  {
    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("BNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* compute new q in star1 */
    compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, pars->grid0, 0);

    /* filter domain shape to keep star1's surface smooth during iterations */
    filter_Coordinates_AnsorgNS_sigma_pm(grid, 0);

    /* see if we keep the new domain shape and C1 */
    restore_grid_pdb_if_change_in_star_is_large(1,grid,pdb, grid_bak,pdb_bak);

    /* average domain shape of star1 on grid and grid_bak */
    average_current_and_old_surfaceshape(Getd("BNSdata_domainshape_weight"),
                                         0, grid,pdb, grid_bak,pdb_bak);

    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("BNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* compute new q in star2 */
    compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, pars->grid0, 3);
  
    /* filter domain shape to keep star2's surface smooth during iterations */
    filter_Coordinates_AnsorgNS_sigma_pm(grid, 3);

    /* see if we keep the new domain shape and C2 */
    restore_grid_pdb_if_change_in_star_is_large(2, grid,pdb, grid_bak,pdb_bak);

    /* average domain shape of star2 on grid and grid_bak */
    average_current_and_old_surfaceshape(Getd("BNSdata_domainshape_weight"),
                                         3, grid,pdb, grid_bak,pdb_bak);

    printf("adjust_C1_C2_q_keep_restmasses: adjusted q, but kept C1/2\n");
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* print new masses */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  printf("     => m01=%.19g m02=%.19g\n", m01, m02);

  /* free grid and pdb for backups */
  free_pdb(pdb_bak, npdb);
  free_grid(grid_bak);
        
  return 0;
}

/* Adjust C1/2: Try adjustment for several Omega and possibly x_CM
   and choose the one with the smallest qmax difference...  */
int adjust_C1_C2_Omega_xCM_q_WT(tGrid *grid, int it, double tol, 
                                double *dOmega)
{
  double Omega;
  int bi1, bi2;
  double qmax1, Xmax1, Ymax1, qmax1sav, diff1_m, diff1_0, diff1_p;
  double qmax2, Xmax2, Ymax2, qmax2sav, diff2_m, diff2_0, diff2_p;
  double dif_m, dif_0, dif_p;

  /* save Omega */
  Omega = Getd("BNSdata_Omega");
  printf("BNSdata_solve step %d: old Omega = %.4e  *dOmega = %.4e\n",
         it, Omega, *dOmega);
         
  /* save old qg in BNSdata_qgold */
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  /* save both qmax */
  qmax1sav = qmax1;
  qmax2sav = qmax2;

  /* compute diff between new qmax1 and qmax1sav for Omega - *dOmega */
  Setd("BNSdata_Omega", Omega - *dOmega);
  printf("BNSdata_solve step %d: get qmax diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);
  diff1_m = qmax1 - qmax1sav;
  diff2_m = qmax2 - qmax2sav;
  dif_m = sqrt(diff1_m*diff1_m + diff2_m*diff2_m);

  /* compute diff between new qmax1 and qmax1sav for Omega + *dOmega */
  Setd("BNSdata_Omega", Omega + *dOmega);
  printf("BNSdata_solve step %d: get qmax diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);
  diff1_p = qmax1 - qmax1sav;
  diff2_p = qmax2 - qmax2sav;
  dif_p = sqrt(diff1_p*diff1_p + diff2_p*diff2_p);

  /* compute diff between new qmax1 and qmax1sav for Omega */
  Setd("BNSdata_Omega", Omega);
  printf("BNSdata_solve step %d: get qmax diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);
  diff1_0 = qmax1 - qmax1sav;
  diff2_0 = qmax2 - qmax2sav;
  dif_0 = sqrt(diff1_0*diff1_0 + diff2_0*diff2_0);

  printf("BNSdata_solve step %d: dif_m=%g\n", it, dif_m);
  printf("BNSdata_solve step %d: dif_0=%g\n", it, dif_0);
  printf("BNSdata_solve step %d: dif_p=%g\n", it, dif_p);
  if(dif_0<=dif_p && dif_0<=dif_m)
  { Setd("BNSdata_Omega", Omega);  *dOmega=*dOmega*0.5; }
  if(dif_p<dif_0  && dif_p<=dif_m) 
    Setd("BNSdata_Omega", Omega+*dOmega);
  if(dif_m<dif_0  && dif_m<dif_p)
    Setd("BNSdata_Omega", Omega-*dOmega);
  printf("BNSdata_solve step %d: new Omega = %.4e  *dOmega = %.4e\n",
         it, Getd("BNSdata_Omega"), *dOmega);

  /* set new Omega, and compute new q */
  if( !dequal(Getd("BNSdata_Omega"), Omega) )
    adjust_C1_C2_q_keep_restmasses(grid, it, tol);

  return 0;
}

/* Adjust C1/2: Try adjustment for several Omega and possibly x_CM
   and choose the one with the smallest L2 q difference...  */
int adjust_C1_C2_Omega_xCM_q_WT_L2(tGrid *grid, int it, double tol, 
                                   double *dOmega)
{
  double Omega;
  double dif_m, dif_0, dif_p;

  /* save Omega */
  Omega = Getd("BNSdata_Omega");
  printf("BNSdata_solve step %d: old Omega = %.4e  *dOmega = %.4e\n",
         it, Omega, *dOmega);
         
  /* save old qg in BNSdata_qgold */
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* compute L2-diff between new q and qold for Omega - *dOmega */
  Setd("BNSdata_Omega", Omega - *dOmega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_m = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega + *dOmega */
  Setd("BNSdata_Omega", Omega + *dOmega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_p = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega */
  Setd("BNSdata_Omega", Omega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  BNS_compute_new_centered_q(grid);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_0 = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* set new Omega */
  printf("BNSdata_solve step %d: dif_m=%g\n", it, dif_m);
  printf("BNSdata_solve step %d: dif_0=%g\n", it, dif_0);
  printf("BNSdata_solve step %d: dif_p=%g\n", it, dif_p);
  if(dif_0<=dif_p && dif_0<=dif_m)
  { Setd("BNSdata_Omega", Omega);  *dOmega=*dOmega*0.5; }
  if(dif_p<dif_0  && dif_p<=dif_m) 
    Setd("BNSdata_Omega", Omega+*dOmega);
  if(dif_m<dif_0  && dif_m<dif_p)
    Setd("BNSdata_Omega", Omega-*dOmega);
  printf("BNSdata_solve step %d: new Omega = %.4e  *dOmega = %.4e\n",
         it, Getd("BNSdata_Omega"), *dOmega);

  /* compute new q */
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);

  return 0;
}

/* Adjust C1/2: Try adjustment for several Omega and possibly x_CM
   and choose the one with the smallest L2 q difference after
   we adjust C1/2...  */
int adjust_C1_C2_Omega_xCM_q_min_qchange(tGrid *grid, int it, double tol, 
                                         double *dOmega)
{
  double Omega;
  double dif_m, dif_0, dif_p;

errorexit("in order to use adjust_C1_C2_Omega_xCM_q_min_qchange you must uncomment\n"
          "//Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind(\"BNSdata_qgold\"));\n"
          "in the function compute_new_q_and_adjust_domainshapes");

  /* save Omega */
  Omega = Getd("BNSdata_Omega");
  printf("BNSdata_solve step %d: old Omega = %.4e  *dOmega = %.4e\n",
         it, Omega, *dOmega);
         
  /* save old q in BNSdata_qgold */
  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* compute L2-diff between new q and qold for Omega - *dOmega */
  Setd("BNSdata_Omega", Omega - *dOmega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_m = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega + *dOmega */
  Setd("BNSdata_Omega", Omega + *dOmega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_p = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* compute L2-diff between new q and qold for Omega */
  Setd("BNSdata_Omega", Omega);
  printf("BNSdata_solve step %d: get q L2-diff. for Omega = %g\n",
         it, Getd("BNSdata_Omega"));
  adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  varadd(grid, Ind("BNSdata_temp1"), 
             1,Ind("BNSdata_qg"), -1,Ind("BNSdata_qgold"));
  dif_0 = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
          varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1"));

  /* set new Omega */
  printf("BNSdata_solve step %d: dif_m=%g\n", it, dif_m);
  printf("BNSdata_solve step %d: dif_0=%g\n", it, dif_0);
  printf("BNSdata_solve step %d: dif_p=%g\n", it, dif_p);
  if(dif_0<=dif_p && dif_0<=dif_m)
  { Setd("BNSdata_Omega", Omega);  *dOmega=*dOmega*0.5; }
  if(dif_p<dif_0  && dif_p<=dif_m) 
    Setd("BNSdata_Omega", Omega+*dOmega);
  if(dif_m<dif_0  && dif_m<dif_p)
    Setd("BNSdata_Omega", Omega-*dOmega);
  printf("BNSdata_solve step %d: new Omega = %.4e  *dOmega = %.4e\n",
         it, Getd("BNSdata_Omega"), *dOmega);

  /* compute new q */
  if( !dequal(Getd("BNSdata_Omega"), Omega) )
    adjust_C1_C2_q_keep_restmasses(grid, it, tol);

  return 0;
}

/* for newton_linesrch_its: compute error in xouts */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void xouts_error_VectorFunc(int n, double *vec, double *fvec)
{
  double xout1, xout2;
  int xind = Ind("x");
  tBox *box1 = xouts_error_VectorFunc__grid->box[0];
  tBox *box3 = xouts_error_VectorFunc__grid->box[3];

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);
  printf("xouts_error_VectorFunc: Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  adjust_C1_C2_q_keep_restmasses(xouts_error_VectorFunc__grid,
                                 xouts_error_VectorFunc__it,
                                 xouts_error_VectorFunc__tol);
  xout1 = box1->v[xind][0];
  xout2 = box3->v[xind][0];
  printf("xouts_error_VectorFunc: Omega=%g x_CM=%g xout1=%g xout2=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"),
         xout1, xout2);  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = xout1 - xouts_error_VectorFunc__xout1;
  if(n>=2) fvec[2] = xout2 - xouts_error_VectorFunc__xout2;
}

/* Adjust Omega and x_CM so that outer edges xout1/2 stay put */
int adjust_C1_C2_Omega_xCM_q_keep_xout(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int xind = Ind("x");
  double OmxCMvec[3];
  double dxout_m0[3];
  double dxout_00[3];
  double dxout_p0[3];
  double dxout_0m[3];
  double dxout_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int keepone   = Getv("BNSdata_adjust", "keep_one_xout");

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_C1_C2_Omega_xCM_q_keep_xout: in BNSdata_solve step %d\n"
         "adjust_C1_C2_Omega_xCM_q_keep_xout: old Omega = %g  x_CM = %g\n",
         it, Omega, x_CM);

  /* set global vars */
  xouts_error_VectorFunc__grid  = grid;
  xouts_error_VectorFunc__it    = it;
  xouts_error_VectorFunc__tol   = tol;
  xouts_error_VectorFunc__xout1 = grid->box[0]->v[xind][0];
  xouts_error_VectorFunc__xout2 = grid->box[3]->v[xind][0];
  printf("adjust_C1_C2_Omega_xCM_q_keep_xout: old xout1 = %g  xout2 = %g\n",
         xouts_error_VectorFunc__xout1, xouts_error_VectorFunc__xout2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    xouts_error_VectorFunc(2, OmxCMvec, dxout_m0);
    OmxCMvec[1] = Omega + dOmega;
    xouts_error_VectorFunc(2, OmxCMvec, dxout_p0);
    OmxCMvec[1] = Omega;
    /* xouts_error_VectorFunc(2, OmxCMvec, dxout_00); */
    printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_m0[1]=%g dxout_m0[2]=%g\n",
           dxout_m0[1], dxout_m0[2]);
    /* printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_00[1]=%g dxout_00[2]=%g\n",
           dxout_00[1], dxout_00[2]); */
    printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_p0[1]=%g dxout_p0[2]=%g\n",
           dxout_p0[1], dxout_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxout_m0[1]*dxout_p0[1]<0.0) && (dxout_m0[2]*dxout_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxout_m0[1]*dxout_p0[1]<0.0) || (dxout_m0[2]*dxout_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xout to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxout_m0[1]*dxout_p0[1]>=0.0)
      {
        xouts_error_VectorFunc(2, OmxCMvec, dxout_00);
        xouts_error_VectorFunc__xout1 = grid->box[0]->v[xind][0];
      }
      if(dxout_m0[2]*dxout_p0[2]>=0.0)
      {
        xouts_error_VectorFunc(2, OmxCMvec, dxout_00);
        xouts_error_VectorFunc__xout2 = grid->box[3]->v[xind][0];
      }
      do_lnsrch = 1;
      printf("adjust_C1_C2_Omega_xCM_q_keep_xout: "
             "changed: old xout1 = %g  xout2 = %g\n",
             xouts_error_VectorFunc__xout1, xouts_error_VectorFunc__xout2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      xouts_error_VectorFunc(2, OmxCMvec, dxout_0m);
      OmxCMvec[2] = x_CM + dx_CM;
      xouts_error_VectorFunc(2, OmxCMvec, dxout_0p);
      OmxCMvec[2] = x_CM;
      /* xouts_error_VectorFunc(2, OmxCMvec, dxout_00); */
      printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_0m[1]=%g dxout_0m[2]=%g\n",
             dxout_m0[1], dxout_m0[2]);
      /* printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_00[1]=%g dxout_00[2]=%g\n",
             dxout_00[1], dxout_00[2]); */
      printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_0p[1]=%g dxout_0p[2]=%g\n",
             dxout_p0[1], dxout_p0[2]);
  
      if( (dxout_m0[1]*dxout_p0[1]<0.0) && (dxout_0m[2]*dxout_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxout_m0[2]*dxout_p0[2]<0.0) && (dxout_0m[1]*dxout_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_C1_C2_Omega_xCM_q_keep_xout: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xout1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_its(OmxCMvec, 2, &check, xouts_error_VectorFunc,
                               1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    xouts_error_VectorFunc(2, OmxCMvec, dxout_00);
    printf("adjust_C1_C2_Omega_xCM_q_keep_xout: dxout_00[1] = %g  dxout_00[2] = %g\n",
           dxout_00[1], dxout_00[2]);
  }
  printf("adjust_C1_C2_Omega_xCM_q_keep_xout: new Omega = %g  x_CM = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_itsP: compute error in xins */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void xins_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double xin1, xin2;
  int xind = Ind("x");
  tBox *box0;
  tBox *box3;
  t_grid_x1_2_struct *pars;
  int n1,n2,n3;

  /* get pars */
  pars = (t_grid_x1_2_struct *) p;
  box0 = pars->grid->box[0];      
  box3 = pars->grid->box[3];

  /* set n1,n2,n3 */
  n1 = box0->n1;
  n2 = box0->n2;
  n3 = box0->n3;

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);
  printf("xins_error_VectorFuncP: Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  adjust_C1_C2_q_keep_restmasses(pars->grid, pars->it, pars->tol);
  xin1 = box0->v[xind][Index(0,n2-1,0)];
  xin2 = box3->v[xind][Index(0,n2-1,0)];
  printf("xins_error_VectorFuncP: Omega=%g x_CM=%g xin1=%g xin2=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"),
         xin1, xin2);  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = xin1 - pars->x1;
  if(n>=2) fvec[2] = xin2 - pars->x2;
}

/* Adjust Omega and x_CM so that inner edges xin1/2 stay put */
int adjust_C1_C2_Omega_xCM_q_keep_xin(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int xind = Ind("x");
  double OmxCMvec[3];
  double dxin_m0[3];
  double dxin_00[3];
  double dxin_p0[3];
  double dxin_0m[3];
  double dxin_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int keepone   = Getv("BNSdata_adjust", "keep_one_xin");
  t_grid_x1_2_struct pars[1]; /* pars array */
  int n1 = grid->box[0]->n1;
  int n2 = grid->box[0]->n2;
  int n3 = grid->box[0]->n3;

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_C1_C2_Omega_xCM_q_keep_xin: in BNSdata_solve step %d\n"
         "adjust_C1_C2_Omega_xCM_q_keep_xin: old Omega = %g  x_CM = %g\n",
         it, Omega, x_CM);

  /* set global vars */
  pars->grid  = grid;
  pars->it    = it;
  pars->tol   = tol;
  pars->x1 = grid->box[0]->v[xind][Index(0,n2-1,0)];
  pars->x2 = grid->box[3]->v[xind][Index(0,n2-1,0)];
  printf("adjust_C1_C2_Omega_xCM_q_keep_xin: old xin1 = %g  xin2 = %g\n",
         pars->x1, pars->x2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    xins_error_VectorFuncP(2, OmxCMvec, dxin_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    xins_error_VectorFuncP(2, OmxCMvec, dxin_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* xins_error_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars); */
    printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_m0[1]=%g dxin_m0[2]=%g\n",
           dxin_m0[1], dxin_m0[2]);
    /* printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_00[1]=%g dxin_00[2]=%g\n",
           dxin_00[1], dxin_00[2]); */
    printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_p0[1]=%g dxin_p0[2]=%g\n",
           dxin_p0[1], dxin_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxin_m0[1]*dxin_p0[1]<0.0) && (dxin_m0[2]*dxin_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxin_m0[1]*dxin_p0[1]<0.0) || (dxin_m0[2]*dxin_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xin to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxin_m0[1]*dxin_p0[1]>=0.0)
      {
        xins_error_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
        pars->x1 = grid->box[0]->v[xind][Index(0,n2-1,0)];
      }
      if(dxin_m0[2]*dxin_p0[2]>=0.0)
      {
        xins_error_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
        pars->x2 = grid->box[3]->v[xind][Index(0,n2-1,0)];
      }
      do_lnsrch = 1;
      printf("adjust_C1_C2_Omega_xCM_q_keep_xin: "
             "changed: old xin1 = %g  xin2 = %g\n",
             pars->x1, pars->x2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      xins_error_VectorFuncP(2, OmxCMvec, dxin_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      xins_error_VectorFuncP(2, OmxCMvec, dxin_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* xins_error_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars); */
      printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_0m[1]=%g dxin_0m[2]=%g\n",
             dxin_m0[1], dxin_m0[2]);
      /* printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_00[1]=%g dxin_00[2]=%g\n",
             dxin_00[1], dxin_00[2]); */
      printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_0p[1]=%g dxin_0p[2]=%g\n",
             dxin_p0[1], dxin_p0[2]);
  
      if( (dxin_m0[1]*dxin_p0[1]<0.0) && (dxin_0m[2]*dxin_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxin_m0[2]*dxin_p0[2]<0.0) && (dxin_0m[1]*dxin_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_C1_C2_Omega_xCM_q_keep_xin: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xin1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, xins_error_VectorFuncP,
                                (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    xins_error_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
    printf("adjust_C1_C2_Omega_xCM_q_keep_xin: dxin_00[1] = %g  dxin_00[2] = %g\n",
           dxin_00[1], dxin_00[2]);
  }
  printf("adjust_C1_C2_Omega_xCM_q_keep_xin: new Omega = %g  x_CM = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_its: compute error in xouts, but without
   changing domain shapes, i.e. for given C1/2 */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void q_at_xout1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double X,Y,Z;
  int bi;
  int blist[2];
  double q1, q2;
  t_grid_x1_2_struct *pars;

  /* get pars */
  pars = (t_grid_x1_2_struct *) p;

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);
  printf("q_at_xout1_2_VectorFuncP: Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));

  /* compute q at xout1 */
  X=0.001;     Y=Z=0;
  blist[0]=0;  blist[1]=1;
  bi = b_X_of_x_forgiven_YZ_inboxlist(pars->grid, blist,2, &X, pars->x1, Y,Z);
  q1 = BNS_compute_new_centered_q_atXYZ(pars->grid, bi, X,Y,Z);
  // printf("  xout1=%g is at bi=%d (X,Y,Z)=(%g,%g,%g)\n", pars->x1, bi,X,Y,Z);

  /* compute q at xout2 */
  X=0.001;     Y=Z=0;
  blist[0]=3;  blist[1]=2;
  bi = b_X_of_x_forgiven_YZ_inboxlist(pars->grid, blist,2, &X, pars->x2, Y,Z);
  q2 = BNS_compute_new_centered_q_atXYZ(pars->grid, bi, X,Y,Z);
  // printf("  xout2=%g is at bi=%d (X,Y,Z)=(%g,%g,%g)\n", pars->x2, bi,X,Y,Z);

  printf("q_at_xout1_2_VectorFuncP: q(xout1)=%g q(xout2)=%g\n", q1, q2);
  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = q1;
  if(n>=2) fvec[2] = q2;
}

/* Adjust Omega and x_CM so that outer edges xout1/2 stay put */
int adjust_Omega_xCM_q_fix_xout(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int xind = Ind("x");
  double OmxCMvec[3];
  double dxout_m0[3];
  double dxout_00[3];
  double dxout_p0[3];
  double dxout_0m[3];
  double dxout_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int keepone   = Getv("BNSdata_adjust", "keep_one_xout");
  t_grid_x1_2_struct pars[1]; /* pars array */

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_q_fix_xout: in BNSdata_solve step %d\n"
         "adjust_Omega_xCM_q_fix_xout: old Omega = %g  x_CM = %g\n",
         it, Omega, x_CM);

  /* set global vars */
  pars->grid = grid;
  pars->x1 = grid->box[0]->v[xind][0];
  pars->x2 = grid->box[3]->v[xind][0];
  printf("adjust_Omega_xCM_q_fix_xout: old xout1 = %g  xout2 = %g\n",
         pars->x1, pars->x2);
  prdivider(0);
  adjust_C1_C2_q_keep_restmasses(grid, it, tol*100.0);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_00, (void *) pars); */
    printf("adjust_Omega_xCM_q_fix_xout: dxout_m0[1]=%g dxout_m0[2]=%g\n",
           dxout_m0[1], dxout_m0[2]);
    /* printf("adjust_Omega_xCM_q_fix_xout: dxout_00[1]=%g dxout_00[2]=%g\n",
           dxout_00[1], dxout_00[2]); */
    printf("adjust_Omega_xCM_q_fix_xout: dxout_p0[1]=%g dxout_p0[2]=%g\n",
           dxout_p0[1], dxout_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxout_m0[1]*dxout_p0[1]<0.0) && (dxout_m0[2]*dxout_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxout_m0[1]*dxout_p0[1]<0.0) || (dxout_m0[2]*dxout_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xout to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxout_m0[1]*dxout_p0[1]>=0.0)
      {
        q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_00, (void *) pars);
        pars->x1 = grid->box[0]->v[xind][0];
      }
      if(dxout_m0[2]*dxout_p0[2]>=0.0)
      {
        q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_00, (void *) pars);
        pars->x2 = grid->box[3]->v[xind][0];
      }
      do_lnsrch = 1;
      printf("adjust_Omega_xCM_q_fix_xout: "
             "changed: old xout1 = %g  xout2 = %g\n",
             pars->x1, pars->x2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_00, (void *) pars); */
      printf("adjust_Omega_xCM_q_fix_xout: dxout_0m[1]=%g dxout_0m[2]=%g\n",
             dxout_m0[1], dxout_m0[2]);
      /* printf("adjust_Omega_xCM_q_fix_xout: dxout_00[1]=%g dxout_00[2]=%g\n",
             dxout_00[1], dxout_00[2]); */
      printf("adjust_Omega_xCM_q_fix_xout: dxout_0p[1]=%g dxout_0p[2]=%g\n",
             dxout_p0[1], dxout_p0[2]);
  
      if( (dxout_m0[1]*dxout_p0[1]<0.0) && (dxout_0m[2]*dxout_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxout_m0[2]*dxout_p0[2]<0.0) && (dxout_0m[1]*dxout_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_q_fix_xout: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xout1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, q_at_xout1_2_VectorFuncP,
                                (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    q_at_xout1_2_VectorFuncP(2, OmxCMvec, dxout_00, (void *) pars);
    printf("adjust_Omega_xCM_q_fix_xout: dxout_00[1] = %g  dxout_00[2] = %g\n",
           dxout_00[1], dxout_00[2]);
  }
  printf("adjust_Omega_xCM_q_fix_xout: new Omega = %g  x_CM = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_itsP: compute error in xins without changing domain
   shapes, but adjust C1/2 using m01/2_guesserror_VectorFuncP */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void q_at_xin1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double X,Y,Z;
  int bi;
  int blist[2];
  double q1, q2;
  t_grid_x1_2_struct *pars; /* pointer to access par p that is passed in */
        
  /* get pars */
  pars = (t_grid_x1_2_struct *) p;

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);
  printf("q_at_xin1_2_VectorFuncP: Omega=%.13g x_CM=%.13g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));

  /* adjust C1/2 with m01/2_guesserror_VectorFuncP */
  if(0)
  {
    t_grid_grid0_m01_m02_struct pars2[1]; /* pars for guesserror_VectorFuncP */
    double Cvec[3];
    int check, stat;

    pars2->grid = pars->grid;
    pars2->m01 = Getd("BNSdata_m01");
    pars2->m02 = Getd("BNSdata_m02");
    Cvec[1] = Getd("BNSdata_C1");
    stat = newton_linesrch_itsP(Cvec, 1, &check, m01_guesserror_VectorFuncP,
                                (void *) pars2, 30, pars->tol*0.05);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    Setd("BNSdata_C1", Cvec[1]);

    Cvec[1] = Getd("BNSdata_C2");
    if(n>=2)
      stat = newton_linesrch_itsP(Cvec, 1, &check, m02_guesserror_VectorFuncP,
                                  (void *) pars2, 30, pars->tol*0.05);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    Setd("BNSdata_C2", Cvec[1]);
    printf("q_at_xin1_2_VectorFuncP: C1=%.13g C2=%.13g\n",
           Getd("BNSdata_C1"), Getd("BNSdata_C2"));
  }
  /* compute q at xin1 */
  X=0.001;
  Y=1.0;
  Z=0.0;
  blist[0]=0;  blist[1]=1;
  bi = b_X_of_x_forgiven_YZ_inboxlist(pars->grid, blist,2, &X, pars->x1, Y,Z);
  q1 = BNS_compute_new_centered_q_atXYZ(pars->grid, bi, X,Y,Z);
  // printf("  xin1=%g is at bi=%d (X,Y,Z)=(%g,%g,%g)\n", pars->x1, bi,X,Y,Z);

  /* compute q at xin2 */
  X=0.001;
  Y=1.0;
  Z=0.0;
  blist[0]=3;  blist[1]=2;
  bi = b_X_of_x_forgiven_YZ_inboxlist(pars->grid, blist,2, &X, pars->x2, Y,Z);
  q2 = BNS_compute_new_centered_q_atXYZ(pars->grid, bi, X,Y,Z);
  // printf("  xin2=%g is at bi=%d (X,Y,Z)=(%g,%g,%g)\n", pars->x2, bi,X,Y,Z);

  printf("q_at_xin1_2_VectorFuncP: q(xin1)=%g q(xin2)=%g\n", q1, q2);
  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = q1;
  if(n>=2) fvec[2] = q2;
}

/* Adjust Omega and x_CM so that inner edges xin1/2 stay put */
int adjust_Omega_xCM_q_fix_xin(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int xind = Ind("x");
  double OmxCMvec[3];
  double dxin_m0[3];
  double dxin_00[3];
  double dxin_p0[3];
  double dxin_0m[3];
  double dxin_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int keepone   = Getv("BNSdata_adjust", "keep_one_xin");
  t_grid_x1_2_struct pars[1];
  int n1 = grid->box[0]->n1;
  int n2 = grid->box[0]->n2;
  int n3 = grid->box[0]->n3;

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_q_fix_xin: in BNSdata_solve step %d\n"
         "adjust_Omega_xCM_q_fix_xin: old Omega = %g  x_CM = %g\n",
         it, Omega, x_CM);

  /* set global vars */
  pars->grid = grid;
  pars->it   = it;
  pars->tol  = tol;
  /* pars->x1 = grid->box[0]->v[xind][Index(0,n2-1,0)];
     pars->x2 = grid->box[3]->v[xind][Index(0,n2-1,0)]; */
  pars->x1 = Getd("BNSdata_xin1");
  pars->x2 = Getd("BNSdata_xin2");
  printf("adjust_Omega_xCM_q_fix_xin: old xin1 = %g  xin2 = %g\n",
         pars->x1, pars->x2);
  prdivider(0);
  adjust_C1_C2_q_keep_restmasses(grid, it, tol*100.0);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars); */
    printf("adjust_Omega_xCM_q_fix_xin: dxin_m0[1]=%g dxin_m0[2]=%g\n",
           dxin_m0[1], dxin_m0[2]);
    /* printf("adjust_Omega_xCM_q_fix_xin: dxin_00[1]=%g dxin_00[2]=%g\n",
           dxin_00[1], dxin_00[2]); */
    printf("adjust_Omega_xCM_q_fix_xin: dxin_p0[1]=%g dxin_p0[2]=%g\n",
           dxin_p0[1], dxin_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxin_m0[1]*dxin_p0[1]<0.0) && (dxin_m0[2]*dxin_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxin_m0[1]*dxin_p0[1]<0.0) || (dxin_m0[2]*dxin_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xin to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxin_m0[1]*dxin_p0[1]>=0.0)
      {
        q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
        pars->x1 = grid->box[0]->v[xind][0];
      }
      if(dxin_m0[2]*dxin_p0[2]>=0.0)
      {
        q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
        pars->x2 = grid->box[3]->v[xind][0];
      }
      do_lnsrch = 1;
      printf("adjust_Omega_xCM_q_fix_xin: "
             "changed: old xin1 = %g  xin2 = %g\n",
             pars->x1, pars->x2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars); */
      printf("adjust_Omega_xCM_q_fix_xin: dxin_0m[1]=%g dxin_0m[2]=%g\n",
             dxin_m0[1], dxin_m0[2]);
      /* printf("adjust_Omega_xCM_q_fix_xin: dxin_00[1]=%g dxin_00[2]=%g\n",
             dxin_00[1], dxin_00[2]); */
      printf("adjust_Omega_xCM_q_fix_xin: dxin_0p[1]=%g dxin_0p[2]=%g\n",
             dxin_p0[1], dxin_p0[2]);
  
      if( (dxin_m0[1]*dxin_p0[1]<0.0) && (dxin_0m[2]*dxin_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxin_m0[2]*dxin_p0[2]<0.0) && (dxin_0m[1]*dxin_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_q_fix_xin: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xin1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, q_at_xin1_2_VectorFuncP,
                                (void *) pars, 1000, tol*0.1);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    q_at_xin1_2_VectorFuncP(2, OmxCMvec, dxin_00, (void *) pars);
    printf("adjust_Omega_xCM_q_fix_xin: dxin_00[1] = %g  dxin_00[2] = %g\n",
           dxin_00[1], dxin_00[2]);
  }
  printf("adjust_Omega_xCM_q_fix_xin: new Omega = %g  x_CM = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}


/* for newton_linesrch_its: compute error in xmaxs */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void xmaxs_error_VectorFunc(int n, double *vec, double *fvec)
{
  int bi1, bi2;
  double Xmax1,Ymax1, Xmax2,Ymax2;
  double xmax1, xmax2;
  tGrid *grid = xmaxs_error_VectorFunc__grid;

//  /* save old q in BNSdata_qgold */
//  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);
  printf("xmaxs_error_VectorFunc: Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,0.0);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,0.0);
  else
    xmax2 = Xmax2;
//  /* compute qmax1 and qmax2 */
//  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
//  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  printf("xmaxs_error_VectorFunc: Omega=%g x_CM=%g xmax1=%g xmax2=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"),
         xmax1, xmax2);  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = xmax1 - xmaxs_error_VectorFunc__xmax1;
  if(n>=2) fvec[2] = xmax2 - xmaxs_error_VectorFunc__xmax2;
}

/* Adjust Omega and x_CM so that point with max q, xmax1/2 stay put */
int adjust_Omega_xCM_keep_xmax(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  double OmxCMvec[3];
  double dxmax_m0[3];
  double dxmax_00[3];
  double dxmax_p0[3];
  double dxmax_0m[3];
  double dxmax_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int keepone   = Getv("BNSdata_adjust", "keep_one_xmax");

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_xmax: in BNSdata_solve step %d\n"
         "adjust_Omega_xCM_keep_xmax: old Omega = %g  x_CM = %g  tol = %g\n",
         it, Omega, x_CM, tol);

  /* set global vars */
  xmaxs_error_VectorFunc__grid  = grid;
  xmaxs_error_VectorFunc__xmax1 = Getd("BNSdata_xmax1");
  xmaxs_error_VectorFunc__xmax2 = Getd("BNSdata_xmax2");
  printf("adjust_Omega_xCM_keep_xmax: old xmax1 = %g  xmax2 = %g\n",
         xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_m0);
    OmxCMvec[1] = Omega + dOmega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_p0);
    OmxCMvec[1] = Omega;
    /* xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00); */
    printf("adjust_Omega_xCM_keep_xmax: dxmax_m0[1]=%g dxmax_m0[2]=%g\n",
           dxmax_m0[1], dxmax_m0[2]);
    /* printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1]=%g dxmax_00[2]=%g\n",
           dxmax_00[1], dxmax_00[2]); */
    printf("adjust_Omega_xCM_keep_xmax: dxmax_p0[1]=%g dxmax_p0[2]=%g\n",
           dxmax_p0[1], dxmax_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxmax_m0[1]*dxmax_p0[1]<0.0) && (dxmax_m0[2]*dxmax_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxmax_m0[1]*dxmax_p0[1]<0.0) || (dxmax_m0[2]*dxmax_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xmax to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxmax_m0[1]*dxmax_p0[1]>=0.0)
      {
        xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
        xmaxs_error_VectorFunc__xmax1 += dxmax_00[1];
      }
      if(dxmax_m0[2]*dxmax_p0[2]>=0.0)
      {
        xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
        xmaxs_error_VectorFunc__xmax2 += dxmax_00[2];
      }
      do_lnsrch = 1;
      printf("adjust_Omega_xCM_keep_xmax: "
             "changed: old xmax1 = %g  xmax2 = %g\n",
             xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_0m);
      OmxCMvec[2] = x_CM + dx_CM;
      xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_0p);
      OmxCMvec[2] = x_CM;
      /* xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00); */
      printf("adjust_Omega_xCM_keep_xmax: dxmax_0m[1]=%g dxmax_0m[2]=%g\n",
             dxmax_m0[1], dxmax_m0[2]);
      /* printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1]=%g dxmax_00[2]=%g\n",
             dxmax_00[1], dxmax_00[2]); */
      printf("adjust_Omega_xCM_keep_xmax: dxmax_0p[1]=%g dxmax_0p[2]=%g\n",
             dxmax_p0[1], dxmax_p0[2]);
  
      if( (dxmax_m0[1]*dxmax_p0[1]<0.0) && (dxmax_0m[2]*dxmax_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxmax_m0[2]*dxmax_p0[2]<0.0) && (dxmax_0m[1]*dxmax_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_xmax: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xmax1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_its(OmxCMvec, 2, &check, xmaxs_error_VectorFunc,
                               1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
    printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1] = %g  dxmax_00[2] = %g\n",
           dxmax_00[1], dxmax_00[2]);
  }
  printf("adjust_Omega_xCM_keep_xmax: new Omega = %g  x_CM = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* check where max are and reset BNSdata_xmax1/2 if they are too far off */
  if(Getv("BNSdata_adjust", "reset_xmax_if_problem"))
  {
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    xmaxs_error_VectorFunc__xmax1 = Getd("BNSdata_xmax1");
    xmaxs_error_VectorFunc__xmax2 = Getd("BNSdata_xmax2");
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
    if(fabs(dxmax_00[1])>5.0*tol) xmaxs_error_VectorFunc__xmax1 += dxmax_00[1];
    if(fabs(dxmax_00[2])>5.0*tol) xmaxs_error_VectorFunc__xmax2 += dxmax_00[2];
    Setd("BNSdata_xmax1", xmaxs_error_VectorFunc__xmax1);
    Setd("BNSdata_xmax2", xmaxs_error_VectorFunc__xmax2);
    /* print current maxs */
    printf("adjust_Omega_xCM_keep_xmax: resetting xmax1 = %g  xmax2 = %g\n",
           xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
    prdivider(0);
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_its: compute derivs of q in domain 0 and 3, or 4 and 5 */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void dqdx_at_Xmax1_2_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = dqdx_at_Xmax1_2_VectorFunc__grid;
  int bi1 = dqdx_at_Xmax1_2_VectorFunc__bi1;
  int bi2 = dqdx_at_Xmax1_2_VectorFunc__bi2;
  double Xmax1 = dqdx_at_Xmax1_2_VectorFunc__Xmax1;
  double Ymax1 = dqdx_at_Xmax1_2_VectorFunc__Ymax1;
  double Xmax2 = dqdx_at_Xmax1_2_VectorFunc__Xmax2;
  double Ymax2 = dqdx_at_Xmax1_2_VectorFunc__Ymax2;
  int b;

//  /* save old q in BNSdata_qgold */
//  varcopy(grid, Ind("BNSdata_qgold"), Ind("BNSdata_qg"));

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* get deriv dq of q in box bi1 and bi2 in BNSdata_temp1
     and dq's coeffs c in BNSdata_temp2 */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *q = box->v[Ind("BNSdata_q")];
    double *dq= box->v[Ind("BNSdata_temp1")];
    double *c = box->v[Ind("BNSdata_temp2")];
    spec_Deriv1(box, 1, q, dq);
    spec_Coeffs(box, dq, c);
  }

  /* find dq/dx at the locations Xmax1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("BNSdata_temp2")], Xmax1,Ymax1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("BNSdata_temp2")], Xmax2,Ymax2,0);

  printf("dqdx_at_Xmax1_2_VectorFunc: Omega=%.13g x_CM=%.13g\n"
         "  => dq/dx(xmax1)=%.13g",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"), fvec[1]);

  if(n>=2) printf("  dq/dx(xmax2)=%.13g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that point with max q, xmax1/2 stay put */
int adjust_Omega_xCM_keep_dqdxmax_eq_0(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int blist[6];
  double OmxCMvec[3];
  double dqdx_m0[3];
  double dqdx_00[3];
  double dqdx_p0[3];
  double dqdx_0m[3];
  double dqdx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double xmax1, xmax2;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  int pr_curxmax =Getv("BNSdata_adjust", "print_current_xmax");

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: in BNSdata_solve step %d\n"
         "adjust_Omega_xCM_keep_dqdxmax_eq_0: old Omega=%g x_CM=%g tol=%g\n",
         it, Omega, x_CM, tol);

  /* save BNSdata_xmax1/2 */
  xmax1 = Getd("BNSdata_xmax1");
  xmax2 = Getd("BNSdata_xmax2");

  /* set global vars */
  dqdx_at_Xmax1_2_VectorFunc__grid  = grid;
  /* for now we assume that the max are in box0/3 at Y=B=0
     or in box4/5 at Y=0 */
  dqdx_at_Xmax1_2_VectorFunc__Ymax1 = 0.0;
  blist[0]=0;  blist[1]=5;
  bi=b_X_of_x_forgiven_YZ_inboxlist(grid, blist, 2,
                                    &dqdx_at_Xmax1_2_VectorFunc__Xmax1,
                                    xmax1, 0.0,0.0);
  dqdx_at_Xmax1_2_VectorFunc__bi1 = bi;
  dqdx_at_Xmax1_2_VectorFunc__Ymax2 = 0.0;
  blist[0]=3;  blist[1]=4;
  bi=b_X_of_x_forgiven_YZ_inboxlist(grid, blist, 2,
                                    &dqdx_at_Xmax1_2_VectorFunc__Xmax2,
                                    xmax2, 0.0,0.0);
  dqdx_at_Xmax1_2_VectorFunc__bi2 = bi;
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: xmax1=%g xmax2=%g\n",
         xmax1, xmax2);
  if(dqdx_at_Xmax1_2_VectorFunc__bi1<0 || dqdx_at_Xmax1_2_VectorFunc__bi2<0)
  {
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: failed to find Xmax1 or "
           "Xmax2\n"
           " dqdx_at_Xmax1_2_VectorFunc__bi1=%d "
           " dqdx_at_Xmax1_2_VectorFunc__bi2=%d\n",
           dqdx_at_Xmax1_2_VectorFunc__bi1, dqdx_at_Xmax1_2_VectorFunc__bi2);
    return -1;
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: Xmax1=%g Xmax2=%g\n",
         dqdx_at_Xmax1_2_VectorFunc__Xmax1, dqdx_at_Xmax1_2_VectorFunc__Xmax2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_m0);
    OmxCMvec[1] = Omega + dOmega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_p0);
    OmxCMvec[1] = Omega;
    /* dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00); */
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_m0[1]=%g dqdx_m0[2]=%g\n",
           dqdx_m0[1], dqdx_m0[2]);
    /* printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]); */
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_p0[1]=%g dqdx_p0[2]=%g\n",
           dqdx_p0[1], dqdx_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_m0[2]*dqdx_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dqdx_m0[1]*dqdx_p0[1]<0.0) || (dqdx_m0[2]*dqdx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_0m);
      OmxCMvec[2] = x_CM + dx_CM;
      dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_0p);
      OmxCMvec[2] = x_CM;
      /* dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00); */
      printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_0m[1]=%g dqdx_0m[2]=%g\n",
             dqdx_m0[1], dqdx_m0[2]);
      /* printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
             dqdx_00[1], dqdx_00[2]); */
      printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_0p[1]=%g dqdx_0p[2]=%g\n",
             dqdx_p0[1], dqdx_p0[2]);
  
      if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_0m[2]*dqdx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dqdx_m0[2]*dqdx_p0[2]<0.0) && (dqdx_0m[1]*dqdx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xmax1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_its(OmxCMvec, 2, &check, dqdx_at_Xmax1_2_VectorFunc,
                               1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00);
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]);
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: new Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* find and print current xmax1/2 */
  if(pr_curxmax)
  {
    double dxmax[3];
    xmaxs_error_VectorFunc__grid  = grid;
    xmaxs_error_VectorFunc__xmax1 = Getd("BNSdata_xmax1");
    xmaxs_error_VectorFunc__xmax2 = Getd("BNSdata_xmax2");
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax);
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* compute q but with formula used for corotation case */
void compute_qcorot_with_corotation_formula(tGrid *grid, int qcor)
{
  char rotationstate1sav[1024], rotationstate2sav[1024];

  /* compute Func=(q we would have for corot) in var with ind qcor */
  snprintf(rotationstate1sav, 1024, "%s", Gets("BNSdata_rotationstate1"));
  snprintf(rotationstate2sav, 1024, "%s", Gets("BNSdata_rotationstate2"));
  Sets("BNSdata_rotationstate1", "corotation");
  Sets("BNSdata_rotationstate2", "corotation");
  BNS_compute_new_q(grid, qcor);
  Sets("BNSdata_rotationstate1", rotationstate1sav);
  Sets("BNSdata_rotationstate2", rotationstate2sav);
}

/* for newton_linesrch_itsP: compute derivs of Func in domain 0 and 3,
   or 4 and 5 */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void dFuncdx_at_Xfm1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int b, bi1, bi2;
  double Xfm1,Yfm1, Xfm2,Yfm2;
  t_grid_bXYZ1_bXYZ2_struct *pars;

  /* get pars */
  pars = (t_grid_bXYZ1_bXYZ2_struct *) p;
  grid = pars->grid;
  bi1 = pars->b1;
  bi2 = pars->b2;
  Xfm1 = pars->X1;
  Yfm1 = pars->Y1;
  Xfm2 = pars->X2;
  Yfm2 = pars->Y2;

  /* set BNSdata_Omega & BNSdata_x_CM */
  Setd("BNSdata_Omega", vec[1]);
  if(n>=2) Setd("BNSdata_x_CM",  vec[2]);

  /* compute Func=(q we would have for corot) in BNSdata_temp4 */
  compute_qcorot_with_corotation_formula(grid, Ind("BNSdata_temp4"));

  /* get deriv dFunc of Func in box bi1 and bi2 in BNSdata_temp1
     and dFunc's coeffs c in BNSdata_temp2 */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *Func = box->v[Ind("BNSdata_temp4")];
    double *dFunc= box->v[Ind("BNSdata_temp1")];
    double *c = box->v[Ind("BNSdata_temp2")];
    spec_Deriv1(box, 1, Func, dFunc);
    spec_Coeffs(box, dFunc, c);
  }

  /* find dFunc/dx at the locations Xfm1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("BNSdata_temp2")], Xfm1,Yfm1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("BNSdata_temp2")], Xfm2,Yfm2,0);

  printf("dFuncdx_at_Xfm1_2_VectorFuncP: Omega=%.13g x_CM=%.13g\n"
         "  => dFunc/dX(Xfm1)=%.13g",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"), fvec[1]);

  if(n>=2) printf("  dFunc/dX(Xfm2)=%.13g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that point with max Func, xfm1/2 stay put */
int adjust_Omega_xCM_keep_dFuncdxfm_eq_0(tGrid *grid, int it, double tol)
{
  int check, stat, bi, bi1,bi2, i;
  double OmxCMvec[3];
  double dqdx_m0[3];
  double dqdx_00[3];
  double dqdx_p0[3];
  double dqdx_0m[3];
  double dqdx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double Xfm1,Yfm1, Xfm2,Yfm2;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  t_grid_bXYZ1_bXYZ2_struct pars[1];

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: in BNSdata_solve step %d\n"
         "  old Omega=%g x_CM=%g tol=%g\n", it, Omega, x_CM, tol);

  /* set Xfm1 and Xfm2, i.e. find max in Func = BNSdata_qcorot */
  bi1=0;  bi2=3;
  find_Varmax_along_x_axis_usingBNSdata_temp123(grid, Ind("BNSdata_qcorot"),
                                                &bi1, &Xfm1, &Yfm1);
  find_Varmax_along_x_axis_usingBNSdata_temp123(grid, Ind("BNSdata_qcorot"),
                                                &bi2, &Xfm2, &Yfm2);
  /* set global vars */
  pars->grid  = grid;
  pars->b1 = bi1;
  pars->X1 = Xfm1;
  pars->Y1 = Yfm1;
  pars->Z1 = 0.0;
  pars->b2 = bi2;
  pars->X2 = Xfm2;
  pars->Y2 = Yfm2;
  pars->Z2 = 0.0;

  if(pars->b1<0 || pars->b2<0)
  {
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: failed to find Xfm1 or "
           "Xfm2\n"
           " pars->b1=%d "
           " pars->b2=%d\n",
           pars->b1, pars->b2);
    return -1;
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: Xfm1=%g Xfm2=%g\n",
         pars->X1, pars->X2);
  prdivider(0);
//  /* adjust C1/2 so that masses are correct */
//  adjust_C1_C2_q_keep_restmasses(grid, it, tol*100.0);
//  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars); */
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_m0[1]=%g dqdx_m0[2]=%g\n",
           dqdx_m0[1], dqdx_m0[2]);
    /* printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]); */
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_p0[1]=%g dqdx_p0[2]=%g\n",
           dqdx_p0[1], dqdx_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_m0[2]*dqdx_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dqdx_m0[1]*dqdx_p0[1]<0.0) || (dqdx_m0[2]*dqdx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars); */
      printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_0m[1]=%g dqdx_0m[2]=%g\n",
             dqdx_m0[1], dqdx_m0[2]);
      /* printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
             dqdx_00[1], dqdx_00[2]); */
      printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_0p[1]=%g dqdx_0p[2]=%g\n",
             dqdx_p0[1], dqdx_p0[2]);
  
      if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_0m[2]*dqdx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dqdx_m0[2]*dqdx_p0[2]<0.0) && (dqdx_0m[1]*dqdx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xfm1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, dFuncdx_at_Xfm1_2_VectorFuncP,
                               (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars);
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]);
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: new Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_itsP: compute derivs of Integrated Euler Eqn 
   in domain 0 and 3, or 4 and 5 */
/* if n=1 only BNSdata_Omega is adjusted
   if n=2 both BNSdata_Omega & BNSdata_x_CM are adjusted */
void dIntegEulerdx_at_X1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int b, bi1, bi2;
  double X1,Y1, X2,Y2;
  double Om, xcm;
  t_grid_bXYZ1_bXYZ2_struct *pars;

  /* get pars */
  pars = (t_grid_bXYZ1_bXYZ2_struct *) p;
  grid = pars->grid;
  bi1 = pars->b1;
  bi2 = pars->b2;
  X1 = pars->X1;
  Y1 = pars->Y1;
  X2 = pars->X2;
  Y2 = pars->Y2;

  /* set Omega & x_CM */
  Om  = vec[1];
  xcm = vec[2];

  /* compute lnIntegEuler in BNSdata_temp4 and 
     dlnIntegEulerdx/y/z in BNSdata_temp1/2/3 */
  BNS_set_dlnIntegEuler(grid, Ind("BNSdata_temp4"), 
                        Ind("BNSdata_temp1"), Om, xcm);

  /* put coeffs of BNSdata_temp1=dlnIntegEulerdx into BNSdata_temp3=c */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *dlnIE = box->v[Ind("BNSdata_temp1")];
    double *c     = box->v[Ind("BNSdata_temp3")];
    spec_Coeffs(box, dlnIE, c);
  }

  /* find dFunc/dx at the locations X1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("BNSdata_temp3")], X1,Y1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("BNSdata_temp3")], X2,Y2,0);

  printf("dIntegEulerdx_at_X1_2_VectorFuncP: Om=%.12g xcm=%.12g\n"
         " => dIntEuler/dX(X1)=%.12g", Om, xcm, fvec[1]);

  if(n>=2) printf(" dIntEuler/dX(X2)=%.12g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that force balance is maintained at xmax */
int adjust_Omega_xCM_forcebalance(tGrid *grid, int it, double tol)
{
  int check, stat, bi, bi1,bi2, i;
  double OmxCMvec[3];
  double dIEdxx_m0[3];
  double dIEdxx_00[3];
  double dIEdxx_p0[3];
  double dIEdxx_0m[3];
  double dIEdxx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double Xqm1,Yqm1, Xqm2,Yqm2;
  int do_lnsrch = Getv("BNSdata_adjust", "always");
  t_grid_bXYZ1_bXYZ2_struct pars[1];

  /* save old Omega, x_CM */
  Omega = Getd("BNSdata_Omega");
  x_CM  = Getd("BNSdata_x_CM");
  dOmega= Omega * Getd("BNSdata_dOmega_fac");
  dx_CM = Getd("BNSdata_b") * Getd("BNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_forcebalance: in BNSdata_solve step %d\n"
         "  old Omega=%g x_CM=%g tol=%g\n", it, Omega, x_CM, tol);

  /* set Xqm1 and Xqm2, i.e. find max in BNSdata_q */
  bi1=0;  bi2=3;
  find_Varmax_along_x_axis_usingBNSdata_temp123(grid, Ind("BNSdata_q"),
                                                &bi1, &Xqm1, &Yqm1);
  find_Varmax_along_x_axis_usingBNSdata_temp123(grid, Ind("BNSdata_q"),
                                                &bi2, &Xqm2, &Yqm2);
  /* set global vars */
  pars->grid  = grid;
  pars->b1 = bi1;
  pars->X1 = Xqm1;
  pars->Y1 = Yqm1;
  pars->Z1 = 0.0;
  pars->b2 = bi2;
  pars->X2 = Xqm2;
  pars->Y2 = Yqm2;
  pars->Z2 = 0.0;

  if(pars->b1<0 || pars->b2<0)
  {
    printf("adjust_Omega_xCM_forcebalance: failed to find Xqm1 or "
           "Xqm2\n"
           " pars->b1=%d "
           " pars->b2=%d\n",
           pars->b1, pars->b2);
    return -1;
  }
  printf("adjust_Omega_xCM_forcebalance:\n");
  printf("  box%d (X1,Y1)=(%g,%g) , box%d (X2,Y2)=(%g,%g)\n",
         pars->b1,pars->X1,pars->Y1, pars->b2,pars->X2,pars->Y2);
  prdivider(0);
//  /* adjust C1/2 so that masses are correct */
//  adjust_C1_C2_q_keep_restmasses(grid, it, tol*100.0);
//  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars); */
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_m0[1]=%g dIEdxx_m0[2]=%g\n",
           dIEdxx_m0[1], dIEdxx_m0[2]);
    /* printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
           dIEdxx_00[1], dIEdxx_00[2]); */
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_p0[1]=%g dIEdxx_p0[2]=%g\n",
           dIEdxx_p0[1], dIEdxx_p0[2]);
    prdivider(0);
    /* check if there is a zero, if so set do_lnsrch=1 */
    if( (dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) && (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0) )
      do_lnsrch = 1;
  }
  if( (do_lnsrch==0) &&
      ((dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) || (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars); */
      printf("adjust_Omega_xCM_forcebalance: dIEdxx_0m[1]=%g dIEdxx_0m[2]=%g\n",
             dIEdxx_m0[1], dIEdxx_m0[2]);
      /* printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
             dIEdxx_00[1], dIEdxx_00[2]); */
      printf("adjust_Omega_xCM_forcebalance: dIEdxx_0p[1]=%g dIEdxx_0p[2]=%g\n",
             dIEdxx_p0[1], dIEdxx_p0[2]);
  
      if( (dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) && (dIEdxx_0m[2]*dIEdxx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0) && (dIEdxx_0m[1]*dIEdxx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_forcebalance: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xfm1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, dIntegEulerdx_at_X1_2_VectorFuncP,
                               (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars);
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
           dIEdxx_00[1], dIEdxx_00[2]);
  }

  /* Nobody has touched the pars BNSdata_Omega and BNSdata_x_CM so far.
     Now set them to what we found */
  Setd("BNSdata_Omega", OmxCMvec[1]);
  Setd("BNSdata_x_CM",  OmxCMvec[2]);
  printf("adjust_Omega_xCM_forcebalance: new Omega=%g x_CM=%g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  return 0;
}

/* try to smoothen data after ell. solve by interpolation */
void smooth_BNSdata_by_Interpolation(tGrid *grid, int nsmooth)
{
  tGrid *grid2;
  int n;

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* initialize coords on grid2 */
  BNSgrid_init_Coords(grid2);

  /* do it a couple of times */
  for(n=1; n<=nsmooth; n++)
  {
    /* interpolate some vars from grid onto new grid2 */
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"));
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"));
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"));
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"));
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"));
    Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"));
    //  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"));
    //  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"));
    //  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"));
  
    /* copy grid2 back into grid */
    copy_grid(grid2, grid, 0);
  }
  /* free grid2 */
  free_grid(grid2);
}

/* adjust m01 and m02 to let them e.g. grow during iterations */
void adjust_BNSdata_m01_m02(void)
{
  double cm01 = Getd("BNSdata_m01");
  double cm02 = Getd("BNSdata_m02");
  double tm01 = Getd("BNSdata_desired_m01");
  double tm02 = Getd("BNSdata_desired_m02");
  double mCh1 = cm01*Getd("BNSdata_m0change");
  double mCh2 = cm02*Getd("BNSdata_m0change");

  /* change BNSdata_m01 and BNSdata_m02 by mCh1/2 to get closer to
     BNSdata_desired_m01 and BNSdata_desired_m02 */
  if(fabs(tm01-cm01) > mCh1) cm01 += mCh1*signum(tm01-cm01);
  else                       cm01  = tm01;
  if(fabs(tm02-cm02) > mCh2) cm02 += mCh2*signum(tm02-cm02);
  else                       cm02  = tm02;
  Setd("BNSdata_m01", cm01);
  Setd("BNSdata_m02", cm02);
  printf("adjust_BNSdata_m01_m02: BNSdata_m0change=%g mCh1=%g mCh2=%g\n",
         Getd("BNSdata_m0change"), mCh1, mCh2);
  printf(" => BNSdata_m01=%.13g BNSdata_m02=%.13g\n", cm01, cm02);
  printf("    BNSdata_desired_m01=%.13g BNSdata_desired_m01=%.13g\n", tm01, tm02);
  /* kappa adjustment ? */
  if(Getv("BNSdata_adjustkappa", "kappa^(-n/2)*m0_EQ_desired_kappa^(-n/2)*m0"))
  {
    double BNSdata_n = Getd("BNSdata_n");
    double kappa;
    kappa = Getd("BNSdata_desired_kappa") *
            pow((tm01+tm02)/(cm01+cm02), -2/BNSdata_n);
    Setd("BNSdata_kappa", kappa);
    printf("    BNSdata_kappa = %g\n", kappa);
    printf("    BNSdata_desired_kappa = %g\n", Getd("BNSdata_desired_kappa"));
  }
}

/* compute weighted average of current and old values,
   and return total error */
double average_current_and_old(double weight, tGrid *grid,
                               tVarList *vlFu, tVarList *vlu,
                               tVarList *vluDerivs, tVarList *vlJdu)
{
  double tm01 = Getd("BNSdata_m01");
  double tm02 = Getd("BNSdata_m02");
  double m01, m02;
  double normresnonlin, L2qdiff, dm01, dm02, m0err, error;

  /* if we iterate over the rest masses the true mass goals are different */
  if(Getv("BNSdata_iterate_m0", "yes"))
  {
    tm01 = Getd("BNSdata_desired_m01");
    tm02 = Getd("BNSdata_desired_m02");
  }

  /* reset new values from ell. solve as average between old and new.
     I.e. do: new = weight*new + (1-weight)*old  */
  varadd(grid, Ind("BNSdata_Psi"),   weight,Ind("BNSdata_Psi"),   (1.0-weight),Ind("BNSdata_Psiold"));
  varadd(grid, Ind("BNSdata_alphaP"),weight,Ind("BNSdata_alphaP"),(1.0-weight),Ind("BNSdata_alphaPold"));
  varadd(grid, Ind("BNSdata_Bx"),    weight,Ind("BNSdata_Bx"),    (1.0-weight),Ind("BNSdata_Boldx"));
  varadd(grid, Ind("BNSdata_By"),    weight,Ind("BNSdata_By"),    (1.0-weight),Ind("BNSdata_Boldy"));
  varadd(grid, Ind("BNSdata_Bz"),    weight,Ind("BNSdata_Bz"),    (1.0-weight),Ind("BNSdata_Boldz"));
  varadd(grid, Ind("BNSdata_Sigma"), weight,Ind("BNSdata_Sigma"), (1.0-weight),Ind("BNSdata_Sigmaold"));

  /* compute masses */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);
  
  /* compute error in masses */
  dm01 = (m01 - tm01)/(tm01+tm02);
  dm02 = (m02 - tm02)/(tm01+tm02);
  m0err = fabs(dm01) + fabs(dm02);
  
  /* evalute residual */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  normresnonlin = GridL2Norm(vlFu);
  printf("average_current_and_old:  weight=%g\n", weight);
  printf(" => m01=%.19g m02=%.19g\n", m01, m02);

  /* compute error in q. Do it like this: */
  /* compute error in new uncentered q, i.e. the diff between saved 
     BNSdata_qnocent and current uncentered q */
  BNS_compute_new_q(grid, Ind("BNSdata_temp2")); /* set temp2 = new unc. q */
  /* set temp1 = temp2 - qnocent = qnocent_new - qnocent_old */
  varadd(grid, Ind("BNSdata_temp1"),
               1,Ind("BNSdata_temp2"), -1,Ind("BNSdata_qnocent"));
  L2qdiff = varBoxL2Norm(grid->box[0], Ind("BNSdata_temp1")) +
            varBoxL2Norm(grid->box[3], Ind("BNSdata_temp1")) +
            varBoxL2Norm(grid->box[4], Ind("BNSdata_temp1")) +
            varBoxL2Norm(grid->box[5], Ind("BNSdata_temp1"));
  
  /* compute total error */
  error = normresnonlin + L2qdiff + m0err;
  printf(" => residual=%.4e L2qdiff=%.4e m0err=%.3e => error=%.4e\n",
         normresnonlin, L2qdiff, m0err, error);

  return error;
}

/* get total residual but without contribution in outer boxes */
double normresnonlin_without_BNSdata_Sigma_inbox12(tGrid *grid)
{
  int b;
  int iSig = Ind("BNSdata_Sigma_Err");

  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  for(b=1; b<=2; b++)
  {
    tBox *box = grid->box[b];
    int i;
    forallpoints(box, i)  box->v[iSig][i] = 0.0;
  }
  return GridL2Norm(vlFu);
}


/* Solve the Equations */
int BNSdata_solve(tGrid *grid)
{
  int    itmax        = Geti("BNSdata_itmax");
  double tol          = Getd("BNSdata_tol");
  double adjusttol    = max2(tol, Getd("BNSdata_adjust_mintol"));
  double esw          = Getd("BNSdata_esw");
  double esw1         = Getd("BNSdata_esw1");
  int    allow_esw1_it= Geti("BNSdata_allow_esw1_first_at");
  double Sigma_esw    = Getd("BNSdata_Sigma_esw");
  double Sigma_esw1   = Getd("BNSdata_Sigma_esw1");
  int    allow_Sigma_esw1_it = Geti("BNSdata_allow_Sigma_esw1_first_at");
  int    Newton_itmax = itmax;
  double NewtTolFac   = Getd("BNSdata_Newton_tolFac");
  double Newton_tol   = tol*NewtTolFac;
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  double realnormres = 1e300;
  double realnormres_old;
  double realSigmares, restres;
  int    itsSinceExtraSigma = 0;
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  tVarList *vldummy;
  int it;
  double dOmega = Getd("BNSdata_Omega")*0.1;
  double totalerr1, totalerr;
  char str[1000];

  /* choose linear solver */
  if(Getv("BNSdata_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("BNSdata_linSolver", "bicgstab_with_fd_UMFPACK_precon"))
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
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

// remove this later:
/*
Setd("BNSdata_C1", -0.832301);
Setd("BNSdata_C2", -0.832301);

Sets("BNSdata_center_new_q", "no");
grid->time  = -1000;
BNS_compute_new_centered_q(grid);
compute_new_q_and_adjust_domainshapes(grid, 0);
//compute_new_q_and_adjust_domainshapes(grid, 3);
write_grid(grid);

Sets("BNSdata_center_new_q", "yes");
grid->time  = -990;
BNS_compute_new_centered_q(grid);
compute_new_q_and_adjust_domainshapes(grid, 0);
//compute_new_q_and_adjust_domainshapes(grid, 3);
write_grid(grid);
exit(12);

grid->time  = -98;
compute_new_q_and_adjust_domainshapes(grid, 0);
compute_new_q_and_adjust_domainshapes(grid, 3);
write_grid(grid);

grid->time  = -97;
BNS_compute_new_centered_q(grid);
write_grid(grid);

exit(11);
*/

//grid->time  = -100;
//F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
//BNSdata_verify_solution(grid);
//printf("calling write_grid(grid)\n");
//write_grid(grid);
////
//BNS_compute_new_centered_q(grid);
//grid->time  = -99;
//F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
//BNSdata_verify_solution(grid);
//printf("calling write_grid(grid)\n");
//write_grid(grid);
////exit(11);
//return 1;

//Yo(1);CheckIfFinite(grid,  "BNSdata_q");

  /* start main iterations */
  prdivider(1);
  printf("BNSdata_solve: starting main iteration loop ...\n");

  /* choose initial Newton_tol */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  normresnonlin = GridL2Norm(vlFu);
  Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

  /* compute diagnostics like ham and mom */
  BNSdata_verify_solution(grid);
  
  /* output grid before any iterations are done */
  grid->time  = -itmax;
  write_grid(grid);
  BNSdata_analyze(grid);
  grid->time += 1.0;

  /* main iteration loop, do it until res is small enough */
  for(it=1; it <= itmax; it++)
  {
    int restart; 

    printf("BNSdata_solve step %d:\n", it);
    prdivider(1);

    /* do checkpointing (works only if checkpoint_restart_it<1) */
    restart = checkpoint(grid);

    /* update pars from file, and write new pars */
    if(parameterio_update_pars(grid) == 1 || restart == 1)
    {
      itmax        = Geti("BNSdata_itmax");
      tol          = Getd("BNSdata_tol");
      adjusttol    = max2(tol, Getd("BNSdata_adjust_mintol"));
      esw          = Getd("BNSdata_esw");
      esw1         = Getd("BNSdata_esw1");
      allow_esw1_it= Geti("BNSdata_allow_esw1_first_at");
      Sigma_esw    = Getd("BNSdata_Sigma_esw");
      Sigma_esw1   = Getd("BNSdata_Sigma_esw1");
      allow_Sigma_esw1_it = Geti("BNSdata_allow_Sigma_esw1_first_at");
      Newton_itmax = itmax;
      NewtTolFac   = Getd("BNSdata_Newton_tolFac");
      linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
      linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
      linSolver_tol    = Getd("BNSdata_linSolver_tol");
      /* reset Newton_tol */
      F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
      normresnonlin = GridL2Norm(vlFu);
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
    }
    parameterio_write_current_pars(grid);

    /* If we read a checkpoint that had finished the solve, i.e.
       that has $time = 0, do not solve again! */
    if(grid->time==0.0 && restart == 1) break;

    /* filter vars if desired */
    BNSdata_filter_Vars(grid, it);

    /* what to do with q at A=0 and q<0 */
    if(Getv("BNSdata_set_negative_q", "zero"))
      set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);
    if(Getv("BNSdata_set_q_atA0", "zero"))
      set_Var_to_Val_atA0(grid, Ind("BNSdata_q"), 0.0);

    /* center q first, then solve and adjust C1/2, Omega, xCM. */
    if(Getv("BNSdata_center_new_q_timebin", "before_ell_solve"))
      BNSdata_center_q_if_desired(grid, it);

    /* center fields around each star and also keep inner edge if desired */
    BNSdata_center_fields_if_desired(grid, it);
    BNSdata_keep_xin_if_desired(grid, it);

    /* save old values before ell. solve */
    varcopy(grid, Ind("BNSdata_Psiold"),    Ind("BNSdata_Psi"));
    varcopy(grid, Ind("BNSdata_alphaPold"), Ind("BNSdata_alphaP"));
    varcopy(grid, Ind("BNSdata_Boldx"),     Ind("BNSdata_Bx"));
    varcopy(grid, Ind("BNSdata_Boldy"),     Ind("BNSdata_By"));
    varcopy(grid, Ind("BNSdata_Boldz"),     Ind("BNSdata_Bz"));
    varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
    varcopy(grid, Ind("BNSdata_qgold"),     Ind("BNSdata_qg"));

    /* set BNSdata_qcorot before ell solve */
    compute_qcorot_with_corotation_formula(grid, Ind("BNSdata_qcorot"));
    
    /* set BNSdata_qnocent before ell solve */
    BNS_compute_new_q(grid, Ind("BNSdata_qnocent"));

    /* set wB before we solve */
    BNS_set_wB(grid, 1, Getd("BNSdata_actual_xmax1"),0.0,0.0); 
    BNS_set_wB(grid, 2, Getd("BNSdata_actual_xmax2"),0.0,0.0); 


    /* check if we do another ell. solve for BNSdata_Sigma */
    realnormres_old = realnormres; /* save realnormres */
    realnormres = normresnonlin_without_BNSdata_Sigma_inbox12(grid);
    itsSinceExtraSigma++;
    prdivider(1);
    printf("BNSdata_solve step %d: itsSinceExtraSigma=%d\n",
           it, itsSinceExtraSigma);
    printf(" realnormres=%g   realnormres_old=%g\n",
           realnormres, realnormres_old);
    if( (realnormres_old <= realnormres*Getd("BNSdata_extraSigmaSolve_fac")) &&
        (itsSinceExtraSigma >= Geti("BNSdata_extraSigmaSolve_every")) )
    {
      /* solve BNSdata_Sigma completely in outer boxes, unless we do it later */
      if(Getv("BNSdata_SmoothSigma", "no"))
      {
        printf("Setting BNSdata_Sigma outside the stars only...\n");
        Sets("BNSdata_KeepInnerSigma", "yes");
        BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
        Sets("BNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        /* reset Sigmaold to take into account new Sigma in outer boxes */
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
      }
      /* reset Newton_tol, use only other vars */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("BNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* solve the ell. eqn for Sigma alone */
      BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
             &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
      totalerr1 = average_current_and_old(Sigma_esw, 
                                          grid,vlFu,vlu,vluDerivs, vlJdu);
      if(Sigma_esw<1.0 && it>=allow_Sigma_esw1_it && allow_Sigma_esw1_it>=0)
      {
        /* complete step */
        totalerr = average_current_and_old(Sigma_esw1/Sigma_esw,
                                           grid,vlFu,vlu,vluDerivs,vlJdu);
        /* but go back to Sigma_esw if totalerr is larger */
        if(totalerr>totalerr1)
          totalerr = average_current_and_old(Sigma_esw/Sigma_esw1, 
                                             grid,vlFu,vlu,vluDerivs,vlJdu);
      }
      /* try to smooth BNSdata_Sigma near the boundary */
      smooth_BNSdata_Sigma_NearBoundary(grid, 2, tol, linear_solver);

      /* reset Newton_tol, use error norm of all vars in vlu */
      F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
      normresnonlin = GridL2Norm(vlFu);
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* reset Sigmaold so that Sigma does not change when we average later */
      varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));

      /* make sure we do not enter this block in the next iteration */
      realnormres = 1e300;
      itsSinceExtraSigma = 0;
    }

    /* How we solve the coupled ell. eqns */
    if(Getv("BNSdata_EllSolver_method", "allatonce"))
    { 
      /* solve BNSdata_Sigma completely in outer boxes at first iteration */
      if(grid->time == 1.0-itmax)
      {
        printf("Setting BNSdata_Sigma outside the stars only (using UMFPACK)...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("BNSdata_KeepInnerSigma", "yes");
        // BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
        //       &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
        BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, UMFPACK_solve_wrapper, 1, "BNSdata_Sigma");
        Sets("BNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* solve the coupled ell. eqns all together */
      vldummy = vlJdu;
      Newton(F_BNSdata, J_BNSdata, vlu, vlFu, vluDerivs, vldummy,
             Newton_itmax, Newton_tol, &normresnonlin, 1,
             linear_solver, Preconditioner_I, vldu, vlJdu, vlduDerivs, vlu,
             linSolver_itmax, linSolver_tolFac, linSolver_tol);
    }
    else if(Getv("BNSdata_EllSolver_method", "BNS_ordered_Eqn_Iterator"))
    {
      // /* reset Newton_tol, so that we always solve for Sigma */
      // normresnonlin = GridL2Norm_of_vars_in_string(grid, "BNSdata_Sigma");
      // Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* solve completely in outer boxes at first iteration */
      if(grid->time == 1.0-itmax)
      {
        printf("Setting BNSdata_Sigma outside the stars only...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("BNSdata_KeepInnerSigma", "yes");
        BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
        Sets("BNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* do we want to solve for Sigma? */
      realSigmares =
       GridL2Norm_of_vars_in_string_withZeroErr_inbox12(grid, "BNSdata_Sigma");
      restres = GridL2Norm_of_vars_in_string(grid,  
                                      Gets("BNSdata_CTS_Eqs_Iteration_order"));
      printf(" realSigmares=%g  restres=%g\n", realSigmares, restres);
      if( realSigmares >= restres * Getd("BNSdata_SigmaSolve_tolFac") )
      {
        /* solve the ell. eqn for Sigma alone */
        BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
        totalerr1 = average_current_and_old(Sigma_esw, 
                                            grid,vlFu,vlu,vluDerivs, vlJdu);
        if(Sigma_esw<1.0 && it>=allow_Sigma_esw1_it && allow_Sigma_esw1_it>=0)
        {
          /* complete step */
          totalerr = average_current_and_old(Sigma_esw1/Sigma_esw,
                                             grid,vlFu,vlu,vluDerivs,vlJdu);
          /* but go back to Sigma_esw if totalerr is larger */
          if(totalerr>totalerr1)
            totalerr = average_current_and_old(Sigma_esw/Sigma_esw1, 
                                               grid,vlFu,vlu,vluDerivs,vlJdu);
        }
        /* try to smooth BNSdata_Sigma near the boundary, recomp Err */
        smooth_BNSdata_Sigma_NearBoundary(grid, 2, tol, linear_solver);
        F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);

        /* reset Sigmaold so that Sigma does not change when we average later */
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
      }

      /* reset Newton_tol */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("BNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* now solve the coupled CTS ell. eqns one after an other */
      BNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method",
                 "BNS_ordered_Eqn_Iterator_SphereSigma"))
    {
      /* set VolAvSigma1/2 to zero for the BNSdata_Sigma solve in this block */
      Sets("BNSdata_desired_VolAvSigma1", "0");
      Sets("BNSdata_desired_VolAvSigma2", "0");
      printf(" setting: BNSdata_desired_VolAvSigma1/2 to zero\n");
  
      // /* reset Newton_tol, so that we always solve for Sigma */
      // normresnonlin = GridL2Norm_of_vars_in_string(grid, "BNSdata_Sigma");
      // Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* solve completely in outer boxes at first iteration */
      if(grid->time == 1.0-itmax)
      {
        printf("Setting BNSdata_Sigma outside the stars only...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("BNSdata_KeepInnerSigma", "yes");
        solve_BNSdata_Sigma_WithSphereBoundary(grid, Newton_itmax, Newton_tol,
                                               &normresnonlin, linear_solver);
        Sets("BNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* solve the ell. eqn for Sigma alone */
      solve_BNSdata_Sigma_WithSphereBoundary(grid, Newton_itmax, Newton_tol,
                                             &normresnonlin, linear_solver);
      /* reset VolAvSigma1/2 after solve */
      set_BNSdata_desired_VolAvSigma12_toMinBCerr(grid, Ind("BNSdata_Sigmaold"));
      totalerr1 = average_current_and_old(Sigma_esw, 
                                          grid,vlFu,vlu,vluDerivs, vlJdu);
      if(Sigma_esw<1.0 && it>=allow_Sigma_esw1_it && allow_Sigma_esw1_it>=0)
      {
        /* complete step */
        totalerr = average_current_and_old(Sigma_esw1/Sigma_esw,
                                           grid,vlFu,vlu,vluDerivs,vlJdu);
        /* but go back to Sigma_esw if totalerr is larger */
        if(totalerr>totalerr1)
          totalerr = average_current_and_old(Sigma_esw/Sigma_esw1, 
                                             grid,vlFu,vlu,vluDerivs,vlJdu);
      }
      /* reset Sigmaold so that Sigma does not change when we average later */
      varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
      /* reset VolAvSigma1/2 again */
      set_BNSdata_desired_VolAvSigma12_toMinBCerr(grid, Ind("BNSdata_Sigma"));

      /* reset Newton_tol */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("BNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* now solve the coupled CTS ell. eqns one after an other */
      BNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "BNS_ordered_Var_Eqn_Iterator"))
    {
      /* solve completely in outer boxes at first iteration */
      if(grid->time == 1.0-itmax)
      {
        printf("Setting BNSdata_Sigma outside the stars only...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("BNSdata_KeepInnerSigma", "yes");
        BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "BNSdata_Sigma");
        Sets("BNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
        /* reset Sigmaold so that Sigma does not change when we average later */
        varcopy(grid, Ind("BNSdata_Sigmaold"),  Ind("BNSdata_Sigma"));
      }
      /* reset Newton_tol */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("BNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* now solve ell. eqns for vars in BNSdata_CTS_Eqs_Iteration_order
         one after an other */
      BNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "BNS_ordered_Eqn_Iterator_old"))
    { /* solve the coupled ell. eqns one after an other */
      if(!Getv("BNSdata_CTS_Eqs_Iteration_order", "BNSdata_Sigma"))
        errorexit("BNSdata_Sigma is missing in BNSdata_CTS_Eqs_Iteration_order");
      BNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "BNS_Eqn_Iterator"))
    { /* solve the coupled ell. eqns one after an other */
      BNS_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                       linear_solver, 1);
//      BNS_Eqn_Iterator(grid, Newton_itmax, tol*0.0001, &normresnonlin,
//                       linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "sequence1"))
    { 
      BNS_Eqn_sequence1(grid, Newton_itmax, tol*0.0001, &normresnonlin,
                       linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "sequence2"))
    {
//      BNS_Eqn_sequence2(grid, Newton_itmax, Newton_tol, &normresnonlin,
//                       linear_solver, 1);
      BNS_Eqn_sequence2(grid, Newton_itmax, tol*0.0001, &normresnonlin,
                       linear_solver, 1);
    }
    else if(Getv("BNSdata_EllSolver_method", "sequence3"))
    {
      BNS_Eqn_sequence3(grid, Newton_itmax, tol*0.0001, &normresnonlin,
                       linear_solver, 1);
    }
    else
      errorexit("BNSdata_solve: unknown BNSdata_EllSolver_method");

    /* if we smoothen data by Interpolation*/
    if(Geti("BNSdata_Interpolation_smooths")>0)
    {
      int n = Geti("BNSdata_Interpolation_smooths");
      printf("calling smooth_BNSdata_by_Interpolation %d times...\n", n);
      fflush(stdout);
      smooth_BNSdata_by_Interpolation(grid, n);
    }

    /* enforce uniqueness on axis:
       set vars equal to val at phi=0 for all phi>0 */
    if(Getv("BNSdata_uniqueness_on_axis", "yes"))
      BNS_enforce_uniqueness_on_axis(vlu);

    /* if we iterate rest masses */
    if(Getv("BNSdata_iterate_m0", "yes")) adjust_BNSdata_m01_m02();

    /* reset new values from ell. solve as average between old and new.
       I.e. do: new = esw*new + (1-esw)*old  */
    totalerr1 = average_current_and_old(esw, grid,vlFu,vlu,vluDerivs, vlJdu);
    if(esw<1.0 && it>=allow_esw1_it && allow_esw1_it>=0)
    {
      /* complete step */
      totalerr = average_current_and_old(esw1/esw, grid,vlFu,vlu,vluDerivs,vlJdu);
      /* but go back to esw if totalerr is larger */
      if(totalerr>totalerr1)
        totalerr = average_current_and_old(esw/esw1, grid,vlFu,vlu,vluDerivs,vlJdu);
    }
    else
      totalerr = totalerr1;

    /* compute diagnostics like ham and mom */
    BNSdata_verify_solution(grid);

    /* break if total error is small enough */
    if(totalerr<tol &&
       Getv("BNSdata_break_if_err_below_tol","after_ell_solve")) break;

    /* write after elliptic solve, but before adjusting q */
    grid->time -= 0.5;
    write_grid(grid);
    grid->time += 0.5;

if(0)
{
    /* use equal mass symmetry and write grid */
    BNSgrid_set_allVars_onLeft_equalmasses(grid);
    F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
    BNSdata_verify_solution(grid);
    grid->time -= 0.4;
    write_grid(grid);
    grid->time += 0.4;
printf(" after symmetry:  => m01=%.19g m02=%.19g\n",
       GetInnerRestMass(grid, 0), GetInnerRestMass(grid, 3));
}
if(0)
{
    /* adjust C1/2 according to Pedro's algorithm. This yields
       a new q as well.  Note: THIS FAILS!!! */
    adjust_C1_C2_Omega_q_Pedro(grid, it, tol);
}
if(0)
{
    /* adjust C1/2, Omega, xCM according to BGM */
    adjust_C1_C2_Omega_q_BGM(grid, it, tol);
}
if(0) /* HYBRID <-- not working... */
{
  int bi1, bi2;
  double Xmax1,Ymax1, Xmax2,Ymax2;
  double qmax1, qmax2;

  /* adjust such that masses are initial masses */
  adjust_C1_C2_Omega_q_Pedro(grid, it, tol);

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
//  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
//  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  /* make both max densities the same */
  Setd("BNSdata_qmax1", qmax1);
  Setd("BNSdata_qmax2", qmax1);          
  adjust_C1_C2_Omega_q_BGM(grid, it, tol);
}
if(0) /* not working */
{
    /* adjust C1/2, Omega, xCM according to WT */
    adjust_C1_C2_Omega_xCM_q_WT(grid, it, tol, &dOmega);
}
if(0) /* not working */
{
    /* adjust C1/2, Omega, xCM according to WT with L2 */
    adjust_C1_C2_Omega_xCM_q_WT_L2(grid, it, tol, &dOmega);
}


    /* reset BNSdata_qmax1/2, BNSdata_xmax1/2 pars if the iteration 
       # it is contained in the list BNSdata_reset_qmax_xmax_pars_at */
    snprintf(str, 999, "%d", it);
    if( Getv("BNSdata_reset_qmax_xmax_pars_at", str) )
      find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars(grid);
      
    /* choose how we adjust C1/2, Omega, xCM: */
    if( it>=Geti("BNSdata_adjust_first_at") &&
        Geti("BNSdata_adjust_first_at")>=0 )
    {
      if(Getv("BNSdata_adjust", "WT_L2_method"))
        adjust_C1_C2_Omega_xCM_q_WT_L2(grid, it, adjusttol, &dOmega);
      else if(Getv("BNSdata_adjust", "keep_xout")) /* old keep_xout */
        adjust_C1_C2_Omega_xCM_q_keep_xout(grid, it, adjusttol);
      else if(Getv("BNSdata_adjust", "keep_xin")) /* old keep_xin */
        adjust_C1_C2_Omega_xCM_q_keep_xin(grid, it, adjusttol);
      else if(Getv("BNSdata_adjust", "fix_xout")) /* alternative keep_xout */
      { /* keep xout1/2 in place */
        adjust_Omega_xCM_q_fix_xout(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "fix_xin")) /* alternative keep_xin */
      { /* keep xin1/2 in place */
        adjust_Omega_xCM_q_fix_xin(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "findandkeep_xmax")) /* old keep_xmax */
      { /* keep xmax1/2 in place */
        adjust_Omega_xCM_keep_xmax(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "keep_xmax"))
      { /* keep xmax1/2 in place, by keeping the pos. of dq/dx=0 in place */
        adjust_Omega_xCM_keep_dqdxmax_eq_0(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "keep_xfm"))
      { /* keep max1/2 of a funtion in place, by keeping the pos. of dFunc/dx=0 in place */
        adjust_Omega_xCM_keep_dFuncdxfm_eq_0(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "forcebalance"))
      { /* keep max1/2 of a funtion in place, by keeping the pos. of dFunc/dx=0 in place */
        adjust_Omega_xCM_forcebalance(grid, it, adjusttol);
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
      }
      else if(Getv("BNSdata_adjust", "min_qchange"))
        adjust_C1_C2_Omega_xCM_q_min_qchange(grid, it, adjusttol, &dOmega);
      else /* adjust C1/2, q while keeping restmasses, Omega and xCM */
        adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol);
    }
    else
    { /* adjust C1/2, q while keeping restmasses, Omega and xCM */
      adjust_C1_C2_q_keep_restmasses(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_restmasses multiplies its tol with 0.01 */
    }

    /* enforce uniqueness on axis:
       set vars equal to val at phi=0 for all phi>0 */
    if(Getv("BNSdata_uniqueness_on_axis", "yes"))
      BNS_enforce_uniqueness_on_axis(vlu);

    /* compute actual max pos of q and center q if needed and 
       if BNSdata_center_new_q is not "no": */
    /* set actual positions of maxima */
    set_BNSdata_actual_xyzmax_pars(grid);
    if(Getv("BNSdata_center_new_q_timebin", "after_adjusting_Omega_xCM"))
      BNSdata_center_q_if_desired(grid, it);

    /* Set VolAvSigma1/2 so that the next ell. solves all try to 
       achieve a certain Volume Average for BNSdata_Sigma. This also
       results in residuals that do not take into account the arbitrary 
       constant that can be added to BNSdata_Sigma. */
    set_BNSdata_desired_VolAvSigma12_toMinBCerr(grid, Ind("BNSdata_Sigma"));

    /* compute diagnostics like ham and mom */
    BNSdata_verify_solution(grid);

    /* evalute residual and break if it is small enough */
    F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
    normresnonlin = GridL2Norm(vlFu);
    printf("BNSdata_solve step %d: residual = %e\n", it, normresnonlin);
    totalerr1 = normresnonlin_without_BNSdata_Sigma_inbox12(grid);
    printf(" with Sigma_Err=0 in box1/2: real residual = %e\n", totalerr1);
    prdivider(1);  fflush(stdout);
    if((normresnonlin<tol || totalerr1<tol) && 
       Getv("BNSdata_break_if_err_below_tol","at_iterationend")) break;

    /* set new tol for Newton */
    Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

    /* write current iteration if we are not done yet and increase counters */
    if(it<=itmax)
    {
      write_grid(grid);
      BNSdata_analyze(grid);
    }
    grid->time += 1.0;
  }
  if(it>itmax)
    printf("BNSdata_solve warning: *** Too many steps! ***\n");

  /* do we want to do a final ell. solve for some vars? */
  if( strlen(Gets("BNSdata_FinalEllSolveVars"))>0 )
  {
    double time = grid->time;
    /* write grid once more (at time=-0.1), since we have not written yet */
    grid->time = -0.1;
    printf("Final elliptic solve at grid->time = -0.1\n");
    write_grid(grid);
    BNSdata_analyze(grid);
    grid->time = time;

    BNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, tol,
                                        &normresnonlin, linear_solver, 1,
                                        Gets("BNSdata_FinalEllSolveVars"));
    BNSdata_verify_solution(grid);
  }

  /* now we have intial data, set time=0 */
  grid->time = 0.0;

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
int BNSdata_verify_solution(tGrid *grid)
{
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

  /* compute constraint and other things that should be zero */
  computeADMconstraints(grid);

  return 0;
}

/* compute absolute error in ANALYSIS */
int BNSdata_analyze(tGrid *grid)
{
  double BNSdata_b = Getd("BNSdata_b");
  double n = Getd("BNSdata_n");
  double Gamma = 1.0 + 1.0/n;
  double kappa = Getd("BNSdata_kappa");
  double Omega = Getd("BNSdata_Omega");
  double x_CM  = Getd("BNSdata_x_CM");
  double xout1, xin1, xin2, xout2;
  double M_ADM, J_ADM, Jx_ADM,Jy_ADM,Jz_ADM, m01, m02;
  double Jx_ADM1,Jy_ADM1,Jz_ADM1, Jx_ADM2,Jy_ADM2,Jz_ADM2;
  double *M_ADM_ofA;
  double M_ADM_in0, M_ADM_in3, M_ADM_in1, M_ADM_in2; 
  int iX = Ind("X");
  int ix = Ind("x");
  int isigma = Ind("Coordinates_AnsorgNS_sigma_pm");
  int b1n1 = grid->box[1]->n1;
  int b1n2 = grid->box[1]->n2;
  int b2n1 = grid->box[2]->n1;
  int b2n2 = grid->box[2]->n2;
  double sigp_00 = grid->box[1]->v[isigma][0]; /* sigma_p at B=phi=0 */
  double sigp_10 = grid->box[1]->v[isigma][Ind_n1n2(0,b1n2-1,0, b1n1,b1n2)]; /* sigma_p at B=1, phi=0 */
  double sigm_00 = grid->box[2]->v[isigma][0]; /* sigma_m at B=phi=0 */
  double sigm_10 = grid->box[2]->v[isigma][Ind_n1n2(0,b2n2-1,0, b2n1,b2n2)]; /* sigma_m at B=1, phi=0 */
  int itemp1 = Ind("BNSdata_temp1");
  int itemp2 = Ind("BNSdata_temp2");
  int itemp3 = Ind("BNSdata_temp3");
  double *temp1 = grid->box[1]->v[itemp1];
  double xmax1, qmax1, xmax2, qmax2;
  int bi1, bi2;
  double Xmax1,Ymax1, Xmax2,Ymax2;
  double Zmax1,Zmax2;
  double glob_xmax1, glob_ymax1, glob_zmax1;
  double glob_xmax2, glob_ymax2, glob_zmax2;
  double global_qmax1, global_qmax2;
  double TOV_rf_surf1, TOV_m1, TOV_Phic1, TOV_Psic1, TOV_m01;  /* for TOV */
  double TOV_rf_surf2, TOV_m2, TOV_Phic2, TOV_Psic2, TOV_m02;  /* for TOV */
  double TOV_r_surf1, TOV_r_surf2, TOV_Psis1, TOV_Psis2;
  FILE *fp;
  char *outdir = Gets("outdir");
  char *name = "BNSdata_properties.txt";
  char *filename;
  int filenamelen;
  int i;
  
  printf("BNSdata_analyze: computing properties of BNS data\n");

  /* get inner and outer edges of both stars */
  xout1= grid->box[0]->x_of_X[1]((void *) grid->box[0], -1, 0.0,0.0,0.0);
  xin1 = grid->box[0]->x_of_X[1]((void *) grid->box[0], -1, 0.0,1.0,0.0);
  xin2 = grid->box[3]->x_of_X[1]((void *) grid->box[3], -1, 0.0,1.0,0.0);
  xout2= grid->box[3]->x_of_X[1]((void *) grid->box[3], -1, 0.0,0.0,0.0);

  /* compute rest masses m01, m02 */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);

  /* compute ADM ang. mom. */
  BNS_set_J_ADM_VolInt_integrand(grid, itemp1,itemp2,itemp3);
  Jx_ADM1 = InnerVolumeIntegral(grid, 0, itemp1);
  Jy_ADM1 = InnerVolumeIntegral(grid, 0, itemp2);
  Jz_ADM1 = InnerVolumeIntegral(grid, 0, itemp3);
  Jx_ADM2 = InnerVolumeIntegral(grid, 3, itemp1);
  Jy_ADM2 = InnerVolumeIntegral(grid, 3, itemp2);
  Jz_ADM2 = InnerVolumeIntegral(grid, 3, itemp3);
  Jx_ADM = Jx_ADM1 + Jx_ADM2;
  Jy_ADM = Jy_ADM1 + Jy_ADM2;
  Jz_ADM = Jz_ADM1 + Jz_ADM2;
  J_ADM = sqrt(Jx_ADM*Jx_ADM + Jy_ADM*Jy_ADM + Jz_ADM*Jz_ADM);
  printf("Jx_ADM1=%.9g Jy_ADM1=%.9g Jz_ADM1=%.9g "
         "Jx_ADM2=%.9g Jy_ADM2=%.9g Jz_ADM2=%.9g\n",
         Jx_ADM1,Jy_ADM1,Jz_ADM1, Jx_ADM2,Jy_ADM2,Jz_ADM2);

  /* compute ADM mass from second deriv of Psi */
  M_ADM = ADMmass_fromPsi_inbox1_at_A1B0(grid, itemp1);
  M_ADM_ofA = dmalloc(b1n1);
  for(i=0; i<b1n1; i++) M_ADM_ofA[i] = grid->box[1]->v[itemp1][i];

  /* compute ADM mass from Volume int */
  BNS_set_M_ADM_VolInt_integrand(grid, itemp1);
  M_ADM_in0 = InnerVolumeIntegral(grid, 0, itemp1);
  M_ADM_in3 = InnerVolumeIntegral(grid, 3, itemp1);
  M_ADM_in1 = VolumeIntegral_inBNSgridBox(grid, 1, itemp1);
  M_ADM_in2 = VolumeIntegral_inBNSgridBox(grid, 2, itemp1);
  M_ADM = M_ADM_in0+M_ADM_in3+M_ADM_in1+M_ADM_in2;

  printf("M_ADM_in0=%.9g M_ADM_in3=%.9g M_ADM_in1=%.9g M_ADM_in2=%.9g\n",
         M_ADM_in0,M_ADM_in3, M_ADM_in1,M_ADM_in2);
  printf("ADM quantities: M_ADM = %.16g  J_ADM = %.16g\n", M_ADM, J_ADM);

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);
  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,0.0);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,0.0);
  else
    xmax2 = Xmax2;
  /* set qmax1/2 */
  Setd("BNSdata_qmax1", qmax1);
  Setd("BNSdata_qmax2", qmax2);
  printf("BNSdata_analyze: BNSdata_qmax1 = %.10g  BNSdata_qmax2=%.10g\n"
         "  at:                    xmax1 = %.10g          xmax2=%.10g\n",
         Getd("BNSdata_qmax1"), Getd("BNSdata_qmax2"), xmax1, xmax2);
  /* set cart positions xmax1/2 of qmax1/2 ? */
  if(Getv("BNSdata_analyze_xmax", "set_BNSdata_xmax"))
  {
    Setd("BNSdata_xmax1", xmax1);
    Setd("BNSdata_xmax2", xmax2);
    printf("  setting:       BNSdata_xmax1 = %.10g  BNSdata_xmax2=%.10g\n",
           Getd("BNSdata_xmax1"), Getd("BNSdata_xmax2"));
  }
  else
  {
    printf("  keeping:       BNSdata_xmax1 = %.10g  BNSdata_xmax2=%.10g\n",
           Getd("BNSdata_xmax1"), Getd("BNSdata_xmax2"));
  }

  /* find global max of q in NS1/2 */
  Zmax1=Zmax2=0.0;
  global_qmax1 = BNSdata_find_position_of_qmax(grid, &bi1, &Xmax1, &Ymax1, &Zmax1);
  global_qmax2 = BNSdata_find_position_of_qmax(grid, &bi2, &Xmax2, &Ymax2, &Zmax2);
  if(grid->box[bi1]->x_of_X[1] != NULL)
  {
    glob_xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
    glob_ymax1 = grid->box[bi1]->x_of_X[2]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
    glob_zmax1 = grid->box[bi1]->x_of_X[3]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
  }
  else
  {
    glob_xmax1 = Xmax1;
    glob_ymax1 = Ymax1;
    glob_zmax1 = Zmax1;
  }
  if(grid->box[bi2]->x_of_X[1] != NULL)
  {
    glob_xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
    glob_ymax2 = grid->box[bi2]->x_of_X[2]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
    glob_zmax2 = grid->box[bi2]->x_of_X[3]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
  }
  else
  {
    glob_xmax2 = Xmax2;
    glob_ymax2 = Ymax2;
    glob_zmax2 = Zmax2;
  }
  printf("BNSdata_analyze: global qmax1=%g in box%d\n"
         "                 at (x,y,z)=(%.11g,%.11g,%.11g)\n",
         global_qmax1, bi1, glob_xmax1,glob_ymax1,glob_zmax1);
  printf("BNSdata_analyze: global qmax2=%g in box%d\n"
         "                 at (x,y,z)=(%.11g,%.11g,%.11g)\n",
         global_qmax2, bi2, glob_xmax2,glob_ymax2,glob_zmax2);

  /* compute TOV */
  TOV_init(P_core1, kappa, Gamma, 0,
           &TOV_rf_surf1, &TOV_m1, &TOV_Phic1, &TOV_Psic1, &TOV_m01);
  TOV_init(P_core2, kappa, Gamma, 0,
           &TOV_rf_surf2, &TOV_m2, &TOV_Phic2, &TOV_Psic2, &TOV_m02);
  TOV_Psis1 = 1.0 + TOV_m1/(2.0*TOV_rf_surf1);
  TOV_Psis2 = 1.0 + TOV_m2/(2.0*TOV_rf_surf2);
  TOV_r_surf1 = TOV_rf_surf1 * TOV_Psis1*TOV_Psis1;
  TOV_r_surf2 = TOV_rf_surf2 * TOV_Psis2*TOV_Psis2;
/*
printf("TOV_m1=%g TOV_r_surf1=%g TOV_Psic1=%g\n",
TOV_m1,TOV_r_surf1, TOV_Psic1);
printf("TOV_m1=%g TOV_r_surf1=%g TOV_Psis1=%g\n",
TOV_m1,TOV_r_surf1, TOV_Psis1);
*/
  /* write into file */
  filenamelen = strlen(outdir) + strlen(name) + 200;
  filename = cmalloc(filenamelen+1);
  snprintf(filename, filenamelen, "%s/%s", outdir, name);
  fp = fopen(filename, "a");
  if(!fp) errorexits("failed opening %s", filename);
  fprintf(fp, "BNS data properties (time = %g):\n", grid->time);
  fprintf(fp, "-------------------\n");
  fprintf(fp, "n\t\t%.19g\n", n);
  fprintf(fp, "Gamma\t\t%.19g\n", Gamma);
  fprintf(fp, "kappa\t\t%.19g\n", kappa);
  fprintf(fp, "x_CM\t\t%.19g\n", x_CM);
  fprintf(fp, "Omega\t\t%.19g\n", Omega);
  fprintf(fp, "m01\t\t%.19g\n", m01);
  fprintf(fp, "m02\t\t%.19g\n", m02);
  fprintf(fp, "J_ADM\t\t%.19g\n", J_ADM);
  fprintf(fp, "M_ADM\t\t%.19g\n", M_ADM);
  for(i=b1n1-1; i>=b1n1-5 && i>=0; i--)
    fprintf(fp, "MpADM[%d]\t%.19g\t(A=%.16g  x=%.16g)\n", i, M_ADM_ofA[i],
            grid->box[1]->v[iX][i], grid->box[1]->v[ix][i]);
  fprintf(fp, "\n");
  fprintf(fp, "(m1/R)_inf\t%.19g\n", TOV_m1/TOV_r_surf1);
  fprintf(fp, "(m01/R)_inf\t%.19g\n", m01/TOV_r_surf1);
  fprintf(fp, "xin1\t\t%+.19g\n", xin1);
  fprintf(fp, "xmax1\t\t%+.19g\n", xmax1);
  fprintf(fp, "xout1\t\t%+.19g\n", xout1);
  fprintf(fp, "qmax1\t\t%.19g\n", qmax1);
  fprintf(fp, "\n");
  fprintf(fp, "(m2/R)_inf\t%.19g\n", TOV_m2/TOV_r_surf2);
  fprintf(fp, "(m02/R)_inf\t%.19g\n", m02/TOV_r_surf2);
  fprintf(fp, "xin2\t\t%+.19g\n", xin2);
  fprintf(fp, "xmax2\t\t%+.19g\n", xmax2);
  fprintf(fp, "xout2\t\t%+.19g\n", xout2);
  fprintf(fp, "qmax2\t\t%.19g\n", qmax2);
  fprintf(fp, "\n");
  fprintf(fp, "BNSdata_b\t%.19g\n", BNSdata_b);
  fprintf(fp, "sigp_00\t\t%+.19g\n", sigp_00);
  fprintf(fp, "sigp_10\t\t%+.19g\n", sigp_10);
  fprintf(fp, "sigm_00\t\t%+.19g\n", sigm_00);
  fprintf(fp, "sigm_10\t\t%+.19g\n", sigm_10);
  fprintf(fp, "\n");
  fclose(fp);
  free(filename);
  free(M_ADM_ofA);
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
  set_BNSdata_Sigma_BC(vlFu,vlu,  vlc2,vlc2,vluDerivs, 1);
}

/* evaluate linearized BNSdata eqns */
void J_BNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  BNS_CTS(vlJdu,vlu,  vlJdu,vldu,vlduDerivs, 0);
        /* ^--not used by BNS_CTS if nonlin=0 */
  /* BCs */
  set_BNSdata_BCs(vlJdu, vldu, vlduDerivs, 0);
  set_BNSdata_Sigma_BC(vlJdu,vlu,  vlJdu,vldu,vlduDerivs, 0);
}


/* evaluate eqn for a SINGLE one comp. var vlw */
void F_oneComp(tVarList *vlFw, tVarList *vlw,
               tVarList *vlwDerivs, tVarList *vlc2)
{
  /* local varlists with the correct grid, with ind copied from global vl */
  tVarList *lvlFu = vlalloc(vlw->grid);
  tVarList *lvlu  = vlalloc(vlw->grid);
  tVarList *lvlJdu = vlalloc(vlw->grid);
  tVarList *lvldu  = vlalloc(vlw->grid);
  tVarList *lvlduDerivs = vlalloc(vlw->grid);
  vlpushvl(lvlFu, vlFu);
  vlpushvl( lvlu,  vlu);
  vlpushvl(lvlJdu, vlJdu);
  vlpushvl( lvldu,  vldu);
  vlpushvl(lvlduDerivs, vlduDerivs);
 
  /* Note: lvlFu,lvlu contains vlFw,vlw */
  BNS_CTS(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 1);
                    /* ^------^-----^--------not used by BNS_CTS if nonlin=1 */
  /* BCs */
  set_BNSdata_BCs(vlFw, vlw, vlwDerivs, 1);
  set_BNSdata_Sigma_BC(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 1);
  
  /* free local varlists */
  vlfree(lvlFu);
  vlfree(lvlu);
  vlfree(lvlJdu);
  vlfree(lvldu);
  vlfree(lvlduDerivs);
}

/* evaluate linearized eqn for a SINGLE one comp. var for vldw */
void J_oneComp(tVarList *vlJdw, tVarList *vldw,
               tVarList *vldwDerivs, tVarList *vlw)
{
  /* local varlists with the correct grid, with ind copied from global vl */
  tVarList *lvlFu = vlalloc(vldw->grid);
  tVarList *lvlu  = vlalloc(vldw->grid);
  tVarList *lvlJdu = vlalloc(vldw->grid);
  tVarList *lvldu  = vlalloc(vldw->grid);
  tVarList *lvlduDerivs = vlalloc(vldw->grid);
  vlpushvl(lvlFu, vlFu);
  vlpushvl( lvlu,  vlu);
  vlpushvl(lvlJdu, vlJdu);
  vlpushvl( lvldu,  vldu);
  vlpushvl(lvlduDerivs, vlduDerivs);

  /* Note: vlJdu,vldu contains vlJdw,vldw */
  BNS_CTS(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 0);
        /* ^--not used by BNS_CTS if nonlin=0 */
  /* BCs */
  set_BNSdata_BCs(vlJdw, vldw, vldwDerivs, 0);
  set_BNSdata_Sigma_BC(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 0);

  /* free local varlists */
  vlfree(lvlFu);
  vlfree(lvlu);
  vlfree(lvlJdu);
  vlfree(lvldu);
  vlfree(lvlduDerivs);
}


/* set BCs for a varlist */
/* to get old version
void set_BNSdata_BCs_old(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
get BNSdata.c from Dec 16 2010
*/


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

/* set BNSdata_A/B/phi from x,y,z in box4/5 */
void set_BNSdata_ABphi(tGrid *grid)
{
  double A, B, phi;
  int b, ind, i03, pl, i,j,k, k_phi;

  for(b=4; b<=5; b++) /* go only over b=4,5 */
  {
    tBox *box = grid->box[b];
    tBox *box03 = grid->box[3*(b==4)]; /* box03 = box0 if b=5, box3 if b=4 */ 
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *pA, *pB, *pphi;

    enablevar_inbox(box, Ind("BNSdata_A"));
    enablevar_inbox(box, Ind("BNSdata_B"));
    enablevar_inbox(box, Ind("BNSdata_phi"));
    pA = box->v[Ind("BNSdata_A")];
    pB = box->v[Ind("BNSdata_B")];
    pphi = box->v[Ind("BNSdata_phi")];
    forallpoints(box, ind)
    {
      phi   = Arg(pY[ind],pZ[ind]);   if(phi<0) phi = 2.0*PI+phi;
      k_phi = box03->n3 * phi/(2.0*PI);
      nearestXYZ_of_xyz_inplane(box03, &i03, &A,&B,&phi,
                                pX[ind],pY[ind],pZ[ind], 3, k_phi);
      compute_ABphi_from_xyz(box03, &A,&B,&phi, pX[ind],pY[ind],pZ[ind]);
      /* do it twice in case it fails and gets A or B = 0 or 1 */
      compute_ABphi_from_xyz(box03, &A,&B,&phi, pX[ind],pY[ind],pZ[ind]);
      pA[ind] = A;
      pB[ind] = B;
      pphi[ind] = phi;
    }
    /* go over outer faces of box4/5 again, but in a better way */
    for(pl=0; pl<n1; pl=pl+n1-1)
    {
      int ind=Index(pl,0,0);
      phi   = Arg(pY[ind],pZ[ind]);   if(phi<0) phi = 2.0*PI+phi;
      k_phi = box03->n3 * phi/(2.0*PI);
      nearestXYZ_of_xyz_inplane(box03, &i03, &A,&B,&phi,
                                pX[ind],pY[ind],pZ[ind], 3, k_phi);
      forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
      {
        ind=Index(i,j,k);
        compute_ABphi_from_xyz(box03, &A,&B,&phi, pX[ind],pY[ind],pZ[ind]);
        pA[ind] = A;
        pB[ind] = B;
        pphi[ind] = phi;
      }
    }
    for(pl=0; pl<n2; pl=pl+n2-1)
    {
      int ind=Index(0,pl,0);
      phi   = Arg(pY[ind],pZ[ind]);   if(phi<0) phi = 2.0*PI+phi;
      k_phi = box03->n3 * phi/(2.0*PI);
      nearestXYZ_of_xyz_inplane(box03, &i03, &A,&B,&phi,
                                pX[ind],pY[ind],pZ[ind], 3, k_phi);
      forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
      {
        ind=Index(i,j,k);
        compute_ABphi_from_xyz(box03, &A,&B,&phi, pX[ind],pY[ind],pZ[ind]);
        pA[ind] = A;
        pB[ind] = B;
        pphi[ind] = phi;
      }
    }
    for(pl=0; pl<n3; pl=pl+n3-1)
    {
      int ind=Index(0,0,pl);
      phi   = Arg(pY[ind],pZ[ind]);   if(phi<0) phi = 2.0*PI+phi;
      k_phi = box03->n3 * phi/(2.0*PI);
      nearestXYZ_of_xyz_inplane(box03, &i03, &A,&B,&phi,
                                pX[ind],pY[ind],pZ[ind], 3, k_phi);
      forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
      {
        ind=Index(i,j,k);
        compute_ABphi_from_xyz(box03, &A,&B,&phi, pX[ind],pY[ind],pZ[ind]);
        pA[ind] = A;
        pB[ind] = B;
        pphi[ind] = phi;
      }
    }


  }
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

/* find residual of vars listed in string */
/* this works only it the string contains at most
   "BNSdata_Psi BNSdata_Bx BNSdata_By BNSdata_Bz BNSdata_alphaP BNSdata_Sigma" 
   */
double GridL2Norm_of_vars_in_string(tGrid *grid, char *str)
{
  int pos;
  char *word;
  double norm;
  double sum=0.0;
    
  word = cmalloc(strlen(str) + 10);
  pos=0;
  while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* make new vlw, ... for var in string word */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      
    F_oneComp(vlFw, vlw, vlwDerivs, NULL);
    norm = GridL2Norm(vlFw);
    sum += norm;
    //printf("%s: residual = %g\n", word, norm);
     
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);
  }
  //printf(" => total residual = %g\n", sum);
  free(word);
  return sum;
}

/* find residual of vars listed in string, but set Err to zero in box1&2 */
/* this works only if the string contains at most
   "BNSdata_Psi BNSdata_Bx BNSdata_By BNSdata_Bz BNSdata_alphaP BNSdata_Sigma" 
   */
double GridL2Norm_of_vars_in_string_withZeroErr_inbox12(tGrid *grid, char *str)
{
  int pos;
  char *word;
  double norm;
  double sum=0.0;
    
  word = cmalloc(strlen(str) + 10);
  pos=0;
  while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
  {
    int b;
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* make new vlw, ... for var in string word */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      
    F_oneComp(vlFw, vlw, vlwDerivs, NULL);
    /* set Error in box1 and 2 to zero */
    for(b=1; b<=2; b++)
    {
      tBox *box = grid->box[b];
      int iVar  = vlFw->index[0];
      int i;
      forallpoints(box, i)  box->v[iVar][i] = 0.0;
    }
    norm = GridL2Norm(vlFw);
    sum += norm;
    //printf("%s: residual = %g\n", word, norm);
     
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);
  }
  //printf(" => total residual = %g\n", sum);
  free(word);
  return sum;
}

/* solve some ell. eqns one after an other in the order given in a string */
int BNS_Eqn_Iterator_for_vars_in_string(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr, char *str)
{
  int    Newton_itmax = itmax;    /* Geti("BNSdata_Newton_itmax"); */
  double Newton_tol   = tol*0.1;  /* Getd("BNSdata_Newton_tol"); */
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int it, pos;
  int prN=pr;
  char *word;

  if(pr)
  { 
    printf("BNS_Eqn_Iterator_for_vars_in_string:\n"); 
    printf("%s\n", str);
    printf("  starting iterations, itmax=%d tol=%g\n", itmax, tol);
  }
  for (it = 0; it < itmax; it++)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* compute error */
    *normres = GridL2Norm_of_vars_in_string(grid, str);

    if(pr)
    {
      prdivider(1);
      printf("BNS_Eqn_Iterator_for_vars_in_string step %d residual = %.4e\n", it, *normres);
      fflush(stdout);
    }
    if (*normres <= tol) break;

    /* set Newton_tol */
    Newton_tol = (*normres)*0.05;

    /* go through BNSdata_Eqn_Iterator_order and solve for all vars in there */
    word = cmalloc(strlen(str) + 10);
    pos=0;
    while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
    {
      /* make new vlw, ... for var in string word */
      make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
               &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      /* call Newton solver for Psi */
      prdivider(1);
      printf("Solving elliptic Eqn for %s:\n", word);
      Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
             Newton_itmax, Newton_tol, &normresnonlin, prN,
             linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
             linSolver_itmax, linSolver_tolFac, linSolver_tol);
      free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                             vldw, vldwDerivs, vlJdw);
    }
    free(word);
  }
  /* warn if we didn't converge */
  if (it >= itmax)
  {
    *normres = GridL2Norm_of_vars_in_string(grid, str);
    prdivider(1);
    printf("BNS_Eqn_Iterator_for_vars_in_string:\n  *** Too many steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
    printf("BNS_Eqn_Iterator_for_vars_in_string:\n  Residual after %d steps:"
           "  residual = %e\n", it, *normres);
  }
  return it;
}

/* solve the coupled ell. eqns one after an other (in a particular order), 
   and iterate */
int BNS_ordered_Eqn_Iterator(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  prdivider(0);
  printf("BNS_ordered_Eqn_Iterator: using BNS_Eqn_Iterator_for_vars_in_string\n");
  return BNS_Eqn_Iterator_for_vars_in_string(grid, itmax, tol, normres,
                 linear_solver, 1, Gets("BNSdata_CTS_Eqs_Iteration_order"));
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
  int    Newton_itmax = itmax;    /* Geti("BNSdata_Newton_itmax"); */
  double Newton_tol   = tol*0.1;  /* Getd("BNSdata_Newton_tol"); */
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int it;
  int prN=pr;

  if(pr) printf("BNS_Eqn_Iterator:  starting iterations, itmax=%d tol=%g\n",
                itmax, tol); 
  for (it = 0; it < itmax; it++)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* compute error vlFu = F(u) */
    F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
    *normres = GridL2Norm(vlFu);
    if(pr)
    {
      prdivider(1);
      printf("BNS_Eqn_Iterator step %d residual = %.4e\n", it, *normres);
      fflush(stdout);
    }
    if (*normres <= tol) break;

    /* set Newton_tol */
    Newton_tol = (*normres)*0.05;

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

  }
  /* warn if we didn't converge */
  if (it >= itmax)
  {
    F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
    *normres = GridL2Norm(vlFu);
    prdivider(1);
    printf("BNS_Eqn_Iterator: *** Too many steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
    printf("BNS_Eqn_Iterator: Residual after %d steps:"
           "  residual = %e\n", it, *normres);
  }
  return it;
}

/* Solve the ell. eqns one after an other, as if they are uncoupled
   (as in Pedro's paper).
   sequence1 order is: Psi, Bx, By, Bz, alphaP, Sigma */
int BNS_Eqn_sequence1(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int prN=pr;
  tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

  /* compute error vlFu = F(u) */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence1: before sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }

  /* make new vlw, ... for Psi */
  make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
           &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Psi");
  /* call Newton solver for Psi */
  prdivider(1);
  if(pr) printf("Solving elliptic Eqn for BNSdata_Psi:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bx:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_By:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bz:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_alphaP:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Sigma:\n");
  Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
         Newton_itmax, Newton_tol, &normresnonlin, prN,
         linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);
  free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                         vldw, vldwDerivs, vlJdw);

  /* check res again */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence1:  after sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }
  return 1;
}

/* Solve the ell. eqns one after an other, as if they are uncoupled
   (as in Pedro's paper).
   sequence2 order is: Psi, alphaP, Bx, By, Bz, Sigma */
int BNS_Eqn_sequence2(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int prN=pr;
  tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

  /* compute error vlFu = F(u) */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence2: before sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }

  /* make new vlw, ... for Psi */
  make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
           &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Psi");
  /* call Newton solver for Psi */
  prdivider(1);
  if(pr) printf("Solving elliptic Eqn for BNSdata_Psi:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_alphaP:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bx:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_By:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bz:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Sigma:\n");
  Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
         Newton_itmax, Newton_tol, &normresnonlin, prN,
         linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);
  free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                         vldw, vldwDerivs, vlJdw);

  /* check res again */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence2:  after sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }
  return 1;
}

/* Solve the ell. eqns one after an other, as if they are uncoupled
   (as in Pedro's paper).
   sequence3 order is: alphaP, Bx, By, Bz, Psi, Sigma */
int BNS_Eqn_sequence3(tGrid *grid, 
  int Newton_itmax, double Newton_tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  int    linSolver_itmax  = Geti("BNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("BNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("BNSdata_linSolver_tol");
  double normresnonlin;
  int prN=pr;
  tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

  /* compute error vlFu = F(u) */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence3: before sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }

  /* make new vlw, ... for alphaP */
  make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
           &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_alphaP");
  /* call Newton solver for alphaP */
  prdivider(1);
  if(pr) printf("Solving elliptic Eqn for BNSdata_alphaP:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bx:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_By:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Bz:\n");
  Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
         Newton_itmax, Newton_tol, &normresnonlin, prN,
         linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);
  free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                         vldw, vldwDerivs, vlJdw);

  /* make new vlw, ... for Psi */
  make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
           &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "BNSdata_Psi");
  /* call Newton solver for Psi */
  prdivider(1);
  if(pr) printf("Solving elliptic Eqn for BNSdata_Psi:\n");
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
  if(pr) printf("Solving elliptic Eqn for BNSdata_Sigma:\n");
  Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
         Newton_itmax, Newton_tol, &normresnonlin, prN,
         linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
         linSolver_itmax, linSolver_tolFac, linSolver_tol);
  free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                         vldw, vldwDerivs, vlJdw);

  /* check res again */
  F_BNSdata(vlFu, vlu, vluDerivs, vlJdu);
  *normres = GridL2Norm(vlFu);
  if(pr)
  {
    printf("BNS_Eqn_sequence3:  after sequence: "
           "residual = %.4e\n", *normres);
    fflush(stdout);
  }
  return 1;
}


/* compute rest mass in box bi = 0 or 3 */
/* rest mass m_0 = \int d^3x \sqrt{g^{(3)}} rho_0 u^0 (-n_0),
   where n_0 = -\alpha, and d^3x \sqrt{g^{(3)}} = d^3x Psi^6  */
double GetInnerRestMass(tGrid *grid, int bi)
{
  int iInteg = Ind("BNSdata_temp1");

  /* set BNSdata_temp1 = Integ = rho_0 u^0 alpha * Psi^6 */
  BNS_set_restmassintegrand(grid, iInteg);
  /* Integ is integrand for rest mass. */
  /* get rest masses */
  return InnerVolumeIntegral(grid, bi, iInteg);
}

/* guess error in m01 from inner Volume int., but without adjusting
   surfaces */
void m01_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double m01;
  t_grid_grid0_m01_m02_struct *pars;
  
  /* get pars */
  pars = (t_grid_grid0_m01_m02_struct *) p;
      
  Setd("BNSdata_C1", vec[1]);
  BNS_compute_new_centered_q(pars->grid);
  m01 = GetInnerRestMass(pars->grid, 0);
  fvec[1] = m01 - pars->m01;
}

/* guess error in m02 from inner Volume int., but without adjusting
   surfaces */
void m02_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double m02;
  t_grid_grid0_m01_m02_struct *pars;

  /* get pars */
  pars = (t_grid_grid0_m01_m02_struct *) p;

  Setd("BNSdata_C2", vec[1]);
  BNS_compute_new_centered_q(pars->grid);
  m02 = GetInnerRestMass(pars->grid, 3);
  fvec[1] = m02 - pars->m02;
}

/* compute differences m01/2 - BNSdata_m01/2 */
void m0_errors_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid2;
  double BNSdata_n = Getd("BNSdata_n");
  double kappa     = Getd("BNSdata_kappa");
  double m01, m02;
  int b, i;
  int n1 = grid->box[1]->n1;
  int n2 = grid->box[1]->n2;
  double *q_b1 = grid->box[1]->v[Ind("BNSdata_q")];
  double *q_b2 = grid->box[2]->v[Ind("BNSdata_q")];

  /* set C1/2 */
  Setd("BNSdata_C1", vec[1]);
  Setd("BNSdata_C2", vec[2]);

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* reset sigma such that q=0 is at A=0 for box0/1 and box3/2 */
  if(q_b1[Index(n1-1,n2-1,0)]<0.0)
    reset_Coordinates_AnsorgNS_sigma_pm(grid, grid2, 0, 1);
  if(q_b2[Index(n1-1,n2-1,0)]<0.0)
    reset_Coordinates_AnsorgNS_sigma_pm(grid, grid2, 3, 2);

  /* initialize coords on grid2 */
  BNSgrid_init_Coords(grid2);

  /* interpolate q (and maybe some other vars) from grid onto new grid2 */
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Psi"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_alphaP"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bx"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_By"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Bz"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_Sigma"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBx"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBy"));
  Interpolate_Var_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_wBz"));
  varcopy(grid2, Ind("BNSdata_q"), Ind("BNSdata_qg"));

//  /* set q to zero if q<0 or in region 1 and 2 */
//  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);

  /***************************************/
  /* compute rest mass error Delta_m01/2 */
  /***************************************/
  /* get rest masses */
  m01 = GetInnerRestMass(grid, 0);
  m02 = GetInnerRestMass(grid, 3);

  printf("m0_errors_VectorFuncP: C1=%g C2=%g  m01=%g m02=%g\n",
         vec[1], vec[2], m01, m02);  fflush(stdout);
//write_grid(grid);

  fvec[1] = m01 - pars->m01;
  fvec[2] = m02 - pars->m02;
}

/* compute the new q and adjust the shape of the boundary between domain0/1
   or domain3/2 accordingly. This func modifies grid. It always
   interpolate vars from grid0 which is left untouched. */
void compute_new_q_and_adjust_domainshapes_InterpFromGrid0(tGrid *grid, 
                                                           tGrid *grid0, 
                                                           int innerdom)
{
  tGrid *grid2;
  int interp_qgold = !Getv("BNSdata_new_q", "FromFields");
  int outerdom;
  void (*Interp_From_Grid1_To_Grid2)(tGrid *grid1, tGrid *grid2, int vind,
                                     int innerdom);

  if(innerdom==0)      outerdom=1;
  else if(innerdom==3) outerdom=2;  
  else errorexit("compute_new_q_and_adjust_domainshapes: "
                 "innerdom is not 0 or 3");

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* make new grid2, which is an exact copy of grid */
  grid2 = make_empty_grid(grid->nvariables, 0);
  copy_grid(grid, grid2, 0);

  /* reset sigma such that q=0 at A=0 */
  reset_Coordinates_AnsorgNS_sigma_pm(grid, grid2, innerdom, outerdom);
  /* NOTE: coords of grid2 are initialized in if clause below. */

  /* do we make changes on both sides of grid? */
  if(Getv("BNSdata_adjust_domainshapes_Grid1_To_Grid2_Interpolator",
          "Interpolate_Var_From_Grid1_To_Grid2_wrapper"))
  {
    /* initialize coords of grid2 on both sides of grid */
    BNSgrid_init_Coords(grid2);
    /* NOTE: next we have to use Interpolate_Var_From_Grid1_To_Grid2_wrapper
       We have to interpolate on both sides of the grid since 
       BNSgrid_init_Coords changes box4/5 on both sides!!! */
    Interp_From_Grid1_To_Grid2 = Interpolate_Var_From_Grid1_To_Grid2_wrapper;
  }
  else
  {
    /* initialize coords of grid2 on side of innerdom */
    BNSgrid_init_Coords_pm(grid2, innerdom);
    /* use interpolator that does only side of innerdom */
    Interp_From_Grid1_To_Grid2 = Interp_Var_From_Grid1_To_Grid2_pm;
  }
  /* interpolate q (and maybe some other vars) from grid onto new grid2 */
  //  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qg"),innerdom);
  //  Interp_From_Grid1_To_Grid2(grid, grid2, Ind("BNSdata_qgold"),innerdom);
  Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_Psi"),innerdom);
  Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_alphaP"),innerdom);
  Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_Bx"),innerdom);
  Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_By"),innerdom);
  Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_Bz"),innerdom);
  if( (innerdom==0 && !Getv("BNSdata_rotationstate1","corotation")) ||
      (innerdom==3 && !Getv("BNSdata_rotationstate2","corotation"))   )
  {
    Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_Sigma"),innerdom);
    Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_wBx"),innerdom);
    Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_wBy"),innerdom);
    Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_wBz"),innerdom);
  }
  if(interp_qgold)
    Interp_From_Grid1_To_Grid2(grid0, grid2, Ind("BNSdata_qgold"),innerdom);

  BNS_compute_new_centered_q(grid2);

//  /* set q to zero if q<0 or in region 1 and 2 */
//  set_Var_to_Val_if_below_limit_or_inbox12(grid, Ind("BNSdata_q"), 0.0, 0.0);

  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);
}
/* compute the new q and adjust the shape of the boundary between domain0/1
   or domain3/2 accordingly. This is a wrapper that modifies grid */
void compute_new_q_and_adjust_domainshapes(tGrid *grid, int innerdom)
{
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid, innerdom);
}

/* compute difference m01 - BNSdata_m01 */
void m01_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  double m01;

  /* set C1 */
  Setd("BNSdata_C1", vec[1]);

  /* adjust grid so that new q=0 is at A=0 */
  //compute_new_q_and_adjust_domainshapes(grid, 0);
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid0, 0);

  /*************************************/
  /* compute rest mass error Delta_m01 */
  /*************************************/
  /* get rest mass */
  m01 = GetInnerRestMass(grid, 0);

  printf("m01_error_VectorFuncP: C1=%.13g  m01=%.13g\n", vec[1], m01);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);

  fvec[1] = m01 - pars->m01;
}

/* compute difference m02 - BNSdata_m02 */
void m02_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  double m02;

  /* set C2 */
  Setd("BNSdata_C2", vec[1]);

  /* adjust grid so that new q=0 is at A=0 */
  //compute_new_q_and_adjust_domainshapes(grid, 3);
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid0, 3);

  /*************************************/
  /* compute rest mass error Delta_m01 */
  /*************************************/
  /* get rest mass */
  m02 = GetInnerRestMass(grid, 3);

  printf("m02_error_VectorFuncP: C2=%.13g  m02=%.13g\n", vec[1], m02);
  fflush(stdout);
//grid->time=-200;
//write_grid(grid);
  
  fvec[1] = m02 - pars->m02;
}


/* find value fvec[1] of dq on x-axis at X=vec[1] by interpolation */
/* Note: before we can use this, we have to set
   grid->box[0]->v[Ind("BNSdata_temp3")][0] = bi;
   grid->box[0]->v[Ind("BNSdata_temp3")][2] = Y; */
void dq_dx_along_x_axis_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = central_q_errors_VectorFunc__grid;
  int bi = grid->box[0]->v[Ind("BNSdata_temp3")][0]; /* I abused BNSdata_temp3 to store bi,Y */
  tBox *box = grid->box[bi];
  double *c = box->v[Ind("BNSdata_temp2")]; /* coeffs of dq_dx */
  double X = vec[1];
  double Y = grid->box[0]->v[Ind("BNSdata_temp3")][2];
  fvec[1] = spec_interpolate(box, c, X,Y,0);
//  printf("     at (bi=%d X=%g Y=%g) dq_dx=%g\n", bi, X,Y, fvec[1]);
}

/* find max q along x-axis */
void find_Varmax_along_x_axis_usingBNSdata_temp123(tGrid *grid, int varind, 
                                                int *bi, double *X, double *Y)
{
  int b;
  int bi_guess = *bi;
  double Xvec[2];
  int stat, check;
  int i, ai, imax;
  double *qa0, *qa1, *qa2;
  int bq0, bq2;
  tBox *box;
  double *q, *dq, *c;
  double *Xp, *c0;

  /* get deriv dq of q in all boxes in BNSdata_temp1
     and dq's coeffs c in BNSdata_temp2, plus q's coeffs c0 in BNSdata_temp3 */
  forallboxes(grid, b)
  {
    box = grid->box[b];
    q = box->v[varind];
    dq= box->v[Ind("BNSdata_temp1")];
    c = box->v[Ind("BNSdata_temp2")];
    c0= box->v[Ind("BNSdata_temp3")];
    spec_Deriv1(box, 1, q, dq);
    spec_Coeffs(box, dq, c);
    spec_Coeffs(box,  q, c0);
  }

  *bi = -1; /* boxindex with dq=0 not yet found */

  /* look at NS1 */
  if(bi_guess<=1 || bi_guess==5)
  {
    /* index of the 2 boxes where we look */
    bq0 = 0;
    bq2 = 5;
  }
  /* look at NS2 */
  else if(bi_guess>=2 && bi_guess<=4)
  {
    /* index of the 2 boxes where we look */
    bq0 = 3;
    bq2 = 4;
  }
  else errorexit("bi_guess is wrong!!!");

  /* make 3 arrays (qa0, qa1, qa2) with q values along x-axis 
     in box0/3 at B=0 and B=1 and in box5/4 at y=z=0 */
  qa0 = dmalloc(grid->box[bq0]->n1);
  qa1 = dmalloc(grid->box[bq0]->n1);
  for(i=0; i<grid->box[bq0]->n1; i++)
  {
    box = grid->box[bq0];
    q = box->v[varind];
    qa0[i] = q[Ind_n1n2(i, 0,         0, box->n1,box->n2)];
    qa1[i] = q[Ind_n1n2(i, box->n2-1, 0, box->n1,box->n2)];
  }
  qa2 = dmalloc(grid->box[bq2]->n1);
  for(i=0; i<grid->box[bq2]->n1; i++)
  {
    box = grid->box[bq2];
    Xp = box->v[Ind("X")];
    c0 = box->v[Ind("BNSdata_temp3")];
    qa2[i] = spec_interpolate(box, c0, Xp[i],0,0);
  }

  /* find max in the 3 arrays */
  max3_in_1d_array(qa0,grid->box[bq0]->n1, qa1,grid->box[bq0]->n1,
                   qa2,grid->box[bq2]->n1, &ai, &imax);

  /* initial guess for max location in q (we need *bi, *X, *Y) */
  if(ai==0)  /* max at B=0 */
  {
    *bi = bq0;
    box = grid->box[*bi];
    Xp  = box->v[Ind("X")];
    /* if we are at the box boundary move one point in */
    if(imax==0)         imax++; 
    if(imax==box->n1-1) imax--;
    *X = Xp[Ind_n1n2(imax, 0, 0, box->n1,box->n2)]; 
    *Y = 0.0; /* <--B=0 */
  }
  else if(ai==1) /* max at B=1 */
  {
    *bi = bq0;
    box = grid->box[*bi];
    Xp  = box->v[Ind("X")];
    /* if we are at the box boundary move one point in */
    if(imax==0)         imax++; 
    if(imax==box->n1-1) imax--;
    *X = Xp[Ind_n1n2(imax, box->n2-1, 0, box->n1,box->n2)]; 
    *Y = 1.0; /* <--B=1 */
  }
  else if(ai==2) /* max in box5/4 */
  {
    *bi = bq2;
    box = grid->box[*bi];
    Xp  = box->v[Ind("X")];
    /* if we are at the box boundary move one point in */
    if(imax==0)         imax++; 
    if(imax==box->n1-1) imax--;
    *X = Xp[Ind_n1n2(imax, 0, 0, box->n1,box->n2)]; 
    *Y = 0.0; /* <--y=0 */
  }
  else errorexit("could not find max of q along x-axis");
//printf("bq0=%d  bq2=%d:   ai=%d  imax=%d\n", bq0,bq2, ai,imax);
//printf("  bi_guess=%d -> qmax is near: *bi=%d *X=%g *Y=%g\n",
//         bi_guess, *bi, *X, *Y);
//for(i=0; i<grid->box[bq0]->n1; i++)
//printf("B=0: i=%d  qa0[i]=%g\n", i,qa0[i]);
//for(i=0; i<grid->box[bq0]->n1; i++)
//printf("B=1: i=%d  qa1[i]=%g\n", i,qa1[i]);
//for(i=0; i<grid->box[bq2]->n1; i++)
//printf("y=0: i=%d  qa2[i]=%g\n", i,qa2[i]);

  /* free 3 arrays (qa0, qa1, qa2) with q values along x-axis */
  free(qa0);
  free(qa1);
  free(qa2);
//double *q = box->v[varind];
//double *Xp = box->v[Ind("X")];
//for(b=0;b<box->n1;b++)
//printf("X=%g q=%g dq=%g\n",
//Xp[Ind_n1n2(b, 0, 0, box->n1,box->n2)],
//q[Ind_n1n2(b, 0, 0, box->n1,box->n2)],
//dq[Ind_n1n2(b, 0, 0, box->n1,box->n2)]);

//  /* print guesses */
//  printf("  bi_guess=%d -> guesses: *bi=%d *X=%g *Y=%g\n",
//         bi_guess, *bi, *X, *Y);

  /* save *bi, *X, *Y in BNSdata_temp3 */
  grid->box[0]->v[Ind("BNSdata_temp3")][0] = *bi;
  grid->box[0]->v[Ind("BNSdata_temp3")][1] = *X;
  grid->box[0]->v[Ind("BNSdata_temp3")][2] = *Y;
  
  /* use newton_linesrch_its to find xmax1 */
  central_q_errors_VectorFunc__grid = grid;
  Xvec[1] = *X;
  stat = newton_linesrch_its(Xvec, 1, &check,
                             dq_dx_along_x_axis_VectorFunc, 1000, dequaleps);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  *X = Xvec[1];
  //printf("  bi_guess=%d -> qmax is at: *bi=%d *X=%g *Y=%g\n",
  //       bi_guess, *bi, *X, *Y);
}

void find_qmax1_along_x_axis(tGrid *grid, int *bi, double *X, double *Y)
{
  find_Varmax_along_x_axis_usingBNSdata_temp123(grid, Ind("BNSdata_q"),
                                                bi, X,Y);
}

/* find max q along x-axis, old version that looks only at dq/dx 
   at endpoints to get initial guess, which is fragile */
void find_qmax1_along_x_axis_old(tGrid *grid, int *bi, double *X, double *Y)
{
  int b;
  int bi_guess = *bi;
  double *dqin, *dqout;
  double dq1,dq2;
  double Xvec[2];
  int stat, check;

  /* get deriv dq of q in all boxes in BNSdata_temp1
     and dq's coeffs c in BNSdata_temp2 */
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *q = box->v[Ind("BNSdata_q")];
    double *dq= box->v[Ind("BNSdata_temp1")];
    double *c = box->v[Ind("BNSdata_temp2")];
    spec_Deriv1(box, 1, q, dq);
    spec_Coeffs(box, dq, c);
  }

  /* look at NS1 */
  if(bi_guess==0)
  {
    tBox *box;
    double *dq, *c;

    *bi = -1; /* boxindex with dq=0 not yet found */

    /* look in box0 at B=0 */
    box = grid->box[0];
    dq = box->v[Ind("BNSdata_temp1")];
    dq1 = dq[Ind_n1n2(2, 0, 0,         box->n1,box->n2)];
    dq2 = dq[Ind_n1n2(box->n1-1, 0, 0, box->n1,box->n2)];
//double *q = box->v[Ind("BNSdata_q")];
//double *Xp = box->v[Ind("X")];
//for(b=0;b<box->n1;b++)
//printf("X=%g q=%g dq=%g\n",
//Xp[Ind_n1n2(b, 0, 0, box->n1,box->n2)],
//q[Ind_n1n2(b, 0, 0, box->n1,box->n2)],
//dq[Ind_n1n2(b, 0, 0, box->n1,box->n2)]);
    if(dq1*dq2<=0)  { *bi = 0;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=0.0; /* <--B=0 */ }
    
    if(*bi<0) /* boxindex *bi with dq=0 not yet found */
    {
      /* look in box0 at B=1 */
      box = grid->box[0];
      dq = box->v[Ind("BNSdata_temp1")];
      dq1 = dq[Ind_n1n2(2, box->n2-1, 0,         box->n1,box->n2)];
      dq2 = dq[Ind_n1n2(box->n1-1, box->n2-1, 0, box->n1,box->n2)];
      if(dq1*dq2<=0)  { *bi = 0;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=1.0; /* <--B=1 */ }
    }
    if(*bi<0) /* boxindex *bi with dq=0 not yet found */
    {
      /* look in box5 */
      box = grid->box[5];
      c = box->v[Ind("BNSdata_temp2")];
      dq1 = spec_interpolate(box, c, box->bbox[0],0,0);
      dq2 = spec_interpolate(box, c, box->bbox[1],0,0);
      if(dq1*dq2<=0)  { *bi = 5;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=0.0; }
    }
    if(*bi<0) errorexit("could not find max of q in NS1");
  }
  /* look at NS2 */
  else if(bi_guess==3)
  {
    tBox *box;
    double *dq, *c;

    *bi = -1; /* boxindex with dq=0 not yet found */

    /* look in box3 at B=0 */
    box = grid->box[3];
    dq = box->v[Ind("BNSdata_temp1")];
    dq1 = dq[Ind_n1n2(2, 0, 0,         box->n1,box->n2)];
    dq2 = dq[Ind_n1n2(box->n1-1, 0, 0, box->n1,box->n2)];
    if(dq1*dq2<=0)  { *bi = 3;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=0.0; /* <--B=0 */ }
    
    if(*bi<0) /* boxindex *bi with dq=0 not yet found */
    {
      /* look in box3 at B=1 */
      box = grid->box[3];
      dq = box->v[Ind("BNSdata_temp1")];
      dq1 = dq[Ind_n1n2(2, box->n2-1, 0,         box->n1,box->n2)];
      dq2 = dq[Ind_n1n2(box->n1-1, box->n2-1, 0, box->n1,box->n2)];
      if(dq1*dq2<=0)  { *bi = 3;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=1.0; /* <--B=1 */ }
    }
    if(*bi<0) /* boxindex *bi with dq=0 not yet found */
    {
      /* look in box4 */
      box = grid->box[4];
      c = box->v[Ind("BNSdata_temp2")];
      dq1 = spec_interpolate(box, c, box->bbox[0],0,0);
      dq2 = spec_interpolate(box, c, box->bbox[1],0,0);
      if(dq1*dq2<=0)  { *bi = 4;  *X=0.5*(box->bbox[1]+box->bbox[0]);  *Y=0.0; }
    }
    if(*bi<0) errorexit("could not find max of q in NS2");
  }
  else errorexit("bi_guess has to be 0 or 3");

//  /* print guesses */
//  printf("  bi_guess=%d -> guesses: *bi=%d *X=%g *Y=%g\n",
//         bi_guess, *bi, *X, *Y);

  /* save *bi, *X, *Y in BNSdata_temp3 */
  grid->box[0]->v[Ind("BNSdata_temp3")][0] = *bi;
  grid->box[0]->v[Ind("BNSdata_temp3")][1] = *X;
  grid->box[0]->v[Ind("BNSdata_temp3")][2] = *Y;
  
  /* use newton_linesrch_its to find xmax1 */
  central_q_errors_VectorFunc__grid = grid;
  Xvec[1] = *X;
  stat = newton_linesrch_its(Xvec, 1, &check,
                             dq_dx_along_x_axis_VectorFunc, 1000, dequaleps);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  *X = Xvec[1];
  printf("  bi_guess=%d -> qmax is at: *bi=%d *X=%g *Y=%g\n",
         bi_guess, *bi, *X, *Y);
}


/* compute deviation in desired central q value and location */
void central_q_errors_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = central_q_errors_VectorFunc__grid;
  int bi1, bi2;
  double qmax1, Xmax1, Ymax1, xmax1;
  double qmax2, Xmax2, Ymax2, xmax2;

  /* set constants */
  Setd("BNSdata_C1", vec[1]);
  Setd("BNSdata_Omega", vec[2]);
  Setd("BNSdata_C2", vec[3]);
  Setd("BNSdata_x_CM", vec[4]);
  printf("central_q_errors_VectorFunc: C1=%.6g Omega=%.6g C2=%.6g xCM=%.6g\n",
         vec[1], vec[2], vec[3], vec[4]);

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);

  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  /* compute Cartesian xmax1 */
  if(bi1==0)
  {
    tBox *box = grid->box[bi1];
    xmax1 = box->x_of_X[1]((void *) box, -1, Xmax1,Ymax1,0);
  }
  else xmax1 = Xmax1;

  /* compute Cartesian xmax2 */
  if(bi2==3)
  {
    tBox *box = grid->box[bi2];
    xmax2 = box->x_of_X[1]((void *) box, -1, Xmax2,Ymax2,0);
  }
  else xmax2 = Xmax2;

  /* set fvec */
  fvec[1] = qmax1 - Getd("BNSdata_qmax1");
  fvec[2] = xmax1 - Getd("BNSdata_xmax1");
  fvec[3] = qmax2 - Getd("BNSdata_qmax2");
  fvec[4] = xmax2 - Getd("BNSdata_xmax2");

  /* print results */
  printf(" qmax1=%.6g xmax1=%g (X=%.6g Y=%.6g) fvec[1]=%g fvec[2]=%g\n",
         qmax1, xmax1, Xmax1, Ymax1, fvec[1], fvec[2]);
  printf(" qmax2=%.6g xmax2=%g (X=%.6g Y=%.6g) fvec[3]=%g fvec[4]=%g\n",
         qmax2, xmax2, Xmax2, Ymax2, fvec[3], fvec[4]);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);
}

/* compute deviation in q value at xmax */
void estimate_q_errors_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = central_q_errors_VectorFunc__grid;
  int bi1, bi2;
  double qmax1, Xmax1, Ymax1, xmax1;
  double qmax2, Xmax2, Ymax2, xmax2;

  /* set constants */
  Setd("BNSdata_C1", vec[1]);
  Setd("BNSdata_C2", vec[2]);
  printf("estimate_q_errors_VectorFunc: C1=%.6g C2=%.6g\n", vec[1], vec[2]);

  /* compute new q */
  BNS_compute_new_centered_q(grid);

  /* find max q locations xmax1/2 in NS1/2 */
  bi1=0;  bi2=3;
  find_qmax1_along_x_axis(grid, &bi1, &Xmax1, &Ymax1);
  find_qmax1_along_x_axis(grid, &bi2, &Xmax2, &Ymax2);

  /* compute qmax1 and qmax2 */
  qmax1 = BNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);
  qmax2 = BNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  /* compute Cartesian xmax1 */
  if(bi1==0)
  {
    tBox *box = grid->box[bi1];
    xmax1 = box->x_of_X[1]((void *) box, -1, Xmax1,Ymax1,0);
  }
  else xmax1 = Xmax1;

  /* compute Cartesian xmax2 */
  if(bi2==3)
  {
    tBox *box = grid->box[bi2];
    xmax2 = box->x_of_X[1]((void *) box, -1, Xmax2,Ymax2,0);
  }
  else xmax2 = Xmax2;

  /* set fvec */
  fvec[1] = qmax1 - Getd("BNSdata_qmax1");
  fvec[2] = qmax2 - Getd("BNSdata_qmax2");

  /* print results */
  printf(" qmax1=%.6g xmax1=%g (X=%.6g Y=%.6g) fvec[1]=%g\n",
         qmax1, xmax1, Xmax1, Ymax1, fvec[1]);
  printf(" qmax2=%.6g xmax2=%g (X=%.6g Y=%.6g) fvec[2]=%g\n",
         qmax2, xmax2, Xmax2, Ymax2, fvec[2]);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);
}


/* find derivative (fvec[1],fvec[2],fvec[3]) of dq 
   at (X,Y,Z)=(vec[1],vec[2],vec[3]) by interpolation if n=3
   otherwise just find dq/dx, dq/dy */
/* Note: before we can use this, we have to set par and "BNSdata_temp2" */
void gradient_q_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tBox *box = (tBox *) p;
  double *cx= box->v[Ind("BNSdata_temp1")]; /* coeffs of dq_dx */
  double *cy= box->v[Ind("BNSdata_temp2")];
  double *cz= box->v[Ind("BNSdata_temp3")];
  double X = vec[1];
  double Y = vec[2];
  double Z = 0;

  if(n==3)
  {
    Z = vec[3];
    fvec[3] = spec_interpolate(box, cz, X,Y,Z);
  }
  fvec[1] = spec_interpolate(box, cx, X,Y,Z);
  fvec[2] = spec_interpolate(box, cy, X,Y,Z);
//  printf("   at (bi=%d X=%g Y=%g Z=%g)\n"
//         "   dq_dx=%.12g dq_dy=%.12g dq_dz=%.12g\n",
//         box->b, X,Y,Z, fvec[1],fvec[2],fvec[3]);
}

/* find max q in or near box *bi, *bi has to be 0 or 3, return qmax */
double BNSdata_find_position_of_qmax(tGrid *grid, int *bi, 
                                     double *X, double *Y, double *Z)
{
  int b;
  int bi_guess = *bi;
  double Xvec[4];
  int stat, check;
  tBox *box;
  double *q, *cx,*cy,*cz, *c0;
  double qmax;

  /* get derivs of q in all boxes in BNSdata_qx/y/z
     and coeffs od derivs in cx/y/z in BNSdata_temp1/2/3, 
     plus q's coeffs c0 in BNSdata_temp4 */
  forallboxes(grid, b)
  {
    box = grid->box[b];
    q = box->v[Ind("BNSdata_q")];
    cx= box->v[Ind("BNSdata_temp1")];
    cy= box->v[Ind("BNSdata_temp2")];
    cz= box->v[Ind("BNSdata_temp3")];
    c0= box->v[Ind("BNSdata_temp4")];
    /* set coeffs of dq in BNSdata_temp1/2/3 */
    FirstDerivsOf_S(box, Ind("BNSdata_q"), Ind("BNSdata_qx"));
    spec_Coeffs(box, box->v[Ind("BNSdata_qx")], cx);
    spec_Coeffs(box, box->v[Ind("BNSdata_qy")], cy);
    spec_Coeffs(box, box->v[Ind("BNSdata_qz")], cz);
    /* set coeffs of q in BNSdata_temp4 */
    spec_Coeffs(box, q, c0);
  }

//  /* print guesses */
//  printf("  bi_guess=%d -> guesses: *bi=%d *X=%.15g *Y=%.15g *Z=%.15g\n",
//         bi_guess, *bi, *X, *Y, *Z);

  /* use newton_linesrch_itsP to find Xmax,Ymax, for Z=0
     use var box as parameters */
/*
  Xvec[1] = *X;
  Xvec[2] = *Y;
  box = grid->box[*bi];
  stat = newton_linesrch_itsP(Xvec, 2, &check,
                              gradient_q_VectorFuncP, (void *) box,
                              1000, dequaleps);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  *X = Xvec[1];
  *Y = Xvec[2];
  if(check || stat<0) printf("  --> X=%.15g *Y=%.15g\n", *X, *Y);
*/
//  /* print estimates */
//  printf("  bi_guess=%d -> estimates: *bi=%d *X=%.15g *Y=%.15g *Z=%.15g\n",
//         bi_guess, *bi, *X, *Y, *Z);


  /* set Xvec */
  Xvec[1] = *X;
  Xvec[2] = *Y;
  Xvec[3] = *Z;

  /* do we actually try to find Xvec[1],Xvec[2],Xvec[3]? */
  if(Getv("BNSdata_find_position_of_qmax","find_XYZ"))
  {
    /* move Y away from boundary in case we are in box0/3 */
    box = grid->box[*bi];
    if(dlesseq(Xvec[2],box->bbox[2]))    Xvec[2] += dequaleps*1000.0;
    if(dgreatereq(Xvec[2],box->bbox[3])) Xvec[2] -= dequaleps*1000.0;

    /* use newton_linesrch_itsP to find max, use var box as parameters */
    stat = newton_linesrch_itsP(Xvec, 3, &check,
                                gradient_q_VectorFuncP, (void *) box,
                                1000, dequaleps);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    if(check || stat<0) printf("  --> Xvec=(%.15g,%.15g,%.15g)\n",
                               Xvec[1],Xvec[2],Xvec[3]);
  }
  
  /* make sure Y is inside bounding box */
  if( Xvec[1]<box->bbox[0] ||
      Xvec[1]>box->bbox[1] ||
      Xvec[2]<box->bbox[2] ||
      Xvec[2]>box->bbox[3] )
  {
    printf("BNSdata_find_position_of_qmax: Xvec=(%.15g,%.15g)\nis outside box%d. "
           "Left *X,*Y,*Z unchanged.\n", Xvec[1],Xvec[2], *bi);
    return -1.0;
  }        
  else
  {     
    *X = Xvec[1];
    *Y = Xvec[2];
    *Z = Xvec[3];
  }

  c0= box->v[Ind("BNSdata_temp4")];
  qmax = spec_interpolate(box, c0, *X,*Y,*Z);
  //printf(" global qmax is at: *bi=%d *X=%.11g *Y=%.11g *Z=%.11g\n",
  //       *bi, *X, *Y, *Z);
  return qmax;
}

/* find cart coords of max of q */
double BNSdata_find_xyz_of_qmax(tGrid *grid, int *bi, 
                                double *x, double *y, double *z)
{
  double qmax, Xmax,Ymax,Zmax;

  /* find global max of q in NS1/2 */
  find_qmax1_along_x_axis(grid, bi, &Xmax, &Ymax);
  //printf("Ymax=%g\n",Ymax);
  Zmax=0.0;
  qmax = BNSdata_find_position_of_qmax(grid, bi, &Xmax, &Ymax, &Zmax);
  if(qmax<0.0)
  {
    if(*bi<=1 || *bi==5)
      *x = Getd("BNSdata_xmax1");
    else
      *x = Getd("BNSdata_xmax2");
    *y = *z = 0.0;
    return qmax;
  }
  
  if(grid->box[*bi]->x_of_X[1] != NULL)
  {
    *x = grid->box[*bi]->x_of_X[1]((void *) grid->box[*bi], -1, Xmax,Ymax,Zmax);
    *y = grid->box[*bi]->x_of_X[2]((void *) grid->box[*bi], -1, Xmax,Ymax,Zmax);
    *z = grid->box[*bi]->x_of_X[3]((void *) grid->box[*bi], -1, Xmax,Ymax,Zmax);
  }
  else
  {
    *x = Xmax;
    *y = Ymax;
    *z = Zmax;
  }
  return qmax;
}

/* set BNSdata_actual_x/y/z/max1/2 pars */
void set_BNSdata_actual_xyzmax_pars(tGrid *grid)
{
  int bi1, bi2;
  double x1,y1,z1, x2,y2,z2;
  double qmax1, qmax2;

  printf("set_BNSdata_actual_xyzmax_pars:\n");
  bi1=0;
  bi2=3;
  qmax1 = BNSdata_find_xyz_of_qmax(grid, &bi1, &x1,&y1,&z1);
  qmax2 = BNSdata_find_xyz_of_qmax(grid, &bi2, &x2,&y2,&z2);
  if(qmax1<0.0 || qmax2<0.0)
  {
    printf("set_BNSdata_actual_xyzmax_pars: BNSdata_find_xyz_of_qmax cannot find one\n");
    printf("of the maxima: qmax1=%g qmax2=%g\n", qmax1, qmax2);
    printf("Setting the following parameters to desired (not actual) values:\n");
    x1 = Getd("BNSdata_xmax1");
    x2 = Getd("BNSdata_xmax2");
    y1 = y2 = z1 = z2 = 0.0;
  }
  Setd("BNSdata_actual_xmax1", x1);
  Setd("BNSdata_actual_ymax1", y1);
  Setd("BNSdata_actual_zmax1", z1);
  Setd("BNSdata_actual_xmax2", x2);
  Setd("BNSdata_actual_ymax2", y2);
  Setd("BNSdata_actual_zmax2", z2);
  printf(" BNSdata_actual_xmax1 = %.13g\n", Getd("BNSdata_actual_xmax1"));
  printf(" BNSdata_actual_ymax1 = %.13g\n", Getd("BNSdata_actual_ymax1"));
  printf(" BNSdata_actual_zmax1 = %.13g\n", Getd("BNSdata_actual_zmax1"));
  printf(" BNSdata_actual_xmax2 = %.13g\n", Getd("BNSdata_actual_xmax2"));
  printf(" BNSdata_actual_ymax2 = %.13g\n", Getd("BNSdata_actual_ymax2"));
  printf(" BNSdata_actual_zmax2 = %.13g\n", Getd("BNSdata_actual_zmax2"));
}
