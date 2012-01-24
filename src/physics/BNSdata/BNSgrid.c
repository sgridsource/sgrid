/* grid_setup.c */
/* Wolfgang Tichy 2008 */


#include "sgrid.h"
#include "BNSdata.h"

#define Power pow

/* struct types used in root finder newton_linesrch_itsP */
typedef struct T_grid_b_struct {
  tGrid *grid; /* grid */
  int b;      /* box1 */
} t_grid_b_struct;

typedef struct T_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct {
  tGrid *grid; /* grid */
  tBox *box;   /* box for DelXR_of_AB_VectorFuncP */
  double X;    /* X for DelXR_of_AB_VectorFuncP */
  double R;    /* R for DelXR_of_AB_VectorFuncP */
  double phi;  /* phi for DelXR_of_AB_VectorFuncP and q_of_sigp_forgiven_BphiP */
  double sigp_1phi; /* sigp_1phi for q_of_sigp_forgiven_BphiP */
  double B;         /* B for q_of_sigp_forgiven_BphiP */
  int    icoeffs;   /* icoeffs for q_of_sigp_forgiven_BphiP */
  int    innerdom;  /* inner domain for q_of_sigp_forgiven_BphiP */
  int    outerdom;  /* outer domain for q_of_sigp_forgiven_BphiP */
} t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct;

/* from Coordinates.c */
extern double  (*Coordinates_AnsorgNS_sigmap)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmap_dB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmap_dphi)(tBox *box, int ind, double B, double phi);
extern double  (*Coordinates_AnsorgNS_sigmam)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dphi)(tBox *box, int ind, double B, double phi);
extern int Coordinates_AnsorgNS_b_ParIndex; /* index of par Coordinates_AnsorgNS_b */

/* global vars in this file */
double sigma1; /* initial sigma_+ */
double sigma2; /* initial sigma_- */
double rf_surf1; /* radius of star1 */
double rf_surf2; /* radius of star2 */
double P_core1;  /* core pressure of star1 */
double P_core2;  /* core pressure of star2 */
/* global vars in this file for root finding */
tBox *BNSdata_q_VectorFunc_box;      /* box for BNdata_q_VectorFunc */
double *BNSdata_q_VectorFunc_coeffs; /* coeffs for BNdata_q_VectorFunc */
double BNSdata_q_VectorFunc_B;       /* B for BNSdata_q_VectorFunc */
double BNSdata_q_VectorFunc_phi;     /* phi for BNSdata_q_VectorFunc */
double BNdata_sigp_VectorFunc_sigp_1phi; /* sigp_1phi for BNdata_sigp_VectorFunc */
double BNdata_sigp_VectorFunc_B;         /* B for BNdata_sigp_VectorFunc */
double BNdata_sigp_VectorFunc_X0;        /* X0 for BNdata_sigp_VectorFunc */
/* the q_of_sigp_forgiven_Bphi__... are needed only for old versions of code still here */
double q_of_sigp_forgiven_Bphi__sigp_1phi; /* sigp_1phi for q_of_sigp_forgiven_Bphi */
double q_of_sigp_forgiven_Bphi__B;         /* B for q_of_sigp_forgiven_Bphi */
double q_of_sigp_forgiven_Bphi__phi;       /* phi for q_of_sigp_forgiven_Bphi */
tGrid *q_of_sigp_forgiven_Bphi__grid;      /* grid for q_of_sigp_forgiven_Bphi */
int    q_of_sigp_forgiven_Bphi__icoeffs;   /* icoeffs for q_of_sigp_forgiven_Bphi */
int    q_of_sigp_forgiven_Bphi__innerdom;  /* inner domain for q_of_sigp_forgiven_Bphi */
int    q_of_sigp_forgiven_Bphi__outerdom;  /* outer domain for q_of_sigp_forgiven_Bphi */


/* funs in this file */
void rf_surf1_VectorFunc(int n, double *vec, double *fvec);
void rf_surf2_VectorFunc(int n, double *vec, double *fvec);
void m01_VectorFunc(int n, double *vec, double *fvec);
void m02_VectorFunc(int n, double *vec, double *fvec);




/* functions for sigma_{+-} */
double return_sigma1(tBox *box, int ind, double B, double phi)
{
  return sigma1;
}
double return_dsigma_zero(tBox *box, int ind, double B, double phi)
{
  return 0.0;
}
double return_sigma2(tBox *box, int ind, double B, double phi)
{
  return sigma2;
}

/* setup initial boxsizes */
int set_boxsizes(tGrid *grid)
{
  double m1, Phic1, Psic1;
  double m2, Phic2, Psic2;
  double kappa     = Getd("BNSdata_kappa");
  double BNSdata_n = Getd("BNSdata_n");
  double BNSdata_b = Getd("BNSdata_b");
  double m01 = Getd("BNSdata_m01");
  double m02 = Getd("BNSdata_m02");
  double Gamma     = 1.0 + 1.0/BNSdata_n;
  double xmin1,xmax1, xmin2,xmax2, xc1, xc2; /* x-positions of stars */
  double DoM, nu; /* distance over total rest mass, and rest mass ratio */
  double DoM3, DoM4, DoM5; /* powers of DoM */
  double xCM, Omega;
  double Fc, qc, Cc, oouzerosqr;
  double vec[2];
  int check;

  printf("set_boxsizes: setting box sizes and coordinates used ...\n");

  /* reset initial BNSdata_m01/2 if needed */
  if(Getv("BNSdata_iterate_m0", "yes"))
  {
    printf(" BNSdata_iterate_m0 = yes : setting:\n");
    /* set BNSdata_desired_m01/2 if needed */
    if(Getd("BNSdata_desired_m01")<0.0) Setd("BNSdata_desired_m01", m01);
    if(Getd("BNSdata_desired_m02")<0.0) Setd("BNSdata_desired_m02", m02);
    if(Getd("BNSdata_desired_kappa")<0.0) Setd("BNSdata_desired_kappa", kappa);
    printf("   BNSdata_desired_m01 = %g\n", Getd("BNSdata_desired_m01"));
    printf("   BNSdata_desired_m02 = %g\n", Getd("BNSdata_desired_m02"));
    printf("   BNSdata_desired_kappa = %g\n", Getd("BNSdata_desired_kappa"));
    /* set initial BNSdata_m01/2 */
    m01 = Getd("BNSdata_init_m01");
    m02 = Getd("BNSdata_init_m02");
    Setd("BNSdata_m01", m01);
    Setd("BNSdata_m02", m02);
    printf("   BNSdata_m01 = BNSdata_init_m01 = m01 = %g\n", m01);
    printf("   BNSdata_m02 = BNSdata_init_m02 = m02 = %g\n", m02);
    printf("   BNSdata_m0change = %g\n", Getd("BNSdata_m0change"));
    /* set initial kappa */
    if(Getv("BNSdata_adjustkappa", "kappa^(-n/2)*m0_EQ_desired_kappa^(-n/2)*m0"))
    {
      double tm01 = Getd("BNSdata_desired_m01");
      double tm02 = Getd("BNSdata_desired_m02");
      kappa = Getd("BNSdata_desired_kappa") *
              pow((tm01+tm02)/(m01+m02), -2/BNSdata_n);
      Setd("BNSdata_kappa", kappa);
      printf("   BNSdata_kappa = %g\n", kappa);
    }
  }

  /* set Coordinates_AnsorgNS_sigmap/m funcs, .s.t. it works with sigma1/2 */
  Coordinates_AnsorgNS_sigmap       = return_sigma1;
  Coordinates_AnsorgNS_dsigmap_dB   = return_dsigma_zero;
  Coordinates_AnsorgNS_dsigmap_dphi = return_dsigma_zero;
  Coordinates_AnsorgNS_sigmam       = return_sigma2;
  Coordinates_AnsorgNS_dsigmam_dB   = return_dsigma_zero;
  Coordinates_AnsorgNS_dsigmam_dphi = return_dsigma_zero;
  Setd("Coordinates_AnsorgNS_b", BNSdata_b);
  Coordinates_AnsorgNS_b_ParIndex = GetParIndex("Coordinates_AnsorgNS_b");

  /* set some box pars */
  if(Getv("BNSdata_grid", "4ABphi_2xyz"))
  {
    int innercube_n = Geti("BNSdata_2xyz_n");
    /* Sets("nboxes", "6");  // <--cannot be set here, needs to be set earlier */
    Seti("box4_n1", innercube_n);
    Seti("box4_n2", innercube_n);
    Seti("box4_n3", innercube_n);
    Seti("box5_n1", innercube_n);
    Seti("box5_n2", innercube_n);
    Seti("box5_n3", innercube_n);

    Sets("Coordinates_AnsorgNS_set_sigma_pm_pointers", "yes");

    Sets("box0_Coordinates", "AnsorgNS0");
    Sets("box0_basis1", "ChebExtrema");
    Sets("box0_min1", "0");
    Sets("box0_basis2", "ChebExtrema");
    Sets("box0_min2", "0");
    Sets("box0_max2", "1");
    Sets("box0_basis3", "Fourier");
    Sets("box0_min3", "0");
    Sets("box0_max3", "2*pi");

    Sets("box1_Coordinates", "AnsorgNS1");
    Sets("box1_basis1", "ChebExtrema");
    Sets("box1_min1", "0");
    Sets("box1_max1", "1");
    Sets("box1_basis2", "ChebExtrema");
    Sets("box1_min2", "0");
    Sets("box1_max2", "1");
    Sets("box1_basis3", "Fourier");
    Sets("box1_min3", "0");
    Sets("box1_max3", "2*pi");

    Sets("box2_Coordinates", "AnsorgNS2");
    Sets("box2_basis1", "ChebExtrema");
    Sets("box2_min1", "0");
    Sets("box2_max1", "1");
    Sets("box2_basis2", "ChebExtrema");
    Sets("box2_min2", "0");
    Sets("box2_max2", "1");
    Sets("box2_basis3", "Fourier");
    Sets("box2_min3", "0");
    Sets("box2_max3", "2*pi");
  
    Sets("box3_Coordinates", "AnsorgNS3");
    Sets("box3_basis1", "ChebExtrema");
    Sets("box3_min1", "0");
    Sets("box3_basis2", "ChebExtrema");
    Sets("box3_min2", "0");
    Sets("box3_max2", "1");
    Sets("box3_basis3", "Fourier");
    Sets("box3_min3", "0");
    Sets("box3_max3", "2*pi");

    Sets("box4_Coordinates", "Cartesian");
    Sets("box4_basis1", "ChebExtrema");
    Sets("box4_basis2", "ChebExtrema");
    Sets("box4_basis3", "ChebExtrema");

    Sets("box5_Coordinates", "Cartesian");
    Sets("box5_basis1", "ChebExtrema");
    Sets("box5_basis2", "ChebExtrema");
    Sets("box5_basis3", "ChebExtrema");
  }

  /* find P_core1, s.t. rest mass is m01 */
  vec[1] = 1e-10;   /* initial guess */
  /* do newton_linesrch_its iterations: */
  newton_linesrch_its(vec, 1, &check, m01_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  if(check) printf(": check=%d\n", check);  
  P_core1 = vec[1];
  printf(" setting: P_core1=%g\n", P_core1);

  /* TOV_init yields m01 for a given P_core1 */
  TOV_init(P_core1, kappa, Gamma, 1, &rf_surf1, &m1, &Phic1, &Psic1, &m01);
  printf(" rf_surf1=%g: m1=%g Phic1=%g Psic1=%g m01=%g\n",
         rf_surf1, m1, Phic1, Psic1, m01);

  if(Getd("BNSdata_m02")>0)
  {
    /* find P_core2, s.t. rest mass is m02 */
    vec[1] = 1e-10;   /* initial guess */
    /* do newton_linesrch_its iterations: */
    newton_linesrch_its(vec, 1, &check, m02_VectorFunc, 
                Geti("Coordinates_newtMAXITS"),
                Getd("Coordinates_newtTOLF") );
    if(check) printf(": check=%d\n", check);  
    P_core2 = vec[1];
    printf(" setting: P_core2=%g\n", P_core2);

    /* TOV_init yields m02 for a given P_core2 */
    TOV_init(P_core2, kappa, Gamma, 1, &rf_surf2, &m2, &Phic2, &Psic2, &m02);
    printf(" rf_surf2=%g: m2=%g Phic2=%g Psic2=%g m02=%g\n",
           rf_surf2, m2, Phic2, Psic2, m02);
  }
  else { P_core2=0.0; m2=m02=Phic2=0.0; Psic2=1.0; rf_surf2=rf_surf1; }

  /* find sigma1, s.t. radius is rf_surf1 */
  vec[1] = 1.0/rf_surf1;   /* initial guess */
  /* do newton_linesrch_its iterations: */
  newton_linesrch_its(vec, 1, &check, rf_surf1_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  if(check) printf(": check=%d\n", check);  
  sigma1 = vec[1];
  printf(" setting: sigma1=%g\n", sigma1);

  /* find sigma2, s.t. radius is rf_surf2 */
  vec[1] = -1.0/rf_surf2;   /* initial guess */
  /* do newton_linesrch_its iterations: */
  newton_linesrch_its(vec, 1, &check, rf_surf2_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  if(check) printf(": check=%d\n", check);  
  sigma2 = vec[1];
  printf(" setting: sigma2=%g\n", sigma2);

  /* find  xmin1,xmax1, xmin2,xmax2, xc1, xc2, and CM for both stars */
  xmax1 = x_of_AnsorgNS0(NULL, -1, 0.0,0.0,0.0);
  xmin1 = x_of_AnsorgNS0(NULL, -1, 0.0,1.0,0.0);
  xmax2 = x_of_AnsorgNS3(NULL, -1, 0.0,1.0,0.0);
  xmin2 = x_of_AnsorgNS3(NULL, -1, 0.0,0.0,0.0);
  printf(" => radius of domain0 = %g,   radius of domain3 = %g\n", 
         0.5*(xmax1-xmin1), 0.5*(xmax2-xmin2) );
  xc1 = 0.5*(xmax1+xmin1);
  xc2 = 0.5*(xmax2+xmin2);
  printf(" => domain0 extends from  xin1 = %g  to  xout1 = %g\n", xmin1,xmax1);
  printf(" => domain3 extends from  xout2 = %g  to  xin2 = %g\n", xmin2,xmax2);
  printf(" => centers are at: xc1 = %.13g,  xc2 = %.13g\n", xc1, xc2);
  printf(" BNSdata_b = %g\n", Getd("BNSdata_b"));
  DoM = fabs(xc1-xc2)/(m01+m02);
  DoM3 = DoM*DoM*DoM;
  DoM4 = DoM3*DoM;
  DoM5 = DoM4*DoM;
  nu = (m01*m02)/pow(m01+m02, 2.0);

  /* set: BNSdata_x_CM  = BNSdata_x_CM_init 
          BNSdata_Omega = BNSdata_Omega_init
     (unless: BNSdata_x_CM_init  = BNSdata_x_CM
              BNSdata_Omega_init = BNSdata_Omega) */
  if(!Getv("BNSdata_x_CM_init", "BNSdata_x_CM"))
    Sets("BNSdata_x_CM", Gets("BNSdata_x_CM_init"));
  if(!Getv("BNSdata_Omega_init", "BNSdata_Omega"))
    Sets("BNSdata_Omega", Gets("BNSdata_Omega_init"));
  printf(" initializing BNSdata_x_CM and BNSdata_Omega using:\n");
  printf(" BNSdata_x_CM  = %s\n", Gets("BNSdata_x_CM"));
  printf(" BNSdata_Omega = %s\n", Gets("BNSdata_Omega"));

  /* set CM and Omega (taken from PN_ADM_2.m) */
  if(Getv("BNSdata_x_CM", "estimate"))
    xCM = (m01*xc1 + m02*xc2)/(m01+m02);
  else
    xCM = Getd("BNSdata_x_CM");
  if(Getv("BNSdata_Omega", "estimate"))
    Omega = sqrt( 64*DoM3/pow(1 + 2*DoM, 6) +nu/DoM4 +
                 (-5*nu + 8*nu*nu)/(8*DoM5)            )/(m01+m02);
  else if(Getv("BNSdata_Omega", "estimate_from_desired_m0"))
  {
    double m1 = Getd("BNSdata_desired_m01");
    double m2 = Getd("BNSdata_desired_m02");
    DoM  = fabs(xc1-xc2)/(m1+m2);
    DoM3 = DoM*DoM*DoM;
    DoM4 = DoM3*DoM;
    DoM5 = DoM4*DoM;
    nu = (m1*m2)/pow(m1+m2, 2.0);
    Omega = sqrt( 64*DoM3/pow(1 + 2*DoM, 6) +nu/DoM4 +
                 (-5*nu + 8*nu*nu)/(8*DoM5)            )/(m01+m02);
  }
  else
    Omega = Getd("BNSdata_Omega");
  if(nu<=0.0) Omega=0.0;
  Setd("BNSdata_x_CM", xCM);
  Setd("BNSdata_Omega", Omega);
  printf(" BNSdata_x_CM and BNSdata_Omega are now set to:\n");
  printf(" BNSdata_x_CM = %g\n", Getd("BNSdata_x_CM"));
  printf(" BNSdata_Omega = %g,  (m01+m02)*BNSdata_Omega = %g\n",
         Getd("BNSdata_Omega"), Getd("BNSdata_Omega")*(m01+m02));

  /* set BNSdata_C1 */
  qc = pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
       pow(P_core1, 1.0/(1.0 + BNSdata_n));
  /* oouzerosqr == alpha2 - 
                   Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]), */
  oouzerosqr = exp(Phic1*2.0) - 
               pow(Psic1+m02/(2*fabs(xc1-xc2)), 4) *
               Omega*Omega*(xc1-xCM)*(xc1-xCM);
  Fc = -sqrt(oouzerosqr);
  /* q == (C1/F - 1.0)/(n+1.0) */
  Cc = Fc*((BNSdata_n+1.0)*qc + 1.0);
  Setd("BNSdata_C1", Cc);
  printf(" BNSdata_C1 = %g\n", Getd("BNSdata_C1"));

  /* set BNSdata_C2 */
  qc = pow(kappa, BNSdata_n/(1.0 + BNSdata_n)) *
       pow(P_core2, 1.0/(1.0 + BNSdata_n));
  /* oouzerosqr == alpha2 - 
                   Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]), */
  oouzerosqr = exp(Phic2*2.0) - 
               pow(Psic2+m01/(2*fabs(xc1-xc2)), 4) *
               Omega*Omega*(xc2-xCM)*(xc2-xCM);
  Fc = -sqrt(oouzerosqr);
  /* q == (C2/F - 1.0)/(n+1.0) */
  Cc = Fc*((BNSdata_n+1.0)*qc + 1.0);
  if(m02<=0.0) Cc = 0.0;
  Setd("BNSdata_C2", Cc);
  printf(" BNSdata_C2 = %g\n", Getd("BNSdata_C2"));

  /* set max A inside stars and adjust boxes4/5 accordingly */
  if(Getv("BNSdata_grid", "4ABphi_2xyz"))
  {
    double box0_max1 = Getd("BNSdata_box0_Amax");
    double box3_max1 = Getd("BNSdata_box3_Amax");
    double scal = 1.05; //1.551108723489246  /* make box4/5 5% larger than needed */
    double xr, xp, xm, xmax, xmin;
    double b = GetCachedNumValByParIndex(Coordinates_AnsorgNS_b_ParIndex);    

    Setd("box0_max1", box0_max1);
    xmin = x_of_AnsorgNS0(NULL, -1, box0_max1,1.0,0.0);
    xmax = x_of_AnsorgNS0(NULL, -1, box0_max1,0.0,0.0);
    xr = scal * 0.5*(xmax-xmin);
    xm = scal * (xmin-b);
    xp = scal * (xmax-b);
    Setd("box5_min1", b + xm);
    Setd("box5_max1", b + xp);
    Setd("box5_min2", -xr);
    Setd("box5_max2",  xr);
    Setd("box5_min3", -xr);
    Setd("box5_max3",  xr);

    Setd("box3_max1", box3_max1);
    xmin = x_of_AnsorgNS3(NULL, -1, box3_max1,0.0,0.0);
    xmax = x_of_AnsorgNS3(NULL, -1, box3_max1,1.0,0.0);
    xr = scal * 0.5*(xmax-xmin);
    xm = scal * (xmin+b);
    xp = scal * (xmax+b);
    Setd("box4_min1", -b + xm);
    Setd("box4_max1", -b + xp);
    Setd("box4_min2", -xr);
    Setd("box4_max2",  xr);
    Setd("box4_min3", -xr);
    Setd("box4_max3",  xr);
  }

  return 0;
}

/* funtion to be passed into newton_linesrch_its by ??? */
void rf_surf1_VectorFunc(int n, double *vec, double *fvec)
{
  double xmax, xmin;
  
  sigma1 = vec[1];
  xmax = x_of_AnsorgNS0(NULL, -1, 0.0,0.0,0.0);
  xmin = x_of_AnsorgNS0(NULL, -1, 0.0,1.0,0.0);
  fvec[1] = (xmax-xmin)-2.0*rf_surf1;
}

/* funtion to be passed into newton_linesrch_its by ??? */
void rf_surf2_VectorFunc(int n, double *vec, double *fvec)
{
  double xmax, xmin;

  sigma2 = vec[1];
  xmin = x_of_AnsorgNS3(NULL, -1, 0.0,0.0,0.0);
  xmax = x_of_AnsorgNS3(NULL, -1, 0.0,1.0,0.0);
  fvec[1] = (xmax-xmin)-2.0*rf_surf2;
}


/* setup initial Coordinates_AnsorgNS_sigma_pm vars */
int set_sigma_pm_vars(tGrid *grid)
{
  int sigma_pm, dsigma_pm_dB, dsigma_pm_dphi;
  int boxindex;

  printf("set_sigma_pm_vars: setting Coordinates_AnsorgNS_sigma_pm ...\n");

  /* make sure that Coordinates uses Coordinates_AnsorgNS_sigma_pm */
  if(!Getv("Coordinates_AnsorgNS_sigma_pm_vars", "yes"))
    errorexit("you need: Coordinates_AnsorgNS_sigma_pm_vars = yes");

  sigma_pm = Ind("Coordinates_AnsorgNS_sigma_pm");  
  dsigma_pm_dB = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  dsigma_pm_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  enablevar(grid, Ind("Temp1")); /* for coeffs in interpolation */

  /* set Coordinates_AnsorgNS_sigma_pm */
  for(boxindex=0; boxindex <=3; boxindex++)
  {
    tBox *box = grid->box[boxindex];
    int ind;

    enablevar_inbox(box, sigma_pm); 
    enablevar_inbox(box, dsigma_pm_dB);
    enablevar_inbox(box, dsigma_pm_dphi);
    forallpoints(box, ind)
    {
      if(boxindex<2) box->v[sigma_pm][ind] = sigma1;
      else           box->v[sigma_pm][ind] = sigma2;
    }
  }
  return 0;
}

/* funtion to be passed into newton_linesrch_its to find P_core1 from m01 */
void m01_VectorFunc(int n, double *vec, double *fvec)
{
  double kappa     = Getd("BNSdata_kappa");
  double BNSdata_n = Getd("BNSdata_n");
  double m01       = Getd("BNSdata_m01");
  double Gamma     = 1.0 + 1.0/BNSdata_n;
  double Pc, m, Phic, Psic, m0;

  Pc = vec[1];
  TOV_init(Pc, kappa, Gamma, 0, &rf_surf1, &m, &Phic, &Psic, &m0);
  fvec[1] = m0-m01;
}
/* funtion to be passed into newton_linesrch_its to find P_core2 from m02 */
void m02_VectorFunc(int n, double *vec, double *fvec)
{
  double kappa     = Getd("BNSdata_kappa");
  double BNSdata_n = Getd("BNSdata_n");
  double m02       = Getd("BNSdata_m02");
  double Gamma     = 1.0 + 1.0/BNSdata_n;
  double Pc, m, Phic, Psic, m0;

  Pc = vec[1];
  TOV_init(Pc, kappa, Gamma, 0, &rf_surf2, &m, &Phic, &Psic, &m0);
  fvec[1] = m0-m02;
}


/*************************************/
/* functions to adjust star surfaces */
/*************************************/

/* get q at A by interpolation along a (B,phi)=(const1,const2) line */
void BNdata_q_VectorFunc(int n, double *vec, double *fvec)
{
  tBox *box = BNSdata_q_VectorFunc_box;
  double *c = BNSdata_q_VectorFunc_coeffs;
  int i;
  double B1;
  double a1=box->bbox[0];
  double b1=box->bbox[1];
  int n1 = box->n1;
  double q = 0.0;
  double A = vec[1];  

  /* interpolate q onto point A */
  for(i=n1-1; i>=0; i--)
  {
    B1 = box->basis1((void *) box, a1,b1, i,n1, A);
    q += c[i] * B1;
  }   
  fvec[1] = q;
}

/* we use this to get sigp_Bphi in reset_Coordinates_AnsorgNS_sigma_pm__old1 */
void BNdata_sigp_VectorFunc(int n, double *vec, double *fvec)
{
  double sigp_Bphi = vec[1];
  double sigp_1phi = BNdata_sigp_VectorFunc_sigp_1phi;
  double B         = BNdata_sigp_VectorFunc_B;
  double X0        = BNdata_sigp_VectorFunc_X0;
  double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
  double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
  double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
  /* double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi); */
  double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
  double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
  double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
  /* double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi); */
  double X,R;

  /* use Eq. (22), (23) or (24) at A=0 */
  X = ReCp_Bphi - B*ReCp_1phi + B*cos(ArgCp_1phi);
  /* R = ImCp_Bphi - B*ImCp_1phi + B*sin(ArgCp_1phi); */

  fvec[1] = (X-X0);
}

// reset_Coordinates_AnsorgNS_sigma_pm__old1 DOES NOT WORK IF 0<B<1 :
/* reset sigma such that the zeros in BNSdata_q are at A=0 */
void reset_Coordinates_AnsorgNS_sigma_pm__old1(tGrid *grid, tGrid *gridnew,
                                         int innerdom,  int outerdom)
{
  int iq = Ind("BNSdata_q");
  int ic = Ind("Temp1"); /* we store coeffs of BNSdata_q in Temp1 */
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int isigma      = Ind("Coordinates_AnsorgNS_sigma_pm");
  int isigma_dB   = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  int isigma_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  double *q_in = grid->box[innerdom]->v[iq];
  double *q_out= grid->box[outerdom]->v[iq];  
  double *c_in = grid->box[innerdom]->v[ic];
  double *c_out= grid->box[outerdom]->v[ic];  
  int n1 = grid->box[innerdom]->n1;
  int n2 = grid->box[innerdom]->n2;
  int n3 = grid->box[innerdom]->n3;
  int i,j,k;
  int ineg_in;  /* q_in<0  at i=ineg_in (and q_in>=0 i=ineg_in+1) */
  int ipos_out; /* q_out>0 at i=ipos_out (and q_out<0 i=ipos_out+1) */
  int i1, i2, dom; /* zero occurs between index i1 and i2 in domain dom */
  double A1, A2;   /* zero occurs between A=A1 and A=A2 in domain dom */
  double A0;      /* q=0 at A=A0 in domain dom */
  double q1, q2;
  double X0, R0;  /* value of X,R at A=A0 */
  double B,phi, x,y,z;
  double ArgCp_1phi, AbsCp_1phi, ReCp_1phi, ImCp_1phi;
  double ArgCp_Bphi, ReCp_Bphi, ImCp_Bphi;
  double sigp_1phi, sigp_Bphi;
  int itmax = Geti("Coordinates_newtMAXITS");
  double tol = Getd("Coordinates_newtTOLF");
  double vec[2];
  int check;

  /* set coeffs of BNSdata_q */
  spec_analysis1(grid->box[innerdom], 1, q_in, c_in);
  spec_analysis1(grid->box[outerdom], 1, q_out, c_out);

  /* loop over j,k i.e. B,phi. 
     NOTE: we assue that n1,n2,n3 are the same in both domains */
  for(j=n2-1; j>=0; j--)
  {
    for(k=0; k<n3; k++)
    {
      /* find indices where q_in and q_out switch sign */
      for(i=n1-1; i>=0; i--) if(q_in[Index(i,j,k)]<0.0) break;
      ineg_in=i;
      for(i=n1-1; i>=0; i--) if(q_out[Index(i,j,k)]>0.0) break;
      ipos_out=i;

      /* if ipos_out=>0, q has zero in outer domain */
      if(ipos_out>=0) { i1=ipos_out; i2=ipos_out+1; dom=outerdom; }
      /* if ineg_in=>0, q has zero in inner domain */
      else if(ineg_in>=0) { i1=ineg_in; i2=ineg_in+1; dom=innerdom; }
           else           { i1=0;       i2=1;         dom=innerdom; }
      A1 = grid->box[dom]->v[iX][i1];
      A2 = grid->box[dom]->v[iX][i2];
      q1 = grid->box[dom]->v[iq][i1];
      q2 = grid->box[dom]->v[iq][i2];
      B   = grid->box[dom]->v[iY][Index(i1,j,k)];
      phi = grid->box[dom]->v[iZ][Index(i1,j,k)];

      /* find zero in q between A1 and A2 */
      if( fabs(q1) < tol || (ineg_in<0 && dom==innerdom) ) A0=A1;
      else if( fabs(q2) < tol ) A0=A2;
      else /* use root finder */
      {
        A0 = A1 - q1*(A2-A1)/(q2-q1); /* initial guess */
        /* use newton_linesrch_its to find A0 */
        BNSdata_q_VectorFunc_box    = grid->box[dom];
        BNSdata_q_VectorFunc_coeffs = grid->box[dom]->v[ic]+Index(0,j,k);
        vec[1] = A0;
        newton_linesrch_its(vec, 1, &check, BNdata_q_VectorFunc, itmax, tol);
        if(check)
          printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d\n", check);  
        A0 = vec[1];
        if(A0<A1 || A0>A2) errorexit("reset_Coordinates_AnsorgNS_sigma_pm: "
                                     "newton_linesrch_its failed!");
      }
printf("\nA0,B,phi=%g,%g,%g  dom=%d A1=%g A2=%g q1=%g q2=%g\n", A0,B,phi, dom, A1,A2, q1,q2);

      /* compute values of X0,R0 at A=A0 */
      xyz_of_AnsorgNS(grid->box[dom], -1, dom, A0,B,phi, &x,&y,&z, &X0,&R0);
printf("old X0=%g R0=%g\n", X0,R0);      
      /* get Cp and sigp at B=1  */
      if(j==n2-1 && k==0) /* B=1 case, but compute it only for phi=0 */
      {
        ArgCp_1phi = acos(X0);
        /* 2 ArgCp = ArcTan[Sin[Pi B/2]/Sinh[sigma/2]] */
        /* Tan[2 ArgCp] = Sin[Pi B/2]/Sinh[sigma/2] */
        sigp_1phi = 2.0 * asinh( (1.0/tan(2.0*ArgCp_1phi)) );
        AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
        ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
        ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
        sigp_Bphi = sigp_1phi;
printf("ArgCp_1phi=%g ReCp_1phi=%g ImCp_1phi=%g\n",ArgCp_1phi, ReCp_1phi,ImCp_1phi);
      }
      if( (j==0 && k==0) || j>0 && j<n2-1 ) /* B<1 case */
      {                              /* if B=0 compute only for phi=0 */
        ReCp_Bphi = X0 + B*(ReCp_1phi - cos(ArgCp_1phi));
        ImCp_Bphi = R0 + B*(ImCp_1phi - sin(ArgCp_1phi));
        ArgCp_Bphi = Arg(ReCp_Bphi, ImCp_Bphi);
        if(B==0) /* Cp_Bphi^2 = tanh(0.25*sigp_Bphi) */
          sigp_Bphi = 4.0 * atanh(ReCp_Bphi*ReCp_Bphi-ImCp_Bphi*ImCp_Bphi);
        else
        {
          sigp_Bphi = 2.0 * asinh( sin(PIh*B)/tan(2.0*ArgCp_Bphi) );
          //or
          /* Note: |tanh(x+iy)|^2 = K  =>  2x =acosh( cos(2y)*(1+K)/(1-K) ) */
          //double AbsCpSqr_Bphi=ReCp_Bphi*ReCp_Bphi+ImCp_Bphi*ImCp_Bphi;
          //double K=AbsCpSqr_Bphi*AbsCpSqr_Bphi;
          //sigp_Bphi = 2.0 * acosh( cos(PIh*B)*(1.0+K)/(1.0-K) );
          //sigp_Bphi = sigp_Bphi*0.5 + acosh( cos(PIh*B)*(1.0+K)/(1.0-K) ) ;
        }
        /* use newton_linesrch_its to find better sigp_Bphi */
        BNdata_sigp_VectorFunc_sigp_1phi = sigp_1phi;
        BNdata_sigp_VectorFunc_B         = B;
        BNdata_sigp_VectorFunc_X0        = X0;
        vec[1] = sigp_Bphi;
        newton_linesrch_its(vec, 1, &check, BNdata_sigp_VectorFunc, itmax, tol);
        if(check)
          printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d\n", check);  
        sigp_Bphi = vec[1];
printf("ReCp_Bphi*ReCp_Bphi+ImCp_Bphi*ImCp_Bphi=%g ", 
sqrt(ReCp_Bphi*ReCp_Bphi+ImCp_Bphi*ImCp_Bphi));
printf("ArgCp_Bphi=%g ReCp_Bphi=%g ImCp_Bphi=%g\n",ArgCp_Bphi, ReCp_Bphi,ImCp_Bphi);
{
double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);

double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
//double Ap=0;
double X = (ReCp_Bphi - B*ReCp_1phi) + B*cos(ArgCp_1phi);
double R = (ImCp_Bphi - B*ImCp_1phi) + B*sin(ArgCp_1phi);
printf("new X=%g R=%g  ArgCp_Bphi=%g AbsCp_Bphi=%g ReCp_Bphi=%g ImCp_Bphi=%g\n", 
X,R, ArgCp_Bphi, AbsCp_Bphi, ReCp_Bphi,ImCp_Bphi); 
}
      }
printf("sigp_Bphi=%g sigp_1phi=%g\n", sigp_Bphi, sigp_1phi);

      /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains */
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
      }
    } /* end for k */
  } /* end for j */
  
  /* compute derivs of sigma */
  spec_Deriv1(gridnew->box[innerdom], 2, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[outerdom], 3, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dphi]);
}


/* find diff. {Xc(A,B,phi),Rc(A,B,phi)}-{X, R} for any A,B,phi*/
void DelXR_of_AB_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double x,y,z, Xc,Rc;
  double A = vec[1];
  double B = vec[2];
  t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *pars;
  tBox *box;
  double phi, X, R;

  /* get pars */
  pars = (t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *) p;
  phi = pars->phi;
  box = pars->box;
  X   = pars->X;
  R   = pars->R;

  /* check if we are within the ranges of A and B */
  if(A>=0.0 && A<=1.0 && B>=0 && B<=1)
  {
    /* get diff between {Xc(A,B,phi),Rc(A,B,phi)}-{X, R} */
    xyz_of_AnsorgNS(box, -1, box->b, A,B,phi, &x,&y,&z, &Xc,&Rc);
    fvec[1] = Xc-X;
    fvec[2] = Rc-R;
  }
  else /* out of range: get fvec by interp. */
  {
    double hi[3];
    double lo[3];
    double vb[3];
    double vi[3];
    double fvb[3];
    double fvi[3];

    lo[1] = lo[2] = 0.0;
    hi[1] = hi[2] = 1.0;

    /* get vb, vi */
    newton_lnsrch_set_vecs_for_lininterp(2, vec,hi,lo, vb, vi);

    /* compute fvb,fvi */
    xyz_of_AnsorgNS(box, -1, box->b, vb[1],vb[2],phi, &x,&y,&z, &Xc,&Rc);
    fvb[1] = Xc-X;
    fvb[2] = Rc-R;
    xyz_of_AnsorgNS(box, -1, box->b, vi[1],vi[2],phi, &x,&y,&z, &Xc,&Rc);
    fvi[1] = Xc-X;
    fvi[2] = Rc-R;
            
    /* set fvec from fvb,fvi */
    newton_lnsrch_get_fvec_by_lininterp(2, vec,vb,vi, fvec,fvb,fvi);

    /* print */
    /*
    printf("DelXR_of_AB_VectorFuncP: box->b=%d, vec=(A,B)=(%g,%g)\n", 
           box->b, A,B);
    printf("vec=(%g,%g) fvec=(%g,%g)\n", vec[1],vec[2], fvec[1],fvec[2]);
    printf(" vb=(%g,%g)  fvb=(%g,%g)\n", vb[1],vb[2], fvb[1],fvb[2]);
    printf(" vi=(%g,%g)  fvi=(%g,%g)\n", vi[1],vi[2], fvi[1],fvi[2]);
    */
  }
}

/* find diff. Xc(A,B,phi)-X for any A 
   calls DelXR_of_AB_VectorFunc above and uses the same pars,
   but ensures that B=0 */
void DelXR_of_A_forB0_VectorFuncP(int n, double *vec1, double *fvec1, void *p)
{
  tBox *box;
  double vec[3];
  double fvec[3];
  t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *pars;

  /* get pars */
  pars = (t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *) p;
  box = pars->box;
  
  vec[1] = vec1[1];
  vec[2] = 0.0; /* since B=0 */
  DelXR_of_AB_VectorFuncP(2, vec, fvec, p);
  if(box->b<2) fvec1[1] = fvec[1]; /* in box0/1 R=0 at B=0 */
  else         fvec1[1] = fvec[2]; /* in box3/2 R=0 at B=0 */
}

/* WE NEED to find sigp_Bphi at B,phi such that q(sigp_Bphi; A=0, B, phi)=0 */
/* q as a func of sigp for a given A=0, B, phi */
void q_of_sigp_forgiven_Bphi__old2(int n, double *sigvec, double *qvec)
{
  double sigp_Bphi = sigvec[1];
  double sigp_1phi = q_of_sigp_forgiven_Bphi__sigp_1phi;
  double B         = q_of_sigp_forgiven_Bphi__B;
  double phi       = q_of_sigp_forgiven_Bphi__phi;
  tGrid *grid      = q_of_sigp_forgiven_Bphi__grid;
  int icoeffs      = q_of_sigp_forgiven_Bphi__icoeffs;
  int innerdom     = q_of_sigp_forgiven_Bphi__innerdom;
  int outerdom     = q_of_sigp_forgiven_Bphi__outerdom;
  int iq = Ind("BNSdata_q");
  double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
  double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
  double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
  double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
  double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
  double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
  double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
  double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
  double X,R;
  double Ac,Bc, Acin,Bcin, Acout,Bcout, Acmax, q;
  double vec[3];
  int i, check, stat,statin,statout, dom;

  /* use Eq. (22), (23) or (24) at A=0 to compute X,R */
  X = ReCp_Bphi - B*ReCp_1phi + B*cos(ArgCp_1phi);
  R = ImCp_Bphi - B*ImCp_1phi + B*sin(ArgCp_1phi);

  /* set Acin,Bcin, Acout,Bcout, statin,statout to invalid values */
  statin=statout=-1;
  Acin=Acout=Bcin=Bcout=-1.0;

  /* find domain and Ac,Bc on current grid, corresponding to X,R */
  dom = innerdom;
  for(i=1; i<=2; i++)
  {
    t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];
    pars->phi = phi;
    pars->box = grid->box[dom];
    pars->X   = X;
    pars->R   = R;
    Acmax  = grid->box[dom]->bbox[1];
    vec[1] = 1e-7; /* initial guess is that Ac,Bc = 0,B*/
    vec[2] = B;
    if(dequal(B,0.0))
      stat = newton_linesrch_itsP(vec, 1, &check, DelXR_of_A_forB0_VectorFuncP, 
                                  (void *) pars, 1000, 1e-10);
    else
      stat = newton_linesrch_itsP(vec, 2, &check, DelXR_of_AB_VectorFuncP, 
                                  (void *) pars, 1000, 1e-10);
    if(check) printf("q_of_sigp_forgiven_Bphi: check=%d\n", check);  
    Ac = vec[1];
    Bc = vec[2];

    /* save vals for later */
    if(dom == innerdom) {  Acin=Ac;  Bcin=Bc;  statin=stat; }
    else                { Acout=Ac; Bcout=Bc; statout=stat; }
//printf("  sigp_Bphi=%g sigp_1phi=%g:\n"
//"       X=%g R=%g: stat=%d dom=%d Ac=%g Bc=%g\n",
//sigp_Bphi,sigp_1phi, X,R, stat,dom, Ac,Bc);
    /* if(stat>=0 && Ac>=0.0 && Ac<=Acmax && Bc>=0.0 && Bc<=1.0) break; */
    if(stat>=0 && dlesseq(0.0,Ac) && dlesseq(Ac,1.0) &&
                  dlesseq(0.0,Bc) && dlesseq(Bc,1.0)   ) break;
    dom = outerdom;
  }
  /* decide which results to use */
  if(dom == outerdom && statin>=0)
  {
    double dA=fabs(Acout)-fabs(Acin);
    double dB=fabs(Bcout)-fabs(Bcin);
    /* switch back to innerdom in some cases */
    if(Ac<0.0)
      if( dA>0.0 && dlesseq(0.0,Bcin) && dlesseq(Bcin,1.0) )
      {dom=innerdom; Ac=0.0; Bc=Bcin; stat=statin;}
    if(Bc<0.0)
      if( dB>0.0 && dlesseq(0.0,Acin) && dlesseq(Acin,1.0) )
      {dom=innerdom; Ac=Acin; Bc=0.0; stat=statin;}
    if(Bc>1.0)
      if(fabs(Bcin)<fabs(Bcout) && dlesseq(0.0,Acin) && dlesseq(Acin,1.0))
      {dom=innerdom; Ac=Acin; Bc=1.0; stat=statin;}
  }
  /* check for failure */
  if(stat<0 || dless(Ac,0.0) || dless(1.0,Ac) ||
               dless(Bc,0.0) || dless(1.0,Bc)   )
  {
    printf("q_of_sigp_forgiven_Bphi: stat=%d dom=%d Ac=%g Bc=%g\n",
           stat,dom, Ac,Bc);
    printf("statin=%d Acin=%g Bcin=%g  statout=%d Acout=%g Bcout=%g\n",
           statin, Acin,Bcin, statout, Acout,Bcout);
    printf("q_of_sigp_forgiven_Bphi: X=%g R=%g for: A=0 B=%g phi=%g\n"
           " sigp_Bphi=%g sigp_1phi=%g\n"
           " ReCp_Bphi=%g ImCp_Bphi=%g ReCp_1phi=%g ImCp_1phi=%g\n",
           X,R, B,phi, sigp_Bphi,sigp_1phi,
           ReCp_Bphi,ImCp_Bphi, ReCp_1phi,ImCp_1phi);
    errorexit("q_of_sigp_forgiven_Bphi: could not find Ac,Bc");
  }
  if(Ac<0.0) Ac=0.0; /* make sure we stay in our box */
  if(Ac>1.0) Ac=1.0;
  if(Bc<0.0) Bc=0.0;
  if(Bc>1.0) Bc=1.0;

  /* if we get point in innerdom with Acmax<Ac<1 we retrun a huge value
     for q, so that newton_linesrch_its thinks it has to backtrack... */
  if(Ac>Acmax && dom==0)      q = 1000.0*grid->box[5]->v[iq][0];
  else if(Ac>Acmax && dom==3) q = 1000.0*grid->box[4]->v[iq][0];
  else
    /* obtain q at Ac,Bc,phi by interpolation */
    /* grid->box[dom]->v[icoeffs]  contains coeffs of q in box */
    q = spec_interpolate(grid->box[dom], grid->box[dom]->v[icoeffs],Ac,Bc,phi);
//printf("         dom=%d Ac=%g Bc=%g  sigp_Bphi=%g q=%g\n",dom,Ac,Bc, sigp_Bphi,q);
  qvec[1] = q;
//qvec[1] = sigp_Bphi-0.9;
}


/* reset sigma such that the zeros in BNSdata_q are at A=0 */
void reset_Coordinates_AnsorgNS_sigma_pm__old2(tGrid *grid, tGrid *gridnew,
                                         int innerdom,  int outerdom)
{
  int iq = Ind("BNSdata_q");
  int ic = Ind("BNSdata_temp1"); /* we store 1D coeffs of BNSdata_q in BNSdata_temp1 */
  int ic3= Ind("BNSdata_temp2"); /* we store 3D coeffs of BNSdata_q in BNSdata_temp2 */
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int isigma      = Ind("Coordinates_AnsorgNS_sigma_pm");
  int isigma_dB   = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  int isigma_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  double *q_in = grid->box[innerdom]->v[iq];
  double *q_out= grid->box[outerdom]->v[iq];  
  double *c_in = grid->box[innerdom]->v[ic];
  double *c_out= grid->box[outerdom]->v[ic];  
  double *c3_in = grid->box[innerdom]->v[ic3];
  double *c3_out= grid->box[outerdom]->v[ic3];  
  int n1 = grid->box[innerdom]->n1;
  int n2 = grid->box[innerdom]->n2;
  int n3 = grid->box[innerdom]->n3;
  int i,j,k, kk;
  int inz_in;   /* q_in<=0  at i=inz_in (and q_in>0 i=inz_in+1) */
  int inz_out;  /* q_out<=0 at i=inz_out (and q_out>0 i=inz_out-1) */
  int i1, i2, dom; /* zero occurs between index i1 and i2 in domain dom */
  double A1, A2;   /* zero occurs between A=A1 and A=A2 in domain dom */
  double A0;      /* q=0 at A=A0 in domain dom */
  double q1, q2;
  double X0, R0;  /* value of X,R at A=A0 */
  double B,phi, x,y,z;
  double ArgCp_1phi, AbsCp_1phi, ReCp_1phi, ImCp_1phi;
  double ArgCp_Bphi, ReCp_Bphi, ImCp_Bphi;
  double sigp_1phi, sigp_0phi, sigp_Bphi;
  int itmax = Geti("Coordinates_newtMAXITS");
  double tol = Getd("Coordinates_newtTOLF");
  double vec[2];
  int check, stat;

  /* set 1D-coeffs of BNSdata_q for interpolation in A-direc. */
  spec_analysis1(grid->box[innerdom], 1, q_in, c_in);
  spec_analysis1(grid->box[outerdom], 1, q_out, c_out);

  /* set 3D-coeffs of BNSdata_q for 3D interpolation */
  spec_Coeffs(grid->box[innerdom], q_in, c3_in);
  spec_Coeffs(grid->box[outerdom], q_out, c3_out);

  /* look at B=1 (j=n2-1) and B=0 (j=0)  
     NOTE: we assue that n1,n2,n3 are the same in both domains */
  for(k=0, j=n2-1; j>=0; j-=n2-1)
  {
    /* find indices where q_in and q_out switch sign */
    for(i=n1-1; i>=0; i--) if(q_in[Index(i,j,k)]<=0.0) break;
    inz_in=i;
    for(i=0; i<n1; i++)   if(q_out[Index(i,j,k)]<=0.0) break;
    inz_out=i;

    /* if inz_in=>0, q has zero in inner domain */
    /* if inz_out<n1, q is negative in outer domain */
    if(inz_in>=0)                   { i1=inz_in;   i2=inz_in+1; dom=innerdom;}
    else if(inz_out==0)             { i1=0;         i2=1;       dom=innerdom;}
    else if(inz_out<n1 && inz_out>0){ i1=inz_out-1; i2=inz_out; dom=outerdom;}
    else
    {
      printf("reset_Coordinates_AnsorgNS_sigma_pm: innerdom=%d  B=%g phi=%g  "
             "inz_in=%d inz_out=%d\n", innerdom, B,phi, inz_in,inz_out);
      printf("q_in[Index(n1-1,j,k)]=%g\n", q_in[Index(n1-1,j,k)]);
      errorexit("reset_Coordinates_AnsorgNS_sigma_pm: q>0 everywhere???");
    }
    A1 = grid->box[dom]->v[iX][Index(i1,j,k)];
    A2 = grid->box[dom]->v[iX][Index(i2,j,k)];
    q1 = grid->box[dom]->v[iq][Index(i1,j,k)];
    q2 = grid->box[dom]->v[iq][Index(i2,j,k)];
    B   = grid->box[dom]->v[iY][Index(i1,j,k)];
    phi = grid->box[dom]->v[iZ][Index(i1,j,k)];
//printf("inz_out=%d inz_in=%d\n", inz_out,inz_in);
//printf("reset_Coordinates_AnsorgNS_sigma_pm: "
//"find zero in q in dom=%d between A1=%g A2=%g\n", dom, A1,A2);

    /* find zero in q between A1 and A2 */
    if( fabs(q1) < tol || (inz_in<0 && dom==innerdom) ) A0=A1;
    else if(q1*q2>=0) A0=0.0;
    else /* use root finder */
    {
      A0 = A1 - q1*(A2-A1)/(q2-q1); /* initial guess */
      /* use newton_linesrch_its to find A0 */
      BNSdata_q_VectorFunc_box    = grid->box[dom];
      BNSdata_q_VectorFunc_coeffs = grid->box[dom]->v[ic]+Index(0,j,k);
      vec[1] = A0;
      stat=newton_linesrch_its(vec, 1, &check, BNdata_q_VectorFunc, itmax, tol);
      if(check || stat<0)
        printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d stat=%d\n",
               check, stat);
      A0 = vec[1];
//printf("stat=%d: zero in q at A=A0=%g\n", stat, A0); Yo(1);
      /* go back to initial guess if we leave A1<A0<A2 in outerdom */
      if( (A0 < A1 || A0 > A2) && dom==outerdom )
        A0 = A1 - q1*(A2-A1)/(q2-q1);
      if(A0 < A1 - 0.05*(A2-A1) || A0 > A2 + 0.05*(A2-A1))
      {
        printf("reset_Coordinates_AnsorgNS_sigma_pm: B=%g phi=%g  "
               "inz_in=%d inz_out=%d\n", B,phi, inz_in,inz_out);
        printf("looked for zero in q in dom=%d between A1=%g A2=%g\n",
               dom, A1,A2);
        printf("q1=q(A1)=%g  q2=q(A2)=%g\n", q1,q2);
        printf("stat=%d: zero in q at A=A0=%g\n", stat, A0);        
        errorexit("reset_Coordinates_AnsorgNS_sigma_pm: "
                  "newton_linesrch_its failed!");
      }
      if(A0<0.0) A0=0.0; /* make sure we stay in box */
      if(A0>1.0) A0=1.0;
    }
//printf("zero in q at A=A0=%g\n", A0);

    /* compute values of X0,R0 at A=A0 */
    xyz_of_AnsorgNS(grid->box[dom], -1, dom, A0,B,phi, &x,&y,&z, &X0,&R0);

    /* get Cp and sigp at B=1  */
    if(j==n2-1) /* B=1 case */
    {
      ArgCp_1phi = acos(X0);
      /* 2 ArgCp = ArcTan[Sin[Pi B/2]/Sinh[sigma/2]] */
      /* Tan[2 ArgCp] = Sin[Pi B/2]/Sinh[sigma/2] */
      sigp_1phi = 2.0 * asinh( (1.0/tan(2.0*ArgCp_1phi)) );
      AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
      ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
      ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
      sigp_Bphi = sigp_1phi;
    }
    if(j==0) /* B=0 case */
    {
      ReCp_Bphi = X0;
      ImCp_Bphi = R0;
      ArgCp_Bphi = Arg(ReCp_Bphi, ImCp_Bphi);
      /* Cp_Bphi^2 = tanh(0.25*sigp_Bphi) */
      sigp_0phi = 4.0 * atanh(ReCp_Bphi*ReCp_Bphi-ImCp_Bphi*ImCp_Bphi);
      sigp_Bphi = sigp_0phi;
    }

    /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains 
       at B=1 and B=0 */
    for(kk=0; kk<n3; kk++)
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,kk)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,kk)] = sigp_Bphi;
      }
  } /* end for j */
//printf("reset_Coordinates_AnsorgNS_sigma_pm: new "
//       "sigp_0phi=%g sigp_1phi=%g\n", sigp_0phi, sigp_1phi);

  /* guess for sigp_Bphi at j=n2-2 */
  sigp_Bphi = sigp_1phi;
  
  /* loop over the remaining j,k i.e. B,phi. 
     NOTE: we assume that n1,n2,n3 are the same in both domains */
  for(j=n2-2; j>0; j--) /* we could include j=0 (B=0) here again, so that most sigp_Bphi are found with the same method */
    for(k=0; k<n3; k++)
    {
      /* find sigp_Bphi at B,phi such that q(sigp_Bphi; A=0, B, phi)=0 */
      B   = grid->box[dom]->v[iY][Index(0,j,k)];
      phi = grid->box[dom]->v[iZ][Index(0,j,k)];
      /* use newton_linesrch_its to find sigp_Bphi */
      q_of_sigp_forgiven_Bphi__sigp_1phi = sigp_1phi;
      q_of_sigp_forgiven_Bphi__B = B;
      q_of_sigp_forgiven_Bphi__phi = phi;
      q_of_sigp_forgiven_Bphi__grid = grid;
      q_of_sigp_forgiven_Bphi__icoeffs = ic3;
      q_of_sigp_forgiven_Bphi__innerdom = innerdom;
      q_of_sigp_forgiven_Bphi__outerdom = outerdom;
      vec[1] = sigp_Bphi;
//vec[1] = grid->box[dom]->v[isigma][Index(0,j,k)];
//printf("itmax=%d tol=%g vec[1]=%g B=%g phi=%g\n",itmax,tol,vec[1], B,phi);
      stat=newton_linesrch_its(vec, 1, &check,
                               q_of_sigp_forgiven_Bphi__old2, itmax, tol);
      /* If q is nowhere negative newton_linesrch_its may not work. In this
         case we should probably search for the zero in (q - 1e-8). */
//printf("stat=%d\n",stat);
      if(check)
        printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d\n", check);  
      sigp_Bphi = vec[1];

      /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains */
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
      }
//printf("B=%g phi=%g  ", B, phi);
//printf("sigp_Bphi=%g sigp_0phi=%g sigp_1phi=%g\n", sigp_Bphi, sigp_0phi, sigp_1phi);
    } /* end for j,k */

  /* make sure that sigma has only one value at B=0 and also at B=1 */
  for(j=0; j<n2; j+=n2-1)
    for(k=1; k<n3; k++)
      for(i=0; i<n1; i++)
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] =
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] =
                        gridnew->box[innerdom]->v[isigma][Index(0,j,0)];

  /* compute derivs of sigma */
  spec_Deriv1(gridnew->box[innerdom], 2, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[outerdom], 3, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dphi]);
}


/* get q at A by direct computation along a (B,phi)=(const1,const2) line */
void BNSdata_q_VectorFunc(int n, double *vec, double *fvec)
{
  tBox *box = BNSdata_q_VectorFunc_box;
  tGrid *grid = box->grid;
  int b = box->b;
  double B   = BNSdata_q_VectorFunc_B;
  double phi = BNSdata_q_VectorFunc_phi;
  double A = vec[1];  

  /* compute q */
  fvec[1] = BNS_compute_new_centered_q_atXYZ(grid,b, A,B,phi);
}

/* go over box and find (A,B) that is closest to (X,R) */
double find_nearest_A_B_given_X_R_phi(tBox *box,
                                      double X, double R, double phi,
                                      double *A, double *B)
{
  int Xind = Ind("X");
  int Yind = Ind("Y");
  int i,j,k, pl;
  double min, dist;
  double vec[3], fvec[3];
  t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];

  if(box->b>3) errorexit("find_nearest_A_B_given_X_R_phi: "
                         "box->b must be between 0 and 3.");

  pars->phi = phi;
  pars->box = box;
  pars->X   = X;
  pars->R   = R;
  pl = find_ind_closest_to_Z0(box, phi);

  min = 1e300;
  forplane3(i,j,k,  box->n1,box->n2,box->n3, pl)
  {
    vec[1] = box->v[Xind][i];
    vec[2] = box->v[Yind][i];
    DelXR_of_AB_VectorFuncP(2, vec, fvec, (void *) pars);
    dist = sqrt(fvec[1]*fvec[1] + fvec[2]*fvec[2]);
    if(dist<min) 
    { dist = min;   *A=vec[1];   *B=vec[2]; }
  }
  return dist;
}

/* WE NEED to find sigp_Bphi at B,phi such that q(sigp_Bphi; A=0, B, phi)=0 */
/* q as a func of sigp for a given A=0, B, phi */
void q_of_sigp_forgiven_BphiP(int n, double *sigvec, double *qvec, void *p)
{
  double sigp_Bphi = sigvec[1];
  /* get pars */
  t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *pars = 
         (t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *) p;
  double sigp_1phi = pars->sigp_1phi;
  double B         = pars->B;
  double phi       = pars->phi;
  tGrid *grid      = pars->grid;
  int innerdom     = pars->innerdom;
  int outerdom     = pars->outerdom;
  double AbsCp_Bphi, ArgCp_Bphi, ReCp_Bphi, ImCp_Bphi;
  double AbsCp_1phi, ArgCp_1phi, ReCp_1phi, ImCp_1phi;
  double X,R;
  double Ac,Bc, Acin,Bcin, Acout,Bcout, Acmax, q;
  double vec[3];
  int i, check, stat,statin,statout, dom;
  int guessmode;

  /* set Cp at B,phi and 1,phi */
  AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
  ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
  ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
  ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
  AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
  ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
  ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
  ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);

  /* use Eq. (22), (23) or (24) at A=0 to compute X,R */
  X = ReCp_Bphi - B*ReCp_1phi + B*cos(ArgCp_1phi);
  R = ImCp_Bphi - B*ImCp_1phi + B*sin(ArgCp_1phi);

  /* try 2 ways of finding an initial guess for Ac,Bc */
  for(guessmode=0; guessmode<=3; guessmode++)
  {
    /* set Acin,Bcin, Acout,Bcout, statin,statout to invalid values */
    statin=statout=-1;
    Acin=Acout=Bcin=Bcout=-1.0;
  
    /* find domain and Ac,Bc on current grid, corresponding to X,R */
    dom = innerdom;
    for(i=1; i<=2; i++)
    {
      t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];
      pars->phi = phi;
      pars->box = grid->box[dom];
      pars->X   = X;
      pars->R   = R;
      Acmax  = grid->box[dom]->bbox[1];
  
      /* get initial guess for Ac,Bc in vec[1],vec[2] */
      if(guessmode==0)
      {
        vec[1] = 1e-7; /* initial guess is that Ac,Bc = 0,B */
        vec[2] = B;
      }
      else if(guessmode==1)
      {
        find_nearest_A_B_given_X_R_phi(grid->box[dom], X,R,phi, &Ac,&Bc);
        if(dlesseq(Ac,0.0)) Ac=1e-7;
        if(dequal(Bc,0.0) || dequal(Bc,1.0)) Bc=B;
        if(dequal(B,0.0))   Bc=0.0;
        vec[1] = Ac;
        vec[2] = Bc;
      }
      else if(guessmode==2)
      {
        vec[1] = 0.01;  /* initial guess is that Ac,Bc = 0.01,B */
        vec[2] = B;
      }
      else
      {
        vec[1] = 0.1;   /* initial guess is that Ac,Bc = 0.1,B */
        vec[2] = B;
      }
//printf("q_of_sigp_forgiven_BphiP: guessmode=%d (A,B)=(%g,%g)\n",
//guessmode, vec[1],vec[2]);
      /* do newton line searches */
      if(dequal(B,0.0))
        stat = newton_linesrch_itsP(vec, 1, &check, DelXR_of_A_forB0_VectorFuncP,
                                    (void *) pars, 1000, 1e-10);
      else
        stat = newton_linesrch_itsP(vec, 2, &check, DelXR_of_AB_VectorFuncP,
                                    (void *) pars, 1000, 1e-10);
      if(check) printf("q_of_sigp_forgiven_BphiP: check=%d\n", check);  
      Ac = vec[1];
      Bc = vec[2];
  
      /* save vals for later */
      if(dom == innerdom) {  Acin=Ac;  Bcin=Bc;  statin=stat; }
      else                { Acout=Ac; Bcout=Bc; statout=stat; }
//printf("  sigp_Bphi=%g sigp_1phi=%g:\n"
//"       X=%g R=%g: stat=%d dom=%d Ac=%g Bc=%g\n",
//sigp_Bphi,sigp_1phi, X,R, stat,dom, Ac,Bc);
      /* if(stat>=0 && Ac>=0.0 && Ac<=Acmax && Bc>=0.0 && Bc<=1.0) break; */
      if(stat>=0 && dlesseq(0.0,Ac) && dlesseq(Ac,1.0) &&
                    dlesseq(0.0,Bc) && dlesseq(Bc,1.0)   ) break;
      dom = outerdom;
    }
    /* decide which results to use */
    if(dom == outerdom && statin>=0)
    {
      double dA=fabs(Acout)-fabs(Acin);
      double dB=fabs(Bcout)-fabs(Bcin);
      /* switch back to innerdom in some cases */
      if(Ac<0.0)
        if( dA>0.0 && dlesseq(0.0,Bcin) && dlesseq(Bcin,1.0) )
        {dom=innerdom; Ac=0.0; Bc=Bcin; stat=statin;}
      if(Bc<0.0)
        if( dB>0.0 && dlesseq(0.0,Acin) && dlesseq(Acin,1.0) )
        {dom=innerdom; Ac=Acin; Bc=0.0; stat=statin;}
      if(Bc>1.0)
        if(fabs(Bcin)<fabs(Bcout) && dlesseq(0.0,Acin) && dlesseq(Acin,1.0))
        {dom=innerdom; Ac=Acin; Bc=1.0; stat=statin;}
    }
    /* leave guessmode loop if we found a resonable result */
    if(stat>=0 && dgreatereq(Ac,0.0) && dgreatereq(1.0,Ac) &&
                  dgreatereq(Bc,0.0) && dgreatereq(1.0,Bc)   )  break;
  }
  /* check for failure */
  if(stat<0 || dless(Ac,0.0) || dless(1.0,Ac) ||
               dless(Bc,0.0) || dless(1.0,Bc)   )
  {
    printf("q_of_sigp_forgiven_BphiP: guessmode=%d\n", guessmode);
    printf("q_of_sigp_forgiven_BphiP: stat=%d dom=%d Ac=%g Bc=%g\n",
           stat,dom, Ac,Bc);
    printf("statin=%d Acin=%g Bcin=%g  statout=%d Acout=%g Bcout=%g\n",
           statin, Acin,Bcin, statout, Acout,Bcout);
    printf("q_of_sigp_forgiven_BphiP: X=%g R=%g for: A=0 B=%g phi=%g\n"
           " sigp_Bphi=%g sigp_1phi=%g\n"
           " ReCp_Bphi=%g ImCp_Bphi=%g ReCp_1phi=%g ImCp_1phi=%g\n",
           X,R, B,phi, sigp_Bphi,sigp_1phi,
           ReCp_Bphi,ImCp_Bphi, ReCp_1phi,ImCp_1phi);
    errorexit("q_of_sigp_forgiven_BphiP: could not find Ac,Bc");
  }
  if(Ac<0.0) Ac=0.0; /* make sure we stay in our box */
  if(Ac>1.0) Ac=1.0;
  if(Bc<0.0) Bc=0.0;
  if(Bc>1.0) Bc=1.0;

  /* if we get point in innerdom with Acmax<Ac<1 we retrun a huge value
     for q, so that newton_linesrch_its thinks it has to backtrack... */
  if(Ac>Acmax && dom==0)      q = 1e6;
  else if(Ac>Acmax && dom==3) q = 1e6;
  else
    /* obtain q at Ac,Bc,phi by direct computation */
    /* grid->box[dom]->v[icoeffs]  contains coeffs of q in box */
    q = BNS_compute_new_centered_q_atXYZ(grid,dom, Ac,Bc,phi);
  qvec[1] = q;
}


/* reset sigma such that the zeros in BNSdata_q are at A=0 */
void reset_Coordinates_AnsorgNS_sigma_pm(tGrid *grid, tGrid *gridnew,
                                         int innerdom,  int outerdom)
{
  int iq = Ind("BNSdata_q");
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int isigma      = Ind("Coordinates_AnsorgNS_sigma_pm");
  int isigma_dB   = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  int isigma_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  double *q_in = grid->box[innerdom]->v[iq];
  double *q_out= grid->box[outerdom]->v[iq];  
  int n1 = grid->box[innerdom]->n1;
  int n2 = grid->box[innerdom]->n2;
  int n3 = grid->box[innerdom]->n3;
  int i,j,k, kk, JK;
  int inz_in;   /* q_in<=0  at i=inz_in (and q_in>0 i=inz_in+1) */
  int inz_out;  /* q_out<=0 at i=inz_out (and q_out>0 i=inz_out-1) */
  int i1, i2, dom; /* zero occurs between index i1 and i2 in domain dom */
  double A1, A2;   /* zero occurs between A=A1 and A=A2 in domain dom */
  double A0;      /* q=0 at A=A0 in domain dom */
  double q1, q2;
  double X0, R0;  /* value of X,R at A=A0 */
  double B,phi, x,y,z;
  double ArgCp_1phi, AbsCp_1phi, ReCp_1phi, ImCp_1phi;
  double ArgCp_Bphi, ReCp_Bphi, ImCp_Bphi;
  double sigp_1phi, sigp_0phi, sigp_Bphi;
  int itmax = Geti("Coordinates_newtMAXITS");
  double tol = Getd("Coordinates_newtTOLF");
  double vec[2];
  int check, stat;
  int use_last_result_as_new_guess;

  /* look at B=1 (j=n2-1) and B=0 (j=0)  
     NOTE: we assue that n1,n2,n3 are the same in both domains */
  for(k=0, j=n2-1; j>=0; j-=n2-1)
  {
    /* find indices where q_in and q_out switch sign */
    for(i=n1-1; i>=0; i--) if(q_in[Index(i,j,k)]<=0.0) break;
    inz_in=i;
    for(i=0; i<n1; i++)   if(q_out[Index(i,j,k)]<=0.0) break;
    inz_out=i;

    /* if q<=0 already at i=n1-1 in innerdom, we quit! */
    if(inz_in==n1-1)
    {
      printf("reset_Coordinates_AnsorgNS_sigma_pm: innerdom=%d  B=%g phi=%g  "
             "inz_in=%d inz_out=%d\n", innerdom, B,phi, inz_in,inz_out);
      printf("q_in[Index(n1-1,j,k)]=%g\n", q_in[Index(n1-1,j,k)]);
      errorexit("reset_Coordinates_AnsorgNS_sigma_pm: q<=0 at i=n1-1 in "
                "inner domain???");
    }

    /* if inz_in=>0, q has zero in inner domain */
    /* if inz_out<n1, q is negative in outer domain */
    if(inz_in>=0)                   { i1=inz_in;   i2=inz_in+1; dom=innerdom;}
    else if(inz_out==0)             { i1=0;         i2=1;       dom=innerdom;}
    else if(inz_out<n1 && inz_out>0){ i1=inz_out-1; i2=inz_out; dom=outerdom;}
    else
    {
      printf("reset_Coordinates_AnsorgNS_sigma_pm: innerdom=%d  B=%g phi=%g  "
             "inz_in=%d inz_out=%d\n", innerdom, B,phi, inz_in,inz_out);
      printf("q_in[Index(n1-1,j,k)]=%g\n", q_in[Index(n1-1,j,k)]);
      errorexit("reset_Coordinates_AnsorgNS_sigma_pm: q>0 everywhere???");
    }
    A1 = grid->box[dom]->v[iX][Index(i1,j,k)];
    A2 = grid->box[dom]->v[iX][Index(i2,j,k)];
    q1 = grid->box[dom]->v[iq][Index(i1,j,k)];
    q2 = grid->box[dom]->v[iq][Index(i2,j,k)];
    B   = grid->box[dom]->v[iY][Index(i1,j,k)];
    phi = grid->box[dom]->v[iZ][Index(i1,j,k)];
//printf("inz_out=%d inz_in=%d\n", inz_out,inz_in);
//printf("reset_Coordinates_AnsorgNS_sigma_pm: "
//"find zero in q in dom=%d between A1=%g A2=%g\n", dom, A1,A2);

    /* find zero in q between A1 and A2 */
    if( fabs(q1) < tol || (inz_in<0 && dom==innerdom) ) A0=A1;
    else if(q1*q2>=0) A0=0.0;
    else /* use root finder */
    {
      A0 = A1 - q1*(A2-A1)/(q2-q1); /* initial guess */
      /* use newton_linesrch_its to find A0 */
      BNSdata_q_VectorFunc_box = grid->box[dom];
      BNSdata_q_VectorFunc_B   = B;
      BNSdata_q_VectorFunc_phi = phi;
      vec[1] = A0;
      stat=newton_linesrch_its(vec, 1, &check, BNSdata_q_VectorFunc, itmax, tol);
      if(check || stat<0)
        printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d stat=%d\n",
               check, stat);
      A0 = vec[1];
//printf("stat=%d: zero in q at A=A0=%g\n", stat, A0); Yo(1);
      /* go back to initial guess if we leave A1<A0<A2 in outerdom */
      if( (A0 < A1 || A0 > A2) && dom==outerdom )
        A0 = A1 - q1*(A2-A1)/(q2-q1);
      if(A0 < A1 - 0.05*(A2-A1) || A0 > A2 + 0.05*(A2-A1))
      {
        printf("reset_Coordinates_AnsorgNS_sigma_pm: B=%g phi=%g  "
               "inz_in=%d inz_out=%d\n", B,phi, inz_in,inz_out);
        printf("looked for zero in q in dom=%d between A1=%g A2=%g\n",
               dom, A1,A2);
        printf("q1=q(A1)=%g  q2=q(A2)=%g\n", q1,q2);
        printf("stat=%d: zero in q at A=A0=%g\n", stat, A0);        
        errorexit("reset_Coordinates_AnsorgNS_sigma_pm: "
                  "newton_linesrch_its failed!");
      }
      if(A0<0.0) A0=0.0; /* make sure we stay in box */
      if(A0>1.0) A0=1.0;
    }
//printf("zero in q at A=A0=%g\n", A0);

    /* compute values of X0,R0 at A=A0 */
    xyz_of_AnsorgNS(grid->box[dom], -1, dom, A0,B,phi, &x,&y,&z, &X0,&R0);

    /* get Cp and sigp at B=1  */
    if(j==n2-1) /* B=1 case */
    {
      ArgCp_1phi = acos(X0);
      /* 2 ArgCp = ArcTan[Sin[Pi B/2]/Sinh[sigma/2]] */
      /* Tan[2 ArgCp] = Sin[Pi B/2]/Sinh[sigma/2] */
      sigp_1phi = 2.0 * asinh( (1.0/tan(2.0*ArgCp_1phi)) );
      AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
      ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
      ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
      sigp_Bphi = sigp_1phi;
    }
    if(j==0) /* B=0 case */
    {
      ReCp_Bphi = X0;
      ImCp_Bphi = R0;
      ArgCp_Bphi = Arg(ReCp_Bphi, ImCp_Bphi);
      /* Cp_Bphi^2 = tanh(0.25*sigp_Bphi) */
      sigp_0phi = 4.0 * atanh(ReCp_Bphi*ReCp_Bphi-ImCp_Bphi*ImCp_Bphi);
      sigp_Bphi = sigp_0phi;
    }

    /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains 
       at B=1 and B=0 */
    for(kk=0; kk<n3; kk++)
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,kk)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,kk)] = sigp_Bphi;
      }
  } /* end for j */
//printf("reset_Coordinates_AnsorgNS_sigma_pm: new "
//       "sigp_0phi=%g sigp_1phi=%g\n", sigp_0phi, sigp_1phi);

#define SERIAL_reset_Coordinates_AnsorgNS_sigma_pm
  /* here we could use SGRID_LEVEL6_Pragma(omp parallel)  
     BUT: if use_last_result_as_new_guess=1:
     the loop below is a serial loop where we use vec[1] from j+1 as 
     initial guess for step j. So we have to set 
     use_last_result_as_new_guess=0 to use SGRID_LEVEL6_Pragma(omp for) .
     ALSO: xyz_of_AnsorgNS and dABphi_dxyz_AnsorgNS get slower if we have
           more than one thread, because static vars cannot be used then. */
#ifndef SERIAL_reset_Coordinates_AnsorgNS_sigma_pm
  SGRID_LEVEL6_Pragma(omp parallel)
  {
    tGrid *grid_p = make_empty_grid(grid->nvariables, 0);
    copy_grid(grid, grid_p, 0);
    /* how we pick initial guess in the loop below */
    use_last_result_as_new_guess=0;
#else
  {
    tGrid *grid_p = grid;
    /* how we pick initial guess in the loop below */
    use_last_result_as_new_guess=1;
#endif

    /* loop over the remaining j,k i.e. B,phi. 
       NOTE: we assume that n1,n2,n3 are the same in both domains */
    /* (j,k)-loop should be:
    for(j=n2-2; j>0; j--) // we could include j=0 (B=0) here again, so that most sigp_Bphi are found with the same method
    for(k=0; k<n3; k++)  */
#ifndef SERIAL_reset_Coordinates_AnsorgNS_sigma_pm
    SGRID_LEVEL6_Pragma(omp for)
#endif
    for(JK=n3*(n2-2)-1; JK>=0; JK--) /* JK = n3*(J-1) + K , J=1,...,n2-2, K=n3-1-k, K=0,...,n3-1 */
    {
      /* compute j,k from JK */
      int j = JK/n3 + 1;
      int k = (n3-1) - (JK%n3);
      double vec[2];
      double B,phi, sigp_Bphi, sigp_old;
      double w0, w1, wold;
      int check, stat;
      t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];

      /* set sigp_Bphi = sigp_1phi when we enter loop at j=n2-2, k=0 */
      if(JK==n3*(n2-2)-1) sigp_Bphi = sigp_1phi;

      /* find sigp_Bphi at B,phi such that q(sigp_Bphi; A=0, B, phi)=0 */
      B   = grid_p->box[dom]->v[iY][Index(0,j,k)];
      phi = grid_p->box[dom]->v[iZ][Index(0,j,k)];
      /* use newton_linesrch_its to find sigp_Bphi */
      pars->sigp_1phi = sigp_1phi;
      pars->B = B;
      pars->phi = phi;
      pars->grid = grid_p;
      pars->innerdom = innerdom;
      pars->outerdom = outerdom;
      if(use_last_result_as_new_guess)
      {
        /* guess for vec[1] is set to sigp_Bphi, the result of previous iteration */
        vec[1] = sigp_Bphi;
      }
      else
      {
        /* old sigp at B,phi and weights */
        sigp_old = grid_p->box[innerdom]->v[isigma][Index(0,j,k)];
        w0 = 1.0-Attenuation01(B*2.0, 2.5, 0.5);
        w1 = Attenuation01(B*2.0-1.0, 2.5, 0.5);
        wold = 1.0 - w0 - w1;
        /* guess for vec[1] is weighted average with weights w0,w1,wold */
        vec[1] = w0*sigp_0phi + w1*sigp_1phi + wold*sigp_old;
      }
//printf("itmax=%d tol=%g vec[1]=%g B=%g phi=%g\n",itmax,tol,vec[1], B,phi);
      stat=newton_linesrch_itsP(vec, 1, &check, q_of_sigp_forgiven_BphiP,
                                (void *) pars, itmax, tol);
      /* If q is nowhere negative newton_linesrch_its may not work. In this
         case we should probably search for the zero in (q - 1e-8). */
//printf("stat=%d\n",stat);
      if(check)
        printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d\n", check);  
      sigp_Bphi = vec[1];

      /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains */
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
      }
//printf("B=%g phi=%g  ", B, phi);
//printf("sigp_Bphi=%g sigp_0phi=%g sigp_1phi=%g\n", sigp_Bphi, sigp_0phi, sigp_1phi);
    } /* end for JK */
#ifdef SERIAL_reset_Coordinates_AnsorgNS_sigma_pm
  }
#else
    /* free local copies */
    free_grid(grid_p);
  }
#endif
  /* make sure that sigma has only one value at B=0 and also at B=1 */
  for(j=0; j<n2; j+=n2-1)
    for(k=1; k<n3; k++)
      for(i=0; i<n1; i++)
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] =
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] =
                        gridnew->box[innerdom]->v[isigma][Index(0,j,0)];

  /* compute derivs of sigma in both domains */
  spec_Deriv1(gridnew->box[innerdom], 2, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[innerdom], 3, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dphi]);
  spec_Deriv1(gridnew->box[outerdom], 2, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[outerdom], 3, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dphi]);
}


/******************************************************************/
/* useful functions                                               */
/******************************************************************/

/* compute ADM mass from point at A=1, B=0, phi=0 */
double ADMmass_fromPsi_inbox1_at_A1B0(tGrid *grid, int iADMmass)
{
  int iPsi   = Ind("BNSdata_Psi");
  int isigma = Ind("Coordinates_AnsorgNS_sigma_pm");
  double sigp_00 = grid->box[1]->v[isigma][0]; /* sigma_pm at A=B=phi=0 */
  double Cp_00   = sqrt( tanh(0.25*sigp_00) );
  double BNSdata_b = Getd("BNSdata_b");
  double *Psi     = grid->box[1]->v[iPsi];
  double *ADMmass = grid->box[1]->v[iADMmass]; 
  double M_ADM;
  int i;
  
  /* compute d_A d_A Psi */
  spec_Deriv2(grid->box[1], 1, Psi, ADMmass);

  /* multiply by 0.5*(BNSdata_b/(Cp_00*Cp_00)) to get ADM mass estimate */
  forallpoints(grid->box[1], i)
    ADMmass[i] *= 0.5*(BNSdata_b/(Cp_00*Cp_00));

  /* M_ADM = 0.5*(BNSdata_b/(Cp_00*Cp_00)) * d_A d_A Psi|A=1,B=0 */
  M_ADM = ADMmass[grid->box[1]->n1-1]; /* ADMmass at A=1, B=phi=0 */
  return M_ADM;
}


/* NOTE: do not use this anymore: */
/* compute volume integral of var with index vind in domain0+5 (if b=0)
   or domain3+4 (if b=3). Note, this version, adds a Psi^6 vol. el. factor. */
double InnerVolumeIntegral_withPsito6(tGrid *grid, int b, int vind)
{
  tGrid *grid2;
  double *var2;
  double *Psi2;
  double *Integ;
  double *Temp3;
  double *pX;
  double *pY;
  double *pZ;
  double *px;
  double *py;
  double *pz;
  double *dXdx;
  double *dXdy;
  double *dXdz;
  double *dYdx;
  double *dYdy;
  double *dYdz;
  double *dZdx;
  double *dZdy;
  double *dZdz;
  double *cv;
  double *cp;
  double box0_max1, box3_max1;  
  int box0_n1, box3_n1;
  char *box0_n1_sav = cmalloc( strlen(Gets("box0_n1"))+10 );
  char *box3_n1_sav = cmalloc( strlen(Gets("box3_n1"))+10 );
  int ib;
  int i,j,k;
  double Xmax;
  double VolInt;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  errorexit("use InnerVolumeIntegral instead of InnerVolumeIntegral_withPsito6");

  /* compute volume integral by making a new grid, where domain b
     covers the entire inside of the star */

  /* adjust box0 to cover the entire iside of star1 */
  box0_max1 = Getd("box0_max1");
  box0_n1   = Geti("box0_n1");
  strcpy(box0_n1_sav, Gets("box0_n1")); /* save box0_n1 */
  Sets("box0_max1", "1");
  Seti("box0_n1", box0_n1+Geti("box5_n1")/2);

  /* adjust box3 to cover the entire iside of star2 */
  box3_max1 = Getd("box3_max1");
  box3_n1   = Geti("box3_n1");
  strcpy(box3_n1_sav, Gets("box3_n1")); /* save box3_n1 */
  Sets("box3_max1", "1");
  Seti("box3_n1", box3_n1+Geti("box4_n1")/2);

  if(b==0) { ib=5; Xmax=box0_max1; }
  else     { ib=4; Xmax=box3_max1; }

  /* make grid with new adjusted boxes */
  grid2 = make_empty_grid(grid->nvariables, 0);
  set_BoxStructures_fromPars(grid2, 0);

  /* enable some vars needed on grid2 */
  enablevar(grid2, vind);
  enablevar(grid2, Ind("BNSdata_Psi"));
  enablevar(grid2, Ind("BNSdata_temp2"));
  enablevar(grid2, Ind("BNSdata_temp3"));
  enablevar(grid2, Ind("x"));
  enablevar(grid2, Ind("y"));
  enablevar(grid2, Ind("z"));
  enablevar(grid2, Ind("dXdx"));
  enablevar(grid2, Ind("dYdx"));
  enablevar(grid2, Ind("dZdx"));
  enablevar(grid2, Ind("Temp1"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_sigma_pm"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dB"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dphi"));  

  /* set Coordinates_AnsorgNS_sigma_pm, ... on grid2. Since grid2
     has different points only in the A-direc. we can copy it.    */  
  {
    int n1 = grid2->box[b]->n1;
    int n2 = grid2->box[b]->n2;
    int n3 = grid2->box[b]->n3;
    int N1 = grid->box[b]->n1;
    int N2 = grid->box[b]->n2;
    double *sigp2      = grid2->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
    double *dsigp_dB2  = grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")];
    double *dsigp_dphi2= grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")];
    double *sigp       = grid->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
    double *dsigp_dB   = grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")];
    double *dsigp_dphi = grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")];
    int i,j,k;

    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)
    for(i=0; i<n1; i++)
    {
      sigp2[(i)+n1*((j)+n2*(k))]       = sigp[0+N1*((j)+N2*(k))];
      dsigp_dB2[(i)+n1*((j)+n2*(k))]   = dsigp_dB[0+N1*((j)+N2*(k))];
      dsigp_dphi2[(i)+n1*((j)+n2*(k))] = dsigp_dphi[0+N1*((j)+N2*(k))];
    } 
  }

  /* initialize coords on grid2 */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid2);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* set var and Psi on grid2 */
  pX    = grid2->box[b]->v[Ind("X")];
  pY    = grid2->box[b]->v[Ind("Y")];
  pZ    = grid2->box[b]->v[Ind("Z")];
  px    = grid2->box[b]->v[Ind("x")];
  py    = grid2->box[b]->v[Ind("y")];
  pz    = grid2->box[b]->v[Ind("z")];
  dXdx  = grid2->box[b]->v[Ind("dXdx")];
  dXdy  = grid2->box[b]->v[Ind("dXdx")+1];
  dXdz  = grid2->box[b]->v[Ind("dXdx")+2];
  dYdx  = grid2->box[b]->v[Ind("dYdx")];
  dYdy  = grid2->box[b]->v[Ind("dYdx")+1];
  dYdz  = grid2->box[b]->v[Ind("dYdx")+2];
  dZdx  = grid2->box[b]->v[Ind("dZdx")];
  dZdy  = grid2->box[b]->v[Ind("dZdx")+1];
  dZdz  = grid2->box[b]->v[Ind("dZdx")+2];
  var2  = grid2->box[b]->v[vind];
  Integ = grid2->box[b]->v[Ind("BNSdata_temp2")];
  Temp3 = grid2->box[b]->v[Ind("BNSdata_temp3")];
  Psi2  = grid2->box[b]->v[Ind("BNSdata_Psi")];
  /* var   = grid->box[b]->v[vind]; */

  /* coeffs of var */
  spec_Coeffs(grid->box[b], grid->box[b]->v[vind], 
                            grid->box[b]->v[Ind("temp1")]);
  spec_Coeffs(grid->box[ib], grid->box[ib]->v[vind], 
                             grid->box[ib]->v[Ind("temp1")]);
  /* coeffs of Psi */
  spec_Coeffs(grid->box[b], grid->box[b]->v[Ind("BNSdata_Psi")], 
                            grid->box[b]->v[Ind("temp2")]);
  spec_Coeffs(grid->box[ib], grid->box[ib]->v[Ind("BNSdata_Psi")], 
                             grid->box[ib]->v[Ind("temp2")]);

  /* set var and Psi on grid2 by interpolation */
  forallpoints(grid2->box[b], i)
    if(pX[i]<Xmax)
    {
      cv = grid->box[b]->v[Ind("temp1")];
      cp = grid->box[b]->v[Ind("temp2")];
      var2[i] = spec_interpolate(grid->box[b], cv, pX[i], pY[i], pZ[i]);
      Psi2[i] = spec_interpolate(grid->box[b], cp, pX[i], pY[i], pZ[i]);
    }
    else
    {
      cv = grid->box[ib]->v[Ind("temp1")];
      cp = grid->box[ib]->v[Ind("temp2")];
      var2[i] = spec_interpolate(grid->box[ib], cv, px[i], py[i], pz[i]);
      Psi2[i] = spec_interpolate(grid->box[ib], cp, px[i], py[i], pz[i]);
    }

  /* set integrand in Integ */
  forallpoints(grid2->box[b], i)
  {
    double q   = var2[i];
    double Psi = Psi2[i];
    double Psi_to6 = Psi*Psi*Psi*Psi*Psi*Psi;
    double det = dXdx[i]*dYdy[i]*dZdz[i] + dXdy[i]*dYdz[i]*dZdx[i] +
                 dXdz[i]*dYdx[i]*dZdy[i] - dXdz[i]*dYdy[i]*dZdx[i] -
                 dXdy[i]*dYdx[i]*dZdz[i] - dXdx[i]*dYdz[i]*dZdy[i];
    double jac;

    if(det!=0.0) jac = 1.0/fabs(det);   else jac = 0.0;
    Integ[i] = q * Psi_to6 * jac;
  }

  /* integrate */
  VolInt = spec_3dIntegral(grid2->box[b], Integ, Temp3);

  /* remove grid2 */
  free_grid(grid2);

  /* reset box0/3 pars */
  Setd("box0_max1", box0_max1);
  Sets("box0_n1", box0_n1_sav);
  Setd("box3_max1", box3_max1);
  Sets("box3_n1", box3_n1_sav);

  free(box0_n1_sav);
  free(box3_n1_sav);

  return VolInt;
}

/* compute volume integral of var with index vind in domain0+5 (if b=0)
   or domain3+4 (if b=3). Here any Psi^6 needs to be already included
   in the var we integrate. */
double InnerVolumeIntegral(tGrid *grid, int b, int vind)
{
  tGrid *grid2;
  double *var2;
  double *Integ;
  double *Temp3;
  double *pX;
  double *pY;
  double *pZ;
  double *px;
  double *py;
  double *pz;
  double *dXdx;
  double *dXdy;
  double *dXdz;
  double *dYdx;
  double *dYdy;
  double *dYdz;
  double *dZdx;
  double *dZdy;
  double *dZdz;
  double *cv;
  double *cv_i;
  double box0_max1, box3_max1;  
  int box0_n1, box3_n1;
  char *box0_max1_sav = cmalloc( strlen(Gets("box0_max1"))+10 );
  char *box3_max1_sav = cmalloc( strlen(Gets("box3_max1"))+10 );
  char *box0_n1_sav = cmalloc( strlen(Gets("box0_n1"))+10 );
  char *box3_n1_sav = cmalloc( strlen(Gets("box3_n1"))+10 );
  int ib;
  int i,j,k;
  double Xmax;
  double VolInt;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  /* compute volume integral by making a new grid, where domain b
     covers the entire inside of the star */

  /* adjust box0 to cover the entire iside of star1 */
  box0_max1 = Getd("box0_max1");
  box0_n1   = Geti("box0_n1");
  strcpy(box0_max1_sav, Gets("box0_max1")); /* save box0_max1 */
  strcpy(box0_n1_sav, Gets("box0_n1"));     /* save box0_n1 */
  Sets("box0_max1", "1");
  Seti("box0_n1", box0_n1+Geti("box5_n1")/2);

  /* adjust box3 to cover the entire iside of star2 */
  box3_max1 = Getd("box3_max1");
  box3_n1   = Geti("box3_n1");
  strcpy(box3_max1_sav, Gets("box3_max1")); /* save box3_max1 */
  strcpy(box3_n1_sav, Gets("box3_n1"));     /* save box3_n1 */
  Sets("box3_max1", "1");
  Seti("box3_n1", box3_n1+Geti("box4_n1")/2);

  if(b==0) { ib=5; Xmax=box0_max1; }
  else     { ib=4; Xmax=box3_max1; }

  /* make grid with new adjusted boxes */
  grid2 = make_empty_grid(grid->nvariables, 0);
  set_BoxStructures_fromPars(grid2, 0);

  /* enable some vars needed on grid2 */
  enablevar(grid2, vind);
  enablevar(grid2, Ind("BNSdata_temp2"));
  enablevar(grid2, Ind("BNSdata_temp3"));
  enablevar(grid2, Ind("x"));
  enablevar(grid2, Ind("y"));
  enablevar(grid2, Ind("z"));
  enablevar(grid2, Ind("dXdx"));
  enablevar(grid2, Ind("dYdx"));
  enablevar(grid2, Ind("dZdx"));
  enablevar(grid2, Ind("Temp1"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_sigma_pm"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dB"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dphi"));  

  /* set Coordinates_AnsorgNS_sigma_pm, ... on grid2. Since grid2
     has different points only in the A-direc. we can copy it.    */  
  {
    int n1 = grid2->box[b]->n1;
    int n2 = grid2->box[b]->n2;
    int n3 = grid2->box[b]->n3;
    int N1 = grid->box[b]->n1;
    int N2 = grid->box[b]->n2;
    double *sigp2      = grid2->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
    double *dsigp_dB2  = grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")];
    double *dsigp_dphi2= grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")];
    double *sigp       = grid->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
    double *dsigp_dB   = grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")];
    double *dsigp_dphi = grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")];
    int i,j,k;

    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)
    for(i=0; i<n1; i++)
    {
      sigp2[(i)+n1*((j)+n2*(k))]       = sigp[0+N1*((j)+N2*(k))];
      dsigp_dB2[(i)+n1*((j)+n2*(k))]   = dsigp_dB[0+N1*((j)+N2*(k))];
      dsigp_dphi2[(i)+n1*((j)+n2*(k))] = dsigp_dphi[0+N1*((j)+N2*(k))];
    } 
  }

  /* initialize coords on grid2 */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid2);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* set var on grid2 */
  pX    = grid2->box[b]->v[Ind("X")];
  pY    = grid2->box[b]->v[Ind("Y")];
  pZ    = grid2->box[b]->v[Ind("Z")];
  px    = grid2->box[b]->v[Ind("x")];
  py    = grid2->box[b]->v[Ind("y")];
  pz    = grid2->box[b]->v[Ind("z")];
  dXdx  = grid2->box[b]->v[Ind("dXdx")];
  dXdy  = grid2->box[b]->v[Ind("dXdx")+1];
  dXdz  = grid2->box[b]->v[Ind("dXdx")+2];
  dYdx  = grid2->box[b]->v[Ind("dYdx")];
  dYdy  = grid2->box[b]->v[Ind("dYdx")+1];
  dYdz  = grid2->box[b]->v[Ind("dYdx")+2];
  dZdx  = grid2->box[b]->v[Ind("dZdx")];
  dZdy  = grid2->box[b]->v[Ind("dZdx")+1];
  dZdz  = grid2->box[b]->v[Ind("dZdx")+2];
  var2  = grid2->box[b]->v[vind];
  Integ = grid2->box[b]->v[Ind("BNSdata_temp2")];
  Temp3 = grid2->box[b]->v[Ind("BNSdata_temp3")];
  /* var   = grid->box[b]->v[vind]; */

  /* coeffs of var */
  spec_Coeffs(grid->box[b], grid->box[b]->v[vind], 
                            grid->box[b]->v[Ind("temp1")]);
  spec_Coeffs(grid->box[ib], grid->box[ib]->v[vind], 
                             grid->box[ib]->v[Ind("temp1")]);

  /* set var and Psi on grid2 by interpolation */
  cv = grid->box[b]->v[Ind("temp1")];
  cv_i = grid->box[ib]->v[Ind("temp1")];
  SGRID_LEVEL6_Pragma(omp parallel for)
  forallpoints(grid2->box[b], i)
  {
    if(pX[i]<Xmax)
      var2[i] = spec_interpolate(grid->box[b], cv, pX[i], pY[i], pZ[i]);
    else
      var2[i] = spec_interpolate(grid->box[ib], cv_i, px[i], py[i], pz[i]);
  }

  /* set integrand in Integ */
  forallpoints(grid2->box[b], i)
  {
    double q   = var2[i];
    double det = dXdx[i]*dYdy[i]*dZdz[i] + dXdy[i]*dYdz[i]*dZdx[i] +
                 dXdz[i]*dYdx[i]*dZdy[i] - dXdz[i]*dYdy[i]*dZdx[i] -
                 dXdy[i]*dYdx[i]*dZdz[i] - dXdx[i]*dYdz[i]*dZdy[i];
    double jac;

    if(det!=0.0) jac = 1.0/fabs(det);   else jac = 0.0;
    Integ[i] = q * jac;
  }

  /* integrate */
  VolInt = spec_3dIntegral(grid2->box[b], Integ, Temp3);

  /* remove grid2 */
  free_grid(grid2);

  /* reset box0/3 pars */
  Sets("box0_max1", box0_max1_sav);
  Sets("box0_n1", box0_n1_sav);
  Sets("box3_max1", box3_max1_sav);
  Sets("box3_n1", box3_n1_sav);

  free(box0_max1_sav);
  free(box3_max1_sav);
  free(box0_n1_sav);
  free(box3_n1_sav);

  return VolInt;
}

/* compute volume integral of var with index vind in domain1 (if b=1)
   or domain2 (if b=2), but it should also work in the other domains.
   Here any Psi^6 needs to be already included in the var we integrate. */
double VolumeIntegral_inBNSgridBox(tGrid *grid, int b, int vind)
{
  double *var;
  double *Integ;
  double *Temp3;
  double *dXdx;
  double *dXdy;
  double *dXdz;
  double *dYdx;
  double *dYdy;
  double *dYdz;
  double *dZdx;
  double *dZdy;
  double *dZdz;
  double VolInt;
  int i;

  var   = grid->box[b]->v[vind];
  dXdx  = grid->box[b]->v[Ind("dXdx")];
  dXdy  = grid->box[b]->v[Ind("dXdx")+1];
  dXdz  = grid->box[b]->v[Ind("dXdx")+2];
  dYdx  = grid->box[b]->v[Ind("dYdx")];
  dYdy  = grid->box[b]->v[Ind("dYdx")+1];
  dYdz  = grid->box[b]->v[Ind("dYdx")+2];
  dZdx  = grid->box[b]->v[Ind("dZdx")];
  dZdy  = grid->box[b]->v[Ind("dZdx")+1];
  dZdz  = grid->box[b]->v[Ind("dZdx")+2];
  Integ = grid->box[b]->v[Ind("BNSdata_temp2")];
  Temp3 = grid->box[b]->v[Ind("BNSdata_temp3")];

  /* set integrand in Integ */
  forallpoints(grid->box[b], i)
  {
    double det = dXdx[i]*dYdy[i]*dZdz[i] + dXdy[i]*dYdz[i]*dZdx[i] +
                 dXdz[i]*dYdx[i]*dZdy[i] - dXdz[i]*dYdy[i]*dZdx[i] -
                 dXdy[i]*dYdx[i]*dZdz[i] - dXdx[i]*dYdz[i]*dZdy[i];
    double jac;

    if(det!=0.0) jac = 1.0/fabs(det);
    /* if det=0 jac should really be infinite, but we hope that the
       integrand goes to zero quickly enough that jac=0 makes difference! */
    else jac = 0.0;

    Integ[i] = var[i] * jac;
  }

  /* integrate */
  VolInt = spec_3dIntegral(grid->box[b], Integ, Temp3);

  return VolInt;
}


/* make a grid with a different number of points in A-,B- and phi-dir */
tGrid *make_grid_with_sigma_pm(tGrid *grid, int nAB, int nphi, int nxyz)
{
  int pr=0;
  tGrid *grid2;
  char *box0_n1_sav = strdup(Gets("box0_n1"));
  char *box1_n1_sav = strdup(Gets("box1_n1"));
  char *box2_n1_sav = strdup(Gets("box2_n1"));
  char *box3_n1_sav = strdup(Gets("box3_n1"));
  char *box0_n2_sav = strdup(Gets("box0_n2"));
  char *box1_n2_sav = strdup(Gets("box1_n2"));
  char *box2_n2_sav = strdup(Gets("box2_n2"));
  char *box3_n2_sav = strdup(Gets("box3_n2"));
  char *box0_n3_sav = strdup(Gets("box0_n3"));
  char *box1_n3_sav = strdup(Gets("box1_n3"));
  char *box2_n3_sav = strdup(Gets("box2_n3"));
  char *box3_n3_sav = strdup(Gets("box3_n3"));
  int b;
  int i,j,k;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  /* adjust n1,n2,n3 in boxes */
  Seti("box0_n1", nAB);
  Seti("box1_n1", nAB);
  Seti("box2_n1", nAB);
  Seti("box3_n1", nAB);
  Seti("box0_n2", nAB);
  Seti("box1_n2", nAB);
  Seti("box2_n2", nAB);
  Seti("box3_n2", nAB);
  Seti("box0_n3", nphi);
  Seti("box1_n3", nphi);
  Seti("box2_n3", nphi);
  Seti("box3_n3", nphi);

  /* make grid with new adjusted boxes */
  grid2 = make_empty_grid(grid->nvariables, pr);
  set_BoxStructures_fromPars(grid2, pr);

  /* enable some vars needed on grid2 */
  enablevar(grid2, Ind("x"));
  enablevar(grid2, Ind("y"));
  enablevar(grid2, Ind("z"));
  enablevar(grid2, Ind("dXdx"));
  enablevar(grid2, Ind("dYdx"));
  enablevar(grid2, Ind("dZdx"));
  enablevar(grid2, Ind("Temp1"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_sigma_pm"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dB"));
  enablevar(grid2, Ind("Coordinates_AnsorgNS_dsigma_pm_dphi"));  

  /* coeffs of Coordinates_AnsorgNS_sigma_pm and its derivs */
  for(b=0; b<=3; b++)
  {
    spec_Coeffs(grid->box[b],
                grid->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")], 
                grid->box[b]->v[Ind("temp1")]);
    spec_Coeffs(grid->box[b],
                grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")], 
                grid->box[b]->v[Ind("temp2")]);
    spec_Coeffs(grid->box[b],
                grid->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")], 
                grid->box[b]->v[Ind("temp3")]);
  }

  /* set var and Psi on grid2 by interpolation */
  for(b=0; b<=3; b++)
  {
    double *c1 = grid->box[b]->v[Ind("temp1")];
    double *c2 = grid->box[b]->v[Ind("temp2")];
    double *c3 = grid->box[b]->v[Ind("temp3")];
    double *sig   =grid2->box[b]->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
    double *sigB  =grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dB")];
    double *sigphi=grid2->box[b]->v[Ind("Coordinates_AnsorgNS_dsigma_pm_dphi")];
    double *X = grid2->box[b]->v[Ind("X")];
    double *Y = grid2->box[b]->v[Ind("Y")];
    double *Z = grid2->box[b]->v[Ind("Z")];
    int n1 = grid2->box[b]->n1;
    int n2 = grid2->box[b]->n2;
    int n3 = grid2->box[b]->n3;
    int ind;

    /* interpolate for i=0 <===> A=0 */
    i=0; /* A=0 */
    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)
    {
      ind = Index(0,j,k);
      sig[ind]  = spec_interpolate(grid->box[b], c1, X[ind], Y[ind], Z[ind]);
      sigB[ind] = spec_interpolate(grid->box[b], c2, X[ind], Y[ind], Z[ind]);
      sigphi[ind]=spec_interpolate(grid->box[b], c3, X[ind], Y[ind], Z[ind]);
      if(pr) 
       printf("(X,Y,Z)=(%g,%g,%g) (i,j,k)=(%d,%d,%d) sig=%g sigB=%g sigphi=%g\n",
              X[ind],Y[ind],Z[ind], i,j,k, sig[ind],sigB[ind],sigphi[ind]);
    }
    /* copy for all other i */
    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)
    for(i=1; i<n1; i++)
    {
      sig[Index(i,j,k)]    = sig[Index(0,j,k)];
      sigB[Index(i,j,k)]   = sigB[Index(0,j,k)];
      sigphi[Index(i,j,k)] = sigphi[Index(0,j,k)];
    }
  }

  /* initialize coords on grid2 */
  if(Coordinates_verbose && pr==0) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid2);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* reset box pars */
  Sets("box0_n1", box0_n1_sav);
  Sets("box1_n1", box1_n1_sav);
  Sets("box2_n1", box2_n1_sav);
  Sets("box3_n1", box3_n1_sav);
  Sets("box0_n2", box0_n2_sav);
  Sets("box1_n2", box1_n2_sav);
  Sets("box2_n2", box2_n2_sav);
  Sets("box3_n2", box3_n2_sav);
  Sets("box0_n3", box0_n3_sav);
  Sets("box1_n3", box1_n3_sav);
  Sets("box2_n3", box2_n3_sav);
  Sets("box3_n3", box3_n3_sav);
  
  /* free saved strings */
  free(box0_n1_sav);
  free(box1_n1_sav);
  free(box2_n1_sav);
  free(box3_n1_sav);
  free(box0_n2_sav);
  free(box1_n2_sav);
  free(box2_n2_sav);
  free(box3_n2_sav);
  free(box0_n3_sav);
  free(box1_n3_sav);
  free(box2_n3_sav);
  free(box3_n3_sav);

  return grid2;
}


/* figure out max A inside star in innerdom and adjust box4/5 accordingly 
   if innerdom<0 we adjust both box 4 and 5. */
void adjust_box4_5_pars_pm(tGrid *grid, int innerdom)
{
  double box0_max1 = grid->box[0]->bbox[1];
  double box3_max1 = grid->box[3]->bbox[1];
  double scal = 1.05; /* make box4/5 5% larger than needed in x-dir */
  double scal2= 1.05; /* make box4/5 5% larger than needed in y/z-dir */
  double Bstep= 0.001;
  double xp, xm, xmax, xmin;
  double ymin,ymax, zmin,zmax, res, B;
  double b = GetCachedNumValByParIndex(Coordinates_AnsorgNS_b_ParIndex);    
  int adjbox4=0;
  int adjbox5=0;

  /* which box needs adjustment */
  if(innerdom==0) adjbox5=1;
  if(innerdom==3) adjbox4=1;
  if(innerdom<0)  { adjbox5=1; adjbox4=1; }

  /* adjust box5 */
  if(adjbox5)
  {
    ymin=ymax=zmin=zmax = 0.0;
    xmin=xmax = b;
    for(B=0.0; B<1.0; B+=Bstep)
    {
      res = x_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,0.0);
      if(res<xmin) xmin=res;
      if(res>xmax) xmax=res;

      res = x_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,PI*0.5);
      if(res<xmin) xmin=res;
      if(res>xmax) xmax=res;

      res = y_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,PI);
      if(res<ymin) ymin=res;

      res = y_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,0.0);
      if(res>ymax) ymax=res;

      res = z_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,PI*1.5);
      if(res<zmin) zmin=res;

      res = z_of_AnsorgNS0(grid->box[0], -1, box0_max1,B,PI*0.5);
      if(res>zmax) zmax=res;
    }
    xm = scal * (xmin-b);
    xp = scal * (xmax-b);
    ymin = scal2 * ymin;
    ymax = scal2 * ymax;
    zmin = scal2 * zmin;
    zmax = scal2 * zmax;

    Setd("box5_min1", b + xm);
    Setd("box5_max1", b + xp);
    Setd("box5_min2", ymin);
    Setd("box5_max2", ymax);
    Setd("box5_min3", zmin);
    Setd("box5_max3", zmax);
  }

  /* adjust box4 */
  if(adjbox4)
  {
    ymin=ymax=zmin=zmax = 0.0;
    xmin=xmax = -b;
    for(B=0.0; B<1.0; B+=Bstep)
    {
      res = x_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,0.0);
      if(res<xmin) xmin=res;
      if(res>xmax) xmax=res;

      res = x_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,PI*0.5);
      if(res<xmin) xmin=res;
      if(res>xmax) xmax=res;

      res = y_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,PI);
      if(res<ymin) ymin=res;

      res = y_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,0.0);
      if(res>ymax) ymax=res;

      res = z_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,PI*1.5);
      if(res<zmin) zmin=res;

      res = z_of_AnsorgNS3(grid->box[3], -1, box3_max1,B,PI*0.5);
      if(res>zmax) zmax=res;
    }
    xm = scal * (xmin+b);
    xp = scal * (xmax+b);
    ymin = scal2 * ymin;
    ymax = scal2 * ymax;
    zmin = scal2 * zmin;
    zmax = scal2 * zmax;

    Setd("box4_min1", -b + xm);
    Setd("box4_max1", -b + xp);
    Setd("box4_min2", ymin);
    Setd("box4_max2", ymax);
    Setd("box4_min3", zmin);
    Setd("box4_max3", zmax);
  }
}
/* figure out max A inside stars and adjust boxes4/5 accordingly */
void adjust_box4_5_pars(tGrid *grid)
{
  adjust_box4_5_pars_pm(grid, -1); /* adjust both box4 and 5 */
}


/* find the box and the coords of the point at the cartesian x,y,z */
/* initially b, *X,*Y,*Z contain the boxindex and the coords of the
   point on the other grid */
/* NOTE: Interp_Var_From_Grid1_To_Grid2_pm has certain smoothing 
         because of BNSgrid_Get_BoxAndCoords_of_xyz (see below). */
int BNSgrid_Get_BoxAndCoords_of_xyz(tGrid *grid1,
                                    double *X1, double *Y1, double *Z1,
                                    int b, double x, double y, double z)
{
  static int BNSdata_box0_Amax_ParIndex=-1;
  static int BNSdata_box3_Amax_ParIndex=-1;
  double Aguess_forb_5;
  double Aguess_forb_4;
  int b1;
  double X = *X1;
  double Y = *Y1;
  double Z = *Z1;
  int blist[6];

//if(dequal(Z, 0.0))
//printf("b =%d  X=%.4g Y=%.4g Z=%.4g  x=%g y=%g z=%g\n",b, X,Y,Z, x,y,z); //Yo(1);

  /* depending on b decide how to obtain X,Y,Z */
  if( (b==0 || b==3 || b==5 || b==4) )
  {
    blist[0]=4;  blist[1]=5;
    b1 = b_XYZ_of_xyz_inboxlist(grid1, blist,2, &X,&Y,&Z, x,y,z);
    if(b1<0) /* not in box4/5 */
    {
      /* set good guesses for Z=phi and X=A, Y=B */
      if(b>=4)
      { 
        Z = Arg(Y,Z); if(Z<0.0) Z+=2.0*PI;
        if(BNSdata_box0_Amax_ParIndex<0 || BNSdata_box3_Amax_ParIndex<0) /* init Aguess once */
        { 
          BNSdata_box0_Amax_ParIndex = GetParIndex("BNSdata_box0_Amax");
          BNSdata_box3_Amax_ParIndex = GetParIndex("BNSdata_box3_Amax");
          //Aguess_forb_5 = Getd("BNSdata_box0_Amax");
          //Aguess_forb_4 = Getd("BNSdata_box3_Amax");
          Aguess_forb_5 = GetCachedNumValByParIndex(BNSdata_box0_Amax_ParIndex);
          Aguess_forb_4 = GetCachedNumValByParIndex(BNSdata_box3_Amax_ParIndex);
          printf("BNSgrid_Get_BoxAndCoords_of_xyz: "
                 "Aguess_forb_5=%g Aguess_forb_4=%g\n"
                 " BNSdata_box0_Amax_ParIndex=%d"
                 " BNSdata_box3_Amax_ParIndex=%d\n",
                  Aguess_forb_5, Aguess_forb_4,
                  BNSdata_box0_Amax_ParIndex, BNSdata_box3_Amax_ParIndex);
        }
        Aguess_forb_5 = GetCachedNumValByParIndex(BNSdata_box0_Amax_ParIndex);
        Aguess_forb_4 = GetCachedNumValByParIndex(BNSdata_box3_Amax_ParIndex);
        if(b==5) X = Aguess_forb_5; /* =0.85, bad guess ??? */
        else     X = Aguess_forb_4; /* =0.85, bad guess ??? */
        Y = 0.5;  /* bad guess ??? */
      }
      if( b<4 && (dequal(Y, 0.0) || dequal(Y, 1.0)) )
      { 
        if(b<2) { blist[0]=0;  blist[1]=1; }
        else    { blist[0]=3;  blist[1]=2; }
        b1 = b_X_of_x_forgiven_YZ_inboxlist(grid1, blist,2, &X, x, Y,Z);
      }
      else if(b==0 || b==5)
      {
        blist[0]=0;  blist[1]=1;
        b1 = b_XYZ_of_xyz_inboxlist(grid1, blist,2, &X,&Y,&Z,
                                    x,y,z);
      }
      else if(b==3 || b==4)
      {
        blist[0]=3;  blist[1]=2;
        b1 = b_XYZ_of_xyz_inboxlist(grid1, blist,2, &X,&Y,&Z,
                                    x,y,z);
      }
    } /* end: not in box4/5 */
  }
  else if(b==1)
  {
    if( (dequal(Y, 0.0) || dequal(Y, 1.0)) )
    {
      blist[0]=0;  blist[1]=1;
      b1 = b_X_of_x_forgiven_YZ_inboxlist(grid1, blist,2, &X, x, Y,Z);
    }
    else
    {
      blist[0]=1;  blist[1]=0;
      b1 = b_XYZ_of_xyz_inboxlist(grid1, blist,2, &X,&Y,&Z,
                                  x,y,z);
    }
  }
  else /* b==2 */
  {
    if( (dequal(Y, 0.0) || dequal(Y, 1.0)) )
    {
      blist[0]=3;  blist[1]=2;
      b1 = b_X_of_x_forgiven_YZ_inboxlist(grid1, blist,2, &X, x, Y,Z);
    }
    else
    {
      blist[0]=2;  blist[1]=3;
      b1 = b_XYZ_of_xyz_inboxlist(grid1, blist,2, &X,&Y,&Z,
                                  x,y,z);
    }
  }

  /* failure? */
  if(b1<0)
  {
    printf("b1=%d  initial guess b=%d *X1=%.4g *Y1=%.4g *Z1=%.4g\n",
           b1, b,*X1,*Y1,*Z1);
    printf("b1=%d  X=%.4g Y=%.4g Z=%.4g  x=%g y=%g z=%g\n",
           b1, X,Y,Z, x,y,z);
  }
  *X1 = X;
  *Y1 = Y;
  *Z1 = Z;
  return b1;
}

/* Interpolate Var with index vind from grid1 to grid2 */
void Interpolate_Var_From_Grid1_To_Grid2(tGrid *grid1, tGrid *grid2, int vind)
{
  Interp_Var_From_Grid1_To_Grid2_pm(grid1, grid2, vind, 0);
  Interp_Var_From_Grid1_To_Grid2_pm(grid1, grid2, vind, 3);
}
/* Interpolate Var with index vind from grid1 to grid2 
   for domains centered on innerdom */
/* NOTE: Interp_Var_From_Grid1_To_Grid2_pm has certain smoothing properties
         since BNSgrid_Get_BoxAndCoords_of_xyz preferably returns coords
         in boxes 5/4 (when possible) where the solution is smoother.
         Thus we interpolate from 5/4 into 0/3, 
         which smooths in box 0/3. */
void Interp_Var_From_Grid1_To_Grid2_pm(tGrid *grid1, tGrid *grid2, int vind,
                                       int innerdom)
{
  int cind = Ind("temp1");
  int Xind = Ind("X");
  int Yind = Ind("Y");
  int Zind = Ind("Z");
  int xind = Ind("x");
  int yind = Ind("y");
  int zind = Ind("z");
  int b,i;

  /* save coeffs of vind on grid1 in cind = Ind("Temp1") */
  forallboxes(grid1, b)
  {
    tBox *box = grid1->box[b];
    /* do nothing if we are on wrong side of grid */
    if( (innerdom==0 && (b>=2 && b<=4)) || (innerdom==3 && (b<=1 || b>=5)) )
      continue;

    spec_Coeffs(box, box->v[vind], box->v[cind]);
  }

  /* loop over grid2 */
  forallboxes(grid2,b)
  {
    tBox *box = grid2->box[b];
    double *pX = box->v[Xind];
    double *pY = box->v[Yind];
    double *pZ = box->v[Zind];
    double *px = box->v[xind];
    double *py = box->v[yind];
    double *pz = box->v[zind];
    double *pv = box->v[vind];

    /* do nothing if we are on wrong side of grid */
    if( (innerdom==0 && (b>=2 && b<=4)) || (innerdom==3 && (b<=1 || b>=5)) )
      continue;

    /* here we can use SGRID_LEVEL6_Pragma(omp parallel) */
#undef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
#ifndef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
    SGRID_LEVEL6_Pragma(omp parallel)
    {
      tGrid *grid1_p = make_empty_grid(grid1->nvariables, 0);
      copy_grid(grid1, grid1_p, 0);
#else
    {
      tGrid *grid1_p = grid1;
#endif
      /* start loop */
#ifndef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
      SGRID_LEVEL6_Pragma(omp for) 
#endif
      forallpoints(box,i)
      {
        double X = pX[i];
        double Y = pY[i];
        double Z = pZ[i];
        int b1;

        /* get b1, X,Y,Z on grid1_p */
        b1 = BNSgrid_Get_BoxAndCoords_of_xyz(grid1_p, &X,&Y,&Z, 
                                             b,px[i],py[i],pz[i]);
        if(b1<0)
        {
          double x,y,z;
          double *sigpm = box->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
          printf("b1=%d grid2: b=%d i=%d X=%g Y=%g Z=%g "
                 "sigpm[i]=%g sigpm[box->n1*(box->n1-1)]=%g\n",
                 b1, b,i,X,Y,Z, sigpm[i], sigpm[box->n1*(box->n1-1)]);
          printf("b1=%d grid2: b=%d x=%g y=%g z=%g\n", b1, b,px[i],py[i],pz[i]);
          box->x_of_X[3]((void *) box, -1,X*0.5,Y*0.5,Z);
          x=box->x_of_X[1]((void *) box, i,X,Y,Z);
          y=box->x_of_X[2]((void *) box, i,X,Y,Z);
          z=box->x_of_X[3]((void *) box, i,X,Y,Z);
          printf("b1=%d grid2: box->x_of_X => x=%g y=%g z=%g\n", b1, x,y,z);
          errorexit("Interp_Var_From_Grid1_To_Grid2_pm: "
                    "could not find X,Y,Z on grid1_p.");
        }

        /* get var at point X,Y,Z by interpolation */
        pv[i] = spec_interpolate(grid1_p->box[b1], grid1_p->box[b1]->v[cind], X,Y,Z);
//if(!finite(pv[i]))
//{
//double x,y,z;
// //double *sigpm = box->v[Ind("Coordinates_AnsorgNS_sigma_pm")];
// //printf("b1=%d grid2: b=%d i=%d X=%g Y=%g Z=%g "
// //       "sigpm[i]=%g sigpm[box->n1*(box->n1-1)]=%g\n",
// //       b1, b,i,X,Y,Z, sigpm[i], sigpm[box->n1*(box->n1-1)]);
//printf("VarName(vind)=%s VarName(cind)=%s\n", VarName(vind), VarName(cind));
//printf("b1=%d grid2: b=%d i=%d X=%g Y=%g Z=%g\n", b1, b,i,X,Y,Z);
//printf("b1=%d grid2: b=%d      x=%g y=%g z=%g\n", b1, b,px[i],py[i],pz[i]);
//CheckIfFinite(grid2, VarName(vind));
//CheckIfFinite(grid2, VarName(cind));
//exit(77);
//}
      } /* end forallpoints loop */
#ifdef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
    }
#else
      /* free local grid copy */
      free_grid(grid1_p);
    }
#endif
  } /* end forallboxes(grid2,b) */
}
/* wrapper for Interpolate_Var_From_Grid1_To_Grid2 with extra dummy argument */
void Interpolate_Var_From_Grid1_To_Grid2_wrapper(tGrid *grid1, tGrid *grid2,
                                                 int vind, int dummy)
{
  Interpolate_Var_From_Grid1_To_Grid2(grid1, grid2, vind);
}


/* copy variable at i=0 from Box1 with index b1 to Box2 with index b2 */
/* This can be used to make a var continues across inner and outer domains. */
void copy_Var_at_i0_from_Box1_Box2(tGrid *grid, int vind, int b1, int b2)
{
  int i,j,k, ijk1, ijk2; 
  int n1_1 = grid->box[b1]->n1;
  int n2_1 = grid->box[b1]->n2;
  int n3_1 = grid->box[b1]->n3;
  int n1_2 = grid->box[b2]->n1;
  int n2_2 = grid->box[b2]->n2;
  int n3_2 = grid->box[b2]->n3;
  double *v1 = grid->box[b1]->v[vind];
  double *v2 = grid->box[b2]->v[vind];
  if(n2_2!=n2_1 || n3_2!=n3_1)
    errorexit("(n2,n3) has to be equal in both boxes");
  forplane1(i,j,k, n1_2,n2_2,n3_2, 0)
  {
    ijk1 = Ind_n1n2(0,j,k,n1_1,n2_1);
    ijk2 = Ind_n1n2(0,j,k,n1_2,n2_2);
    v2[ijk2] = v1[ijk1];
  }
}


/* enforce uniqueness on axis by copying from phi=0 to all other phi values
   on the axis */
void BNS_enforce_uniqueness_on_axis(tVarList *vlu)
{
  tGrid *grid = vlu->grid;
  int b, i;
  int vind;

  for(vind=0; vind<vlu->n; vind++)
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;
    double *var  = box->v[vlu->index[vind]];
 
    if(b<=3) 
    {
      char str[1000];
      snprintf(str, 999, "box%d_basis2", b);
      if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
      {
        /* loop over rho=0 boundary */
        for(k=1; k<n3; k++)       /* <-- all phi>0 */
        for(j=0; j<n2; j=j+n2-1)  /* <-- B=0 and B=1 */
        for(i=0; i<n1; i++)       /* <-- all A */
          var[Index(i,j,k)] = var[Index(i,j,0)];
      }
    }
  } /* end forallboxes */
}


/* compute weighted average of the new q2 on grid2 and the old q1 on grid1 
   at a point X2,Y2,Z2 in grid2 coords */
double BNS_update_q_atXYZ(tGrid *grid2, 
                          int b2, double X2, double Y2, double Z2,
                          double w, tGrid *grid1)
{
  double x,y,z, Xp,Rp;
  int b1;
  double X1,Y1,Z1;
  double q2, q1;

  /* get q on grid 2 at X2,Y2,Z2 */
  q2 = BNS_compute_new_centered_q_atXYZ(grid2,b2, X2,Y2,Z2);

  /* get q on grid 1 at the same point */
  if(b2<4)
    xyz_of_AnsorgNS(grid2->box[b2], -1, b2, X2,Y2,Z2, &x,&y,&z, &Xp,&Rp);
  else
    { x=X2;  y=Y2;  z=Z2; }
  X1 = X2;  Y1 = Y2;  Z1 = Z2;
  b1 = BNSgrid_Get_BoxAndCoords_of_xyz(grid1, &X1,&Y1,&Z1, b2,x,y,z);
  q1 = BNS_compute_new_centered_q_atXYZ(grid1,b1, X1,Y1,Z1);

  /* return weighted average */
  return w*q2 + (1.0-w)*q1;
}

/* compute weighted average of the new q2 on grid2 and the old q1 on grid1 
   at a point X2,Y2,Z2 in grid2 coords */
void BNS_update_q(tGrid *grid2, double w, tGrid *grid1)
{
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int iq = Ind("BNSdata_q");
  int b2, i;

  forallboxes(grid2,b2)
  {
    tBox *box = grid2->box[b2];
    double *X2 = box->v[iX];
    double *Y2 = box->v[iY];
    double *Z2 = box->v[iZ];
    double *q = box->v[iq];

    forallpoints(box, i)
      q[i] = BNS_update_q_atXYZ(grid2,b2, X2[i],Y2[i],Z2[i], w, grid1);
  }
}


/************************************************************************/
/* utilities to load checkpoints with different resolution and pars */
/************************************************************************/
/* for newton_linesrch_itsP: compute error in rest mass */
void m0errOFsigmafac_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int dom;
  double fac = vec[1];
  double m0;
  t_grid_b_struct *pars;

  /* get pars */
  pars = (t_grid_b_struct *) p;
  grid = pars->grid;
  dom  = pars->b;

  /* scale sigma */
  BNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, fac, dom);
  BNSgrid_init_Coords(grid);
  /* set wB */
  if(dom==0) BNS_set_wB(grid, 1, Getd("BNSdata_xmax1"),0.0,0.0);
  else       BNS_set_wB(grid, 2, Getd("BNSdata_xmax2"),0.0,0.0); 
  /* get new mass */
  m0 =  GetInnerRestMass(grid, dom);
  printf("m0errOFsigmafac_VectorFuncP: fac=%g => m0%d=%g\n",
         fac, (dom>=2)+1, m0);
  fflush(stdout);
  /* set sigma back to what it was */
  BNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, 1.0/fac, dom);
  BNSgrid_init_Coords(grid);
  /* set wB */
  if(dom==0) BNS_set_wB(grid, 1, Getd("BNSdata_xmax1"),0.0,0.0);
  else       BNS_set_wB(grid, 2, Getd("BNSdata_xmax2"),0.0,0.0); 
  /* set fvec */
  if(dom==0)
    fvec[1] = Getd("BNSdata_m01") - m0;
  else
    fvec[1] = Getd("BNSdata_m02") - m0;
  //printf("  => m0_err=%g\n", fvec[1]);
  //fflush(stdout);
}

/* load data from a previous result in a checkpoint file 
   (with e.g. lower res) */
void BNSgrid_load_initial_guess_from_checkpoint(tGrid *grid, char *filename)
{
  tVarList *varlist;
  char *parlist;
  int j, bi;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");
  char *BNSdata_b_sav = strdup(Gets("BNSdata_b"));
  char *BNSdata_x_CM_sav = strdup(Gets("BNSdata_x_CM"));
  char *BNSdata_xmax1_sav = strdup(Gets("BNSdata_xmax1"));
  char *BNSdata_xmax2_sav = strdup(Gets("BNSdata_xmax2"));

  /* make list of pars we want to read */  
  parlist = "BNSdata_b BNSdata_m01 BNSdata_m02 BNSdata_kappa "
            "BNSdata_Omega BNSdata_x_CM BNSdata_C1 BNSdata_C2 "
            "BNSdata_qmax1 BNSdata_qmax2 BNSdata_xmax1 BNSdata_xmax2 "
            "BNSdata_actual_xmax1 BNSdata_actual_xmax2";

  /* make list of vars we want to read */
  varlist = vlalloc(grid);
  vlpush(varlist, Ind("BNSdata_Psi"));
  vlpush(varlist, Ind("BNSdata_Bx"));
  vlpush(varlist, Ind("BNSdata_alphaP"));
  vlpush(varlist, Ind("BNSdata_Sigma"));
  vlpush(varlist, Ind("BNSdata_q"));
  vlpush(varlist, Ind("Coordinates_AnsorgNS_sigma_pm"));

//printVarList(varlist);

  /* read vars in varlist and pars in parlist*/
  checkpoint_interpolate_Vars_get_Pars(filename, varlist, parlist);

  /* once we have a new Coordinates_AnsorgNS_sigma_pm, we need re-init
     the Coords again. */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* set q to zero if q<0, and also in region 1 & 2 */
  forallboxes(grid, bi)
  {
    int i;
    double *BNSdata_q = grid->box[bi]->v[Ind("BNSdata_q")];
    forallpoints(grid->box[bi], i)
      if( BNSdata_q[i]<0.0 || bi==1 || bi==2 )  BNSdata_q[i] = 0.0;
  }

  if(strcmp(BNSdata_b_sav, Gets("BNSdata_b"))==0)
  {
    /* set values of A,B,phi in box4/5 */
    set_BNSdata_ABphi(grid);

    /* set wB */
    BNS_set_wB(grid, 1, Getd("BNSdata_actual_xmax1"),0.0,0.0);
    BNS_set_wB(grid, 2, Getd("BNSdata_actual_xmax2"),0.0,0.0); 
  }
  else /* set BNSdata_b=BNSdata_b_sav and adjust other pars */
  {
    double m01 = Getd("BNSdata_m01");
    double m02 = Getd("BNSdata_m02");
    double r, rc, dr, dOm;
    double facvec[2];
    int stat, check;
    t_grid_b_struct pars[1];

    /* rc = (distance from checkpoint) */
    rc = Getd("BNSdata_xmax1")-Getd("BNSdata_xmax2");
    Sets("BNSdata_b", BNSdata_b_sav); /* set b back to saved value */
    /* set some things back to saved values */
    Sets("BNSdata_x_CM", BNSdata_x_CM_sav);
    Sets("BNSdata_xmax1", BNSdata_xmax1_sav);
    Sets("BNSdata_xmax2", BNSdata_xmax2_sav);
    r = Getd("BNSdata_xmax1")-Getd("BNSdata_xmax2");
    dr = r-rc;
    /* adjust Omega according to Kepler's law 
       Om = sqrt(M/r^3)  ==> dOm = -(3/2) Om/r * dr. Here r=2*bc dr=2dx */
    dOm = -1.5*(Getd("BNSdata_Omega")/rc)*dr;
    Setd("BNSdata_Omega", Getd("BNSdata_Omega")+dOm);

    printf("Setting:\n");
    printf(" BNSdata_b = %s\n", Gets("BNSdata_b"));
    printf(" BNSdata_xmax1 = %s\n", Gets("BNSdata_xmax1"));
    printf(" BNSdata_xmax2 = %s\n", Gets("BNSdata_xmax2"));
    printf(" BNSdata_Omega = %s\n", Gets("BNSdata_Omega"));
    printf(" BNSdata_x_CM = %s\n", Gets("BNSdata_x_CM"));
    /* adjust coords such that all is ok, e.g. box5/4 need to cover
       holes in box0/3 */
    BNSgrid_init_Coords(grid);
    /* set wB */
    BNS_set_wB(grid, 1, Getd("BNSdata_xmax1"),0.0,0.0);
    BNS_set_wB(grid, 2, Getd("BNSdata_xmax2"),0.0,0.0); 
    /* Ok now we have complete data, but the var 
       Coordinates_AnsorgNS_sigma_pm should probabaly be changed as well! */

    /* pick Coordinates_AnsorgNS_sigma_pm such that the masses are what we
       want */
    /***********************************************************************************************/
    /* do newton_linesrch_itsP iterations until sigma+ and sigma- are such that masses are correct */
    /***********************************************************************************************/
    printf("Rescaling Coordinates_AnsorgNS_sigma_pm so that masses are correct:\n");
    /* star1 */
    pars->grid = grid;
    pars->b = 0;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    BNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, facvec[1], 0);
    BNSgrid_init_Coords(grid);
    /* set wB */
    BNS_set_wB(grid, 1, Getd("BNSdata_xmax1"),0.0,0.0);

    /* star2 */
    pars->grid = grid;
    pars->b = 3;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m02);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    BNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, facvec[1], 3);
    BNSgrid_init_Coords(grid);
    /* set wB */
    BNS_set_wB(grid, 2, Getd("BNSdata_xmax2"),0.0,0.0); 
  }

  /* check the vars we set for NANs */
  vlpush(varlist, Ind("BNSdata_wBx"));
  printf("checking for NANs in:\n");
  for(j = 0; j < varlist->n; j++) printf("%s ", VarName(varlist->index[j]));
  printf("\n");  
  NumberChecker_CheckIfFinite_VarList(varlist);
  fflush(stdout);

  /* free varlist */
  vlfree(varlist);
  
  /* free saved strings */
  free(BNSdata_b_sav);
  free(BNSdata_x_CM_sav);
  free(BNSdata_xmax1_sav);
  free(BNSdata_xmax2_sav);
}

/************************************************************************/
/* utilities to manipulate the grid */
/************************************************************************/
/* given a new 
   (Coordinates_AnsorgNS_sigma_pm, 
    Coordinates_AnsorgNS_dsigma_pm_dB, Coordinates_AnsorgNS_dsigma_pm_dphi)
   on one side of grid initialize Coordinates.
   Do both sides if innerdom<0 */
void BNSgrid_init_Coords_pm(tGrid *grid, int innerdom)
{
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  /* avoid too much printf */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");

  /* initialize coords on grid */
  init_CoordTransform_And_Derivs(grid);

  /* reset box5/4 boundaries so that A=Amax in box0/3 will be inside box5/4 */
  adjust_box4_5_pars_pm(grid, innerdom);
  set_BoxStructures_fromPars(grid, 0);

  /* reset x,y,z, dXdx and such */
  init_CoordTransform_And_Derivs(grid);

  /* set values of A,B,phi in box4/5 */
  set_BNSdata_ABphi(grid);

  /* put back original Coordinates_verbose */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");
}
/* given a new 
   (Coordinates_AnsorgNS_sigma_pm, 
    Coordinates_AnsorgNS_dsigma_pm_dB, Coordinates_AnsorgNS_dsigma_pm_dphi)
   initialize Coordinates on both sides of grid */
void BNSgrid_init_Coords(tGrid *grid)
{
  BNSgrid_init_Coords_pm(grid, -1); /* init both sides of grid */
}



/* scale Coordinates_AnsorgNS_sigma and its deriv on one side
   by a factor fac. For star1: ibd=0 , for star2 ibd=3 */
void BNSgrid_scale_Coordinates_AnsorgNS_sigma(tGrid *grid, double fac, int ibd)
{
  int b,i;
  int isigma      = Ind("Coordinates_AnsorgNS_sigma_pm");
  int isigma_dB   = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  int isigma_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  int obd; /* outer box */

  if(ibd==3)      obd = 2;
  else if(ibd==0) obd = 1;
  else errorexit("BNSgrid_scale_Coordinates_AnsorgNS_sigma: set ibd=3 or ibd=0.");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *sigma      =  box->v[isigma];
    double *sigma_dB   =  box->v[isigma_dB];
    double *sigma_dphi =  box->v[isigma_dphi];

    if(b!=ibd && b!=obd) continue;
    forallpoints(box, i)
    {
      sigma[i] *= fac;
      sigma_dB[i] *= fac;
      sigma_dphi[i] *= fac;
    }
  }
}

/*****************************************************************/
/* useful funcs for equal mass symmetry                           */
/*****************************************************************/

/* make DomainShape on left equal to the one on the right:
   if we copy from right into left, set inner box destination ibd=3 
   if we copy from left into right, set inner box destination ibd=0 */
void BNSgrid_copy_DomainShape(tGrid *grid, int ibd)
{
  int isigma      = Ind("Coordinates_AnsorgNS_sigma_pm");
  int isigma_dB   = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");
  int isigma_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dphi");
  int i,j,k, kd, ijks, ijkd;
  int n1 = grid->box[ibd]->n1;
  int n2 = grid->box[ibd]->n2;
  int n3 = grid->box[ibd]->n3;
  double sigp_Bphi;
  int obd; /* outer box on destination */
  int ibs; /* inner box on source */
  int obs; /* outer box on source */

errorexit("BNSgrid_copy_DomainShape is untested. But you can take out "
          "this line to try it!");

  if(ibd==3)      { obd = 2;  ibs = 0;  obs = 1; }
  else if(ibd==0) { obd = 1;  ibs = 3;  obs = 2; }
  else errorexit("BNSgrid_copy_DomainShape: set ibd=3 or ibd=0.");
  printf("BNSgrid_copy_DomainShape:  from: box%d/%d  to: box%d/%d\n",
         ibs,obs, ibd,obd);
  if(n3 % 2)
    errorexit("BNSgrid_copy_DomainShape: box[0]->n3 has to be divisible by 2.");

  /* use the fact that (A,B,phi)    on right 
     corresponds to    (A,B,Pi-phi) on the left */
  /* set Coordinates_AnsorgNS_sigma_pm on dest. side from 
         Coordinates_AnsorgNS_sigma_pm on source side     */ 
  for(k=1; k<n3; k++) /* assume n1,n2,n3 are equal on left and right */
  {
    if(k<=n3/2) kd=n3/2 - k;
    else        kd=3*n3/2 - k;
    for(j=0; j<n2; j++)
    for(i=0; i<n1; i++)
    {
      ijks = Index(i,j,k);
      ijkd = Index(i,j,kd);
      sigp_Bphi = grid->box[ibs]->v[isigma][ijks];
      grid->box[ibd]->v[isigma][ijkd] = -sigp_Bphi;
      grid->box[obd]->v[isigma][ijkd] = -sigp_Bphi;
    }
  }

  /* compute derivs of sigma in both domains */
  spec_Deriv1(grid->box[ibd], 2, grid->box[ibd]->v[isigma],
              grid->box[ibd]->v[isigma_dB]);
  spec_Deriv1(grid->box[ibd], 3, grid->box[ibd]->v[isigma],
              grid->box[ibd]->v[isigma_dphi]);
  spec_Deriv1(grid->box[obd], 2, grid->box[obd]->v[isigma],
              grid->box[obd]->v[isigma_dB]);
  spec_Deriv1(grid->box[obd], 3, grid->box[obd]->v[isigma],
              grid->box[obd]->v[isigma_dphi]);

  /* initialize coords */
  BNSgrid_init_Coords(grid);
}

/* copy a var between right and left in case of equal masses:
   if we copy from right into left, set inner box destination ibd=3
   if we copy from left into right, set inner box destination ibd=0
   we should call BNSgrid_copy_DomainShape before */
void BNSgrid_set_Var_equalmasses_sym(tGrid *grid, int ibd, int iv, int sym)
{
  int i,j,k, kd, ijks, ijkd;
  int n1 = grid->box[ibd]->n1;
  int n2 = grid->box[ibd]->n2;
  int n3 = grid->box[ibd]->n3;
  double var_in, var_out;
  int obd; /* outer box on destination */
  int ibs; /* inner box on source */
  int obs; /* outer box on source */

errorexit("BNSgrid_set_Var_equalmasses_sym is untested. But you can take out "
          "this line to try it!");

  if(ibd==3)      { obd = 2;  ibs = 0;  obs = 1; }
  else if(ibd==0) { obd = 1;  ibs = 3;  obs = 2; }
  else errorexit("BNSgrid_set_Var_equalmasses_sym: set ibd=3 or ibd=0.");
  if(n3 % 2)
    errorexit("BNSgrid_copy_DomainShape: box[0/1/2/3]->n3 has to be divisible by 2.");
  //  BNSgrid_copy_DomainShape(grid, 3);

  printf("BNSgrid_set_Var_equalmasses_sym: %s from: box%d/%d to: box%d/%d  "
         "sym=%d\n", VarName(iv), ibs,obs, ibd,obd, sym);
        
  /* use the fact that (A,B,phi)    on right 
     corresponds to    (A,B,Pi-phi) on the left */
  /* set var on left from var on right */
  for(k=1; k<n3; k++) /* assume n1,n2,n3 are equal on left and right */
  {
    if(k<=n3/2) kd=n3/2 - k;
    else        kd=3*n3/2 - k;
    for(j=0; j<n2; j++)
    for(i=0; i<n1; i++)
    {
      ijks = Index(i,j,k);
      ijkd = Index(i,j,kd);
      var_in = grid->box[ibs]->v[iv][ijks];
      var_out= grid->box[obs]->v[iv][ijks];
      grid->box[ibd]->v[iv][ijkd] = var_in*sym;
      grid->box[obd]->v[iv][ijkd] = var_out*sym;
    }
  }
}

/* copy all vars from right to left for equal masses */
void BNSgrid_set_allVars_onLeft_equalmasses(tGrid *grid)
{
  printf("BNSgrid_set_allVars_onLeft_equalmasses...\n");

errorexit("BNSgrid_set_allVars_onLeft_equalmasses is untested. "
          "I'm not sure about the syms. But you can take out "
          "this line to try it!");

  /* adjust Domain Shape on left */
  BNSgrid_copy_DomainShape(grid, 3);

  /* set scalars on left */
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_Psi"), 1);
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_alphaP"), 1);
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_Sigma"), 1);
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_q"), 1);
  /* set vectors on left */
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_Bx"), -1);
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_By"), -1);
  BNSgrid_set_Var_equalmasses_sym(grid, 3, Ind("BNSdata_Bz"), 1);
}
