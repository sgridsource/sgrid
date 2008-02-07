/* grid_setup.c */
/* Wolfgang Tichy 2008 */


#include "sgrid.h"
#include "BNSdata.h"

#define Power pow


/* from Coordinates.c */
extern double  (*Coordinates_AnsorgNS_sigmap)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmap_dB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmap_dphi)(tBox *box, int ind, double B, double phi);
extern double  (*Coordinates_AnsorgNS_sigmam)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dphi)(tBox *box, int ind, double B, double phi);
extern double Coordinates_AnsorgNS_b; /* value of x if A=1 in AnsorgNS0/3 */

/* global vars in this file */
double sigma1; /* initial sigma_+ */
double sigma2; /* initial sigma_- */
double rf_surf1; /* radius of star1 */
double rf_surf2; /* radius of star2 */

/* funs in this file */
void rf_surf1_VectorFunc(int n, double *vec, double *fvec);
void rf_surf2_VectorFunc(int n, double *vec, double *fvec);




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
  double m1, Phic1, Psic1, m01;
  double m2, Phic2, Psic2, m02;
  double kappa     = Getd("BNSdata_kappa");
  double BNSdata_n = Getd("BNSdata_n");
  double BNSdata_b = Getd("BNSdata_b");
  double Gamma     = 1.0 + 1.0/BNSdata_n;
  double vec[2];
  int check;

  double Pc1=0.01; // change this: find Pc1 s.t. m01 = Getd("BNSdata_m01");
  double Pc2=0.01; // ???

  printf("set_boxsizes: setting box sizes and coordinates used ...\n");
  
  /* set Coordinates_AnsorgNS_sigmap/m funcs, .s.t. it works with sigma1/2 */
  Coordinates_AnsorgNS_sigmap       = return_sigma1;
  Coordinates_AnsorgNS_dsigmap_dB   = return_dsigma_zero;
  Coordinates_AnsorgNS_dsigmap_dphi = return_dsigma_zero;
  Coordinates_AnsorgNS_sigmam       = return_sigma2;
  Coordinates_AnsorgNS_dsigmam_dB   = return_dsigma_zero;
  Coordinates_AnsorgNS_dsigmam_dphi = return_dsigma_zero;
  Setd("Coordinates_AnsorgNS_b", BNSdata_b);
  Coordinates_AnsorgNS_b = BNSdata_b;

  /* set some box pars */
  if(Getv("BNSdata_grid", "4ABphi_2xyz"))
  {
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

  /* TOV_init yields m01 for a given Pc1 */
  TOV_init(Pc1, kappa, Gamma, &rf_surf1, &m1, &Phic1, &Psic1, &m01);
  sigma1 = 1.0; // rf_surf1; //check??? not true!!!


  printf(" rf_surf1=%g m1=%g Phic1=%g Psic1=%g m01=%g\n",
         rf_surf1, m1, Phic1, Psic1, m01);
{
double m,P,Phi,Psi,m0;
double rf=1.02822;
TOV_m_P_Phi_Psi_m0_OF_rf(rf, rf_surf1, kappa, Gamma,
                      Pc1, Phic1, Psic1,
                      &m, &P, &Phi, &Psi, &m0);
printf(" check rf=%g: m=%g P=%g Phi=%g Psi=%g m0=%g\n", rf,m,P,Phi,Psi,m0);
}
  TOV_init(Pc2, kappa, Gamma, &rf_surf2, &m2, &Phic2, &Psic2, &m02);
  sigma2 = rf_surf2; //check???
// set box sizes
// ...
//Sets("Coordinates_AnsorgNS_set_sigma_pm_pointers", "no");


  /* find sigma1, s.t. radius is rf_surf1 */
  vec[1] = 1.0/rf_surf1;   /* initial guess */
  /* do newton_lnsrch iterations: */
  newton_lnsrch(vec, 1, &check, rf_surf1_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  if(check) printf(": check=%d\n", check);  
  sigma1 = vec[1];
  printf(" setting: sigma1=%g\n", sigma1);

  /* find sigma2, s.t. radius is rf_surf2 */
  vec[1] = -1.0/rf_surf2;   /* initial guess */
  /* do newton_lnsrch iterations: */
  newton_lnsrch(vec, 1, &check, rf_surf2_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  if(check) printf(": check=%d\n", check);  
  sigma2 = vec[1];
  printf(" setting: sigma2=%g\n", sigma2);

//remove this
sigma1=1;
sigma2=-1;

  printf(" => radius of domain0 = %g,   radius of domain3 = %g\n", 
         0.5*(x_of_AnsorgNS0(NULL, -1, 0.0,0.0,0.0)-
              x_of_AnsorgNS0(NULL, -1, 0.0,1.0,0.0)),
         0.5*(x_of_AnsorgNS3(NULL, -1, 0.0,1.0,0.0)-
              x_of_AnsorgNS3(NULL, -1, 0.0,0.0,0.0)) );

  return 0;
}

/* funtion to be passed into newton_lnsrch by ??? */
void rf_surf1_VectorFunc(int n, double *vec, double *fvec)
{
  double xmax, xmin;
  
  sigma1 = vec[1];
  xmax = x_of_AnsorgNS0(NULL, -1, 0.0,0.0,0.0);
  xmin = x_of_AnsorgNS0(NULL, -1, 0.0,1.0,0.0);
  fvec[1] = (xmax-xmin)-2.0*rf_surf1;
}

/* funtion to be passed into newton_lnsrch by ??? */
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
  dsigma_pm_dphi = Ind("Coordinates_AnsorgNS_dsigma_pm_dB");

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
