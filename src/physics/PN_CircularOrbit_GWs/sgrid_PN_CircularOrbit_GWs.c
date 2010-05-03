/* sgrid_PN_CircularOrbit_GWs.c */
/* Wolfgang Tichy, 4/2010 */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"


int sgrid_PN_CircularOrbit_GWs() 
{
  if (!Getv("physics", "PN_CircularOrbit_GWs")) return 0;
  printf("Adding PN_CircularOrbit_GWs\n");

  /* functions: */
  AddFun(PRE_GRID, PN_CircularOrbit_GWs_set_boxsize, "setup box size");
  AddFun(INITIALDATA, PN_CircularOrbit_GWs_startup, "startup");
  AddFun(ANALYZE, PN_CircularOrbit_GWs, "PN_CircularOrbit_GWs: read data and compute modes");

  /* variables */
  AddVar("PN_CircularOrbit_GWs_hplus" ,  "", "h+ on a sphere, note H = (h+ - ihx)");
  AddVar("PN_CircularOrbit_GWs_hcross",  "", "hx on a sphere, note H = (h+ - ihx)");
  AddVar("PN_CircularOrbit_GWs_Re_Hmode","", "real part of H = (h+ - ihx) mode");
  AddVar("PN_CircularOrbit_GWs_Im_Hmode","", "imag part of H = (h+ - ihx) mode");
  AddVar("PN_CircularOrbit_GWs_Re_Psi4",     "", "real part of Psi4");
  AddVar("PN_CircularOrbit_GWs_Im_Psi4",     "", "imag part of Psi4");
  AddVar("PN_CircularOrbit_GWs_Re_Psi4mode", "", "real part of Psi4 mode");
  AddVar("PN_CircularOrbit_GWs_Im_Psi4mode", "", "imag part of Psi4 mode");
  AddVar("PN_CircularOrbit_GWs_Re_sYlm", "", "real part of spin-weighted Ylm");
  AddVar("PN_CircularOrbit_GWs_Im_sYlm", "", "imaginary part spin-weighted Ylm");

  /* Parameters: */
  AddPar("PN_CircularOrbit_GWs_hfile_prefix", "h_",
         "prefix for GW file containing h");
  AddPar("PN_CircularOrbit_GWs_computePsi4", "no", "if we compute Psi4 [no,yes]");
  AddPar("PN_CircularOrbit_GWs_Psi4file_prefix", "psi4_",
         "prefix for GW file containing Psi4");
  AddPar("PN_CircularOrbit_GWs_orbitfile", "orbit.t",
         "imaginary part of data on sphere");
  AddPar("PN_CircularOrbit_GWs_lmax", "6", "max l mode we compute");
  AddPar("PN_CircularOrbit_GWs_sphere_Lmax", "13",
         "max l mode we can represent on grid. We need Lmax>=lmax");
  AddPar("PN_CircularOrbit_GWs_spinweight", "-2", "spin weight we use in sYlm");
  AddPar("PN_CircularOrbit_GWs_HmodeOutputFormat", "plus_cross",
         "Output format for Hmodes [plus_cross,Re_Im] (NINJA uses plus_cross)");

  AddPar("PN_CircularOrbit_GWs_match_NR", "no", "if we match NR Psi4 [no,yes]");
  if(Getv("PN_CircularOrbit_GWs_match_NR", "yes"))
  {
    AddFun(POST_INITIALDATA, minimize_PN_NR_diff,
           "minimize_PN_NR_diff: minimize PN-NR difference");
    
    AddVar("PN_CircularOrbit_GWs_Re_NRPsi4", "", "real part of NR Psi4");
    AddVar("PN_CircularOrbit_GWs_Im_NRPsi4", "", "imag part of NR Psi4");
    AddVar("PN_CircularOrbit_GWs_PN_NR_diff", "", "PN - NR");

    AddPar("PN_CircularOrbit_GWs_NRPsi4file_prefix", "psi4_",
           "prefix for file containing numerical Psi4 modes");
  }

  AddPar("PN_CircularOrbit_GWs_omega", "0.01", "initial orbital frequency * m");
  AddPar("PN_CircularOrbit_GWs_m1", "0.5", "mass of particle 1");
  AddPar("PN_CircularOrbit_GWs_m2", "0.5", "mass of particle 2");
  AddPar("PN_CircularOrbit_GWs_chi1x", "0", "initial S1x/m1^2");
  AddPar("PN_CircularOrbit_GWs_chi1y", "0", "initial S1y/m1^2");
  AddPar("PN_CircularOrbit_GWs_chi1z", "0", "initial S1z/m1^2");
  AddPar("PN_CircularOrbit_GWs_chi2x", "0", "initial S2x/m2^2");
  AddPar("PN_CircularOrbit_GWs_chi2y", "0", "initial S2y/m2^2");
  AddPar("PN_CircularOrbit_GWs_chi2z", "0", "initial S2z/m2^2");
  AddPar("PN_CircularOrbit_GWs_Lnx", "0", "initial L_hat_Newton_x");
  AddPar("PN_CircularOrbit_GWs_Lny", "0", "initial L_hat_Newton_y");
  AddPar("PN_CircularOrbit_GWs_Lnz", "1", "initial L_hat_Newton_z");
  AddPar("PN_CircularOrbit_GWs_Phi", "0", "initial orbital phase Phi");
  AddPar("PN_CircularOrbit_GWs_D", "1.26", "distance from source");
  AddPar("PN_CircularOrbit_GWs_t1", "0",    "initial time");
  AddPar("PN_CircularOrbit_GWs_t2", "1000", "final time");
  AddPar("PN_CircularOrbit_GWs_dt", "10",   "time step");               

  return 0;
}
