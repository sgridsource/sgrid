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
         "filename for evolved orbital parameters");
  AddPar("PN_CircularOrbit_GWs_orbitfile_format", "S/m^2",
         "format used in orbitfile [S/m^2,yvec,Add_E_J,Add_Eb/mu_L/Mmu]");
  AddPar("PN_CircularOrbit_GWs_lmax", "6", "max l mode we compute");
  AddPar("PN_CircularOrbit_GWs_sphere_Lmax", "13",
         "max l mode we can represent on grid. We need Lmax>=lmax");
  AddPar("PN_CircularOrbit_GWs_spinweight", "-2", "spin weight we use in sYlm");
  AddPar("PN_CircularOrbit_GWs_HmodeOutputFormat", "plus_cross",
         "Output format for Hmodes [plus_cross,Re_Im] (NINJA uses plus_cross)");
  AddPar("PN_CircularOrbit_GWs_theta_h", "-1", "add h+,hx to orbitfile at "
         "theta_h,phi_h if theta_h>=0");
  AddPar("PN_CircularOrbit_GWs_phi_h",   "0", "add h+,hx to orbitfile at "
         "theta_h,phi_h if theta_h>=0");

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
    AddPar("PN_CircularOrbit_GWs_mintol", "1e-6", "tol for powell minimizer");
    AddPar("PN_CircularOrbit_GWs_min_p_format", "3",
           "parameter format [1,3,5,7,109,9,11,13] "
           "see cases in func_to_minimize_for_numrec(p)");
  }

  AddPar("PN_CircularOrbit_GWs_omega", "0.01", "initial orbital frequency * M");
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
  AddPar("PN_CircularOrbit_GWs_verbose", "yes", "verbose [yes,no]");

  /* pick EOM used */
  AddPar("PN_CircularOrbit_GWs_OrbitEOMtype", "BuonannoEtAl2003",
         "[Kidder1995,BuonannoEtAl2003,TaylorT4,Kidder1995_v2cut,"
         "BuonannoEtAl2003_v2cut]");

  /* flags for trajectories: if 1 include this PN term in EOM: */
  AddPar("PN_CircularOrbit_GWs_OrbitOv1", "1", "if 1 include O(v)^1 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv2", "1", "if 1 include O(v)^2 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv3", "1", "if 1 include O(v)^3 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv4", "1", "if 1 include O(v)^4 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv5", "1", "if 1 include O(v)^5 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv6", "1", "if 1 include O(v)^6 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOv7", "1", "if 1 include O(v)^7 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOLS1", "1", "if 1 include O(LS)^1 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOLS2", "1", "if 1 include O(LS)^2 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOLS3", "1", "if 1 include O(LS)^3 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_OrbitOSS1", "1", "if 1 include O(SS)^1 terms in EOM");
  AddPar("PN_CircularOrbit_GWs_Orbitv2cutoff", "0.6",
         "value of v^2 at which domega/dt=0");
  AddPar("PN_CircularOrbit_GWs_Orbitv6cutoff", "0.6",
         "value of v^6 at which domega/dt=0");
  
  /* flags for wave amplitudes: if 1 include this PN term in h_ij: */
  AddPar("PN_CircularOrbit_GWs_AmpOv1", "1", "if 1 include O(v)^1 in Amp");
  AddPar("PN_CircularOrbit_GWs_AmpOv2", "1", "if 1 include O(v)^2 in Amp");
  AddPar("PN_CircularOrbit_GWs_AmpOv3", "1", "if 1 include O(v)^3 in Amp");
  AddPar("PN_CircularOrbit_GWs_AmpOv4", "0", "if 1 include O(v)^4 in Amp");

  return 0;
}
