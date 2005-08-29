/* sgrid_ScalarWave.c */
/* Wolfgang Tichy  8/2005 */

#include "sgrid.h"
#include "ScalarWave.h"




int sgrid_ScalarWave(void) 
{
  if (!Getv("physics", "ScalarWave")) return 0;
  printf("Adding ScalarWave\n");

  /* functions */
  AddFun(POST_INITIALDATA, ScalarWave_startup, 
	 "initialize ScalarWave system from adm initial data");
  AddFun(ANALYZE, ScalarWave_analyze, "compute ScalarWave energy density");

  /* variables */
  AddVar("ScalarWave_psi",    "", "scalar");
  AddVar("ScalarWave_psidot", "", "time deriv of scalar");

  /* derivatives which need to be precomputed before each evo steo */
  /* Note: ADMvars_ddg, ADMvars_dg, ADMvars_dA are overwritten and used 
           in place of ScalarWave_ddg, ScalarWave_dg, ScalarWave_dA */
  AddVar("ScalarWave_dpsi",   "i",    "1st deriv of scalar");
  AddVar("ScalarWave_ddpsi",  "(ij)", "2nd deriv of scalar");

  /* energy density of psi */
  AddVar("ScalarWave_rho",      "", "energy density of scalar");
  AddVar("ScalarWave_2dInt_rho","", "surface integral of energy density");
  AddVar("ScalarWave_temp",     "", "temp var");
    
  /* parameters */
  AddPar("ScalarWave_useDD", "yes",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("ScalarWave_reset_doubleCoveredPoints", "no",
         "whether we reset double covered points after each evo step [no,yes]");
  AddPar("ScalarWave_coordinateDependentFilter", "no",
         "whether coordinate dependent filters (e.g. for SphereicalDF) "
         "are used [no,yes]");

  AddPar("ScalarWave_A",    "1", "amplitude of wave");
  AddPar("ScalarWave_sigma","1", "sigma of Gaussian wavepacket");
  AddPar("ScalarWave_r0",   "1", "location of spherical Gaussian wavepacket");
  AddPar("ScalarWave_x0",   "1", "x-location Gaussian wavepacket");
  AddPar("ScalarWave_y0",   "0", "y-location Gaussian wavepacket");
  AddPar("ScalarWave_z0",   "0", "z-location Gaussian wavepacket");
  AddPar("ScalarWave_waveform", "spherical", "form a wave");
  AddPar("ScalarWave_nonlinearity",  "0", "non-linear potential");
         
  return 0;
}
