/* sgrid_TestDerivs.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "TestDerivs.h"


int sgrid_TestDerivs() 
{
  if (!Getv("physics", "TestDerivs")) return 0;
  printf("Adding TestDerivs\n");

  /* functions */
  AddFun(PRE_COORDINATES, Set_Test_Coordinates_AnsorgNS_sigma_pm,
         "setup the sigma_{+-} vars from AnsorgNS");
  AddFun(PRE_COORDINATES, Setup_CubedSpheres_initboxes,
         "setup some cubed spheres");
  AddFun(INITIALDATA, TestDerivs_startup, "initialize test");
  AddFun(ANALYZE, TestDerivs_analyze, "compute test results");

  /* variables */
  AddVar("TestDerivs_u",       "",     "test function");
  AddVar("TestDerivs_Err_du",  "i",    "Error in 1st deriv of test function");
  AddVar("TestDerivs_Err_ddu", "(ij)", "Error in 2nd deriv of test function");
   
  /* parameters */
  AddPar("TestDerivs_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("TestDerivs_A",     "1",   "amplitude of wave");
  AddPar("TestDerivs_sigmax","1",   "sigmax of Gaussian wavepacket");
  AddPar("TestDerivs_sigmay","1.3", "sigmay of Gaussian wavepacket");
  AddPar("TestDerivs_sigmaz","0.7", "sigmaz of Gaussian wavepacket");
  AddPar("TestDerivs_x0",    "1",   "x-location of Gaussian wavepacket");
  AddPar("TestDerivs_y0",    "0.3", "y-location of Gaussian wavepacket");
  AddPar("TestDerivs_z0",    "0.1", "z-location of Gaussian wavepacket");
  AddPar("TestDerivs_grid",  "", "grid we use [CubedSpheres]");
	     	   	   	 
  return 0;
}
