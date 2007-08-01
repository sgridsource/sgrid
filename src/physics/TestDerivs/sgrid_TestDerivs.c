/* sgrid_TestDerivs.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "TestDerivs.h"


int sgrid_TestDerivs() 
{
  if (!Getv("physics", "TestDerivs")) return;
  printf("Adding TestDerivs\n");

  /* functions */
  AddFun(INITIALDATA, TestDerivs_startup, "initialize test");
  AddFun(ANALYZE, TestDerivs_analyze, "compute test results");

  /* variables */
  AddVar("TestDerivs_u",       "",     "test function");
  AddVar("TestDerivs_Err_du",  "i",    "Error in 1st deriv of test function");
  AddVar("TestDerivs_Err_ddu", "(ij)", "Error in 2nd deriv of test function");
   
  /* parameters */
  AddPar("TestDerivs_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("TestDerivs_A",     "1", "amplitude of wave");
  AddPar("TestDerivs_sigmax","1", "sigmax of Gaussian wavepacket");
  AddPar("TestDerivs_sigmay","1", "sigmay of Gaussian wavepacket");
  AddPar("TestDerivs_sigmaz","1", "sigmaz of Gaussian wavepacket");
  AddPar("TestDerivs_x0",    "1", "x-location of Gaussian wavepacket");
  AddPar("TestDerivs_y0",    "0", "y-location of Gaussian wavepacket");
  AddPar("TestDerivs_z0",    "0", "z-location of Gaussian wavepacket");
	     	   	   	 
  return 0;
}
