/* sgrid_ScalarOnKerr.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "ScalarOnKerr.h"




int sgrid_ScalarOnKerr(void) 
{
  if (!Getv("physics", "ScalarOnKerr")) return 0;
  printf("Adding ScalarOnKerr\n");

  /* functions */
  AddFun(POST_INITIALDATA, ScalarOnKerr_startup, 
	 "initialize ScalarOnKerr system from adm initial data");
  AddFun(ANALYZE, ScalarOnKerr_analyze, "compute something useful");

  /* variables */
  AddVar("ScalarOnKerr_psi", "", "scalar");
  AddVar("ScalarOnKerr_pi",  "", "time deriv of scalar");

  /* derivatives which need to be precomputed before each evo steo */
  AddVar("ScalarOnKerr_dpsi",  "i",    "1st spatial deriv of scalar");
  AddVar("ScalarOnKerr_ddpsi", "(ij)", "2nd spatial deriv of scalar");
  AddVar("ScalarOnKerr_dpi",   "i",    "1st spatial deriv of pi");

  /* Kerr background 4-metric */
  AddConstantVar("ScalarOnKerr_g",     "(ab)",  "Kerr metric");
  AddConstantVar("ScalarOnKerr_gup",   "(AB)",  "inverse Kerr metric");
  AddConstantVar("ScalarOnKerr_Gamma", "A(bc)", "Christoffel symbol of Kerr metric");

  /* Kerr background in 3+1 split */
  AddConstantVar("ScalarOnKerr3d_g",     "(ij)",  "Kerr 3-metric");
  AddConstantVar("ScalarOnKerr3d_alpha", "",  	  "Kerr lapse");
  AddConstantVar("ScalarOnKerr3d_beta",  "I",     "Kerr shift");
  AddConstantVar("ScalarOnKerr3d_K",     "(ij)",  "extrinsic curv.");
  AddConstantVar("ScalarOnKerr3d_TrK",   "",      "trace of K");
  AddConstantVar("ScalarOnKerr3d_gup",   "(IJ)",  "inverse Kerr 3-metric");
  AddConstantVar("ScalarOnKerr3d_Gamma", "I(jk)", "Christoffel symbol of 3-metric");
  AddConstantVar("ScalarOnKerr3d_dalpha","i",     "1st spatial deriv of lapse");

  if(!Getv("physics", "ADMvars"))
  {
    AddVar("temp1",    "", "temporary storage 1");
    AddVar("temp2",    "", "temporary storage 2");
    AddVar("temp3",    "", "temporary storage 3");
  }

  /* parameters */
  AddPar("BHmass", "1.0", "mass of black hole");
  AddPar("BHsx", "0.0", "spin_x of black hole");
  AddPar("BHsy", "0.0", "spin_y of black hole");
  AddPar("BHsz", "0.0", "spin_z of black hole");
  AddPar("ScalarOnKerr_reset_doubleCoveredPoints", "yes",
         "whether we reset double covered points after each evo step [no,yes]");
         
  return 0;
}
