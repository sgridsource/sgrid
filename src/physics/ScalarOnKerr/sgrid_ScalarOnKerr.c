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
  AddVar("ScalarOnKerr_Pi",  "", "time deriv of scalar");
  AddVar("ScalarOnKerr_rho", "", "Source");

  /* derivatives which need to be precomputed before each evo steo */
  AddVar("ScalarOnKerr_dpsi",  "i",    "1st spatial deriv of scalar");
  AddVar("ScalarOnKerr_ddpsi", "(ij)", "2nd spatial deriv of scalar");
  AddVar("ScalarOnKerr_dPi",   "i",    "1st spatial deriv of Pi");

  /* Kerr background 4-metric */
  AddConstantVar("ScalarOnKerr_g",     "(ab)",  "Kerr metric");
  AddConstantVar("ScalarOnKerr_gup",   "(AB)",  "inverse Kerr metric");
  AddConstantVar("ScalarOnKerr_Gamma", "A(bc)", "Christoffel symbol of Kerr metric");
  AddConstantVar("ScalarOnKerr_G",     "A",     "G^A = Gamma^A_bc g^bc");

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
         "whether we reset double covered points in each evo substep [no,yes]");
  AddPar("ScalarOnKerr_filter_unew", "simple",
         "whether we filter all unew in each evo substep "
         "[no, simple, naive_Ylm]");
  AddPar("ScalarOnKerr_filter_n2frac", "0.66666666667",
         "fraction of the n2 coeffs to keep, when filtering Y-direc.");
  AddPar("ScalarOnKerr_filter_n3frac", "0.66666666667",
         "fraction of the n3 coeffs to keep, when filtering Z-direc.");
  AddPar("ScalarOnKerr_filter_shift2", "0",
         "shift index of last coeff to keep, when filtering Y-direc.");
  AddPar("ScalarOnKerr_filter_shift3", "0",
         "shift index of last coeff to keep, when filtering Z-direc.");
  AddPar("ScalarOnKerr_filter_YZregion", "square",
         "region outside which we apply filters [ellipse,square]");
  AddPar("ScalarOnKerr_special_nPi_filter", "simple",
         "whether we filter all new Pi in each evo substep "
         "[no, simple, naive_Ylm]");

  return 0;
}
