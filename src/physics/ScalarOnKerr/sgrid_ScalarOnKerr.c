/* sgrid_ScalarOnKerr.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "ScalarOnKerr.h"




int sgrid_ScalarOnKerr(void) 
{
  if (!Getv("physics", "ScalarOnKerr")) return 0;
  printf("Adding ScalarOnKerr\n");

  /* functions */
  AddFun(PRE_GRID, ScalarOnKerr_setup_boxes, "setup initial box sizes");
  AddFun(POST_INITIALDATA, ScalarOnKerr_startup, 
	 "initialize ScalarOnKerr system from adm initial data");
  AddFun(ANALYZE, ScalarOnKerr_analyze, "compute something useful");

  /* variables */
  AddVar("ScalarOnKerr_psi", "", "scalar");
  AddVar("ScalarOnKerr_Pi",  "", "time deriv of scalar");
  AddVar("ScalarOnKerr_rho", "", "Source");

  /* char. variables */
  AddVar("ScalarOnKerr_Up",  "", "characteristic var. U_{+}");
  AddVar("ScalarOnKerr_Um",  "", "characteristic var. U_{-}");

  /* derivatives which need to be precomputed before each evo step */
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
  AddConstantVar("ScalarOnKerr3d_G",     "I",     "G^I = Gamma^I_jk g^jk");
  AddConstantVar("ScalarOnKerr3d_dalpha","i",     "1st spatial deriv of lapse");
  AddConstantVar("ScalarOnKerr3d_dbeta","Ij",     "1st spatial deriv of shift");

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
  AddPar("ScalarOnKerr_sourcetype","Y22test",
         "source we use [Y22test,DV_source]");
  AddPar("ScalarOnKerr_reset_doubleCoveredPoints", "yes",
         "whether we reset double covered points in each evo substep [no,yes]");
//  AddPar("ScalarOnKerr_filter_unew", "simple",
//         "whether we filter all unew in each evo substep "
//         "[no, simple, naive_Ylm]");
  AddPar("ScalarOnKerr_filter_vars", "no",
         "which vars we filter in each evo substep [no,psi,Pi,phi]");
  AddPar("ScalarOnKerr_filter_type", "simple",
         "how we filter in each evo substep [simple,naive_Ylm,Ylm_lmshift,X2/2]");
  AddPar("ScalarOnKerr_filter_time", "afterRHS afterBC",
         "when we filter in evo substep [afterRHS,afterBC] or [POST_EVOLVE]");
  if(Getv("ScalarOnKerr_filter_time", "POST_EVOLVE"))
    AddFun(POST_EVOLVE, ScalarOnKerr_filter, "filter in POST_EVOLVE");
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
//  AddPar("ScalarOnKerr_special_nPi_filter", "simple",
//         "whether we filter all new Pi in each evo substep "
//         "[no, simple, naive_Ylm]");
//  AddPar("ScalarOnKerr_special_nPi_nphi_filter", "naive_Ylm",
//         "whether we filter all new Pi and phi_i in each evo substep "
//         "[no, simple, naive_Ylm]");
  AddPar("ScalarOnKerr_overlap_shells", "no",
         "whether we use overlapping shells [no,yes]");
  AddPar("ScalarOnKerr_Pi_def", "psidot",
         "def. of Pi [psidot,ScheelsPi]");
  AddPar("ScalarOnKerr_1stOrder_inSpace", "no",
         "introduce extra vars to make system 1st order in space [no,yes]");
  if(Getv("ScalarOnKerr_1stOrder_inSpace", "yes"))
  {
    /* add phix=dpsi/dx, phiy=dpsi/dy, phiz=dpsi/dz, for 1st order reduction */
    AddVar("ScalarOnKerr_phi",  "i", "1st spatial derivs of scalar");
    /* add 0 speed char. vars: phim =m^j phix^j, phil =l^j phix^j 
       [n.m=0, n.l=0, m.l=0, n.n=1, m.m=1, l.l=1] */
    AddVar("ScalarOnKerr_U0",  "i", "transverse char. var. U0");
    AddVar("ScalarOnKerr_dphi", "ij", "spatial derivs of phix,...");
  }
  AddPar("ScalarOnKerr_constraints", "no", "compute constraints [no,yes]");
  if(Getv("ScalarOnKerr_constraints", "yes") && 
     Getv("ScalarOnKerr_1stOrder_inSpace", "yes"))
  {
    AddVar("ScalarOnKerr_C",  "i", "Cx=phix - d/dx psi");
  }
  AddPar("ScalarOnKerr_modes", "no", "compute modes [no,yes]");
  if(Getv("ScalarOnKerr_modes", "yes"))
  {
    AddVar("ScalarOnKerr_rhomodes",  "",  "coeffs of rho");
    AddVar("ScalarOnKerr_psimodes",  "",  "coeffs of psi");
    AddVar("ScalarOnKerr_Pimodes",   "",  "coeffs of Pi");
    if(Getv("ScalarOnKerr_1stOrder_inSpace", "yes"))
      AddVar("ScalarOnKerr_phimodes",  "i", "coeffs of phi_i");
  }

  return 0;
}
