/* sgrid_BNSdata.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "BNSdata.h"


int sgrid_BNSdata() 
{
  if (!Getv("physics", "BNSdata")) return 0;
  printf("Adding BNSdata\n");

  /* functions */
  AddFun(PRE_GRID, set_boxsizes, "setup initial box sizes");
  AddFun(POST_GRID, set_sigma_pm_vars, "setup the sigma_{+-} vars from AnsorgNS");
  AddFun(PRE_INITIALDATA, BNSdata_startup, "initialize BNSdata");
  AddFun(INITIALDATA, setBNSdata, "set the BNS data");
  AddFun(ANALYZE, BNSdata_analyze, "compute properties of BNS data");

  /* variables */
  AddVar("BNSdata_Psi",     "",     "new conf. factor");
  AddVar("BNSdata_Psi",     "i",    "1st deriv of Psi");
  AddVar("BNSdata_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("BNSdata_B",       "I",    "beta - Omega cross r");
  AddVar("BNSdata_B",       "Ij",   "1st deriv of B^i");
  AddVar("BNSdata_B",       "I(jk)","2nd deriv of B^i");
  AddVar("BNSdata_alphaP",  "",     "lapse times Psi");
  AddVar("BNSdata_alphaP",  "i",    "1st deriv of alphaP");
  AddVar("BNSdata_alphaP",  "(ij)", "2nd deriv of alphaP");
  AddVar("BNSdata_Sigma",   "",     "Sigma");
  AddVar("BNSdata_Sigma",   "i",    "1st deriv of Sigma");
  AddVar("BNSdata_Sigma",   "(ij)", "2nd deriv of Sigma");

  AddVar("BNSdata_vRS", "I", "solenoidal velocity in rotating frame");
  AddVar("BNSdata_q",   "",  "q := P/rho0");
  AddVar("BNSdata_vRS", "Ij","1st deriv of vRS");
  AddVar("BNSdata_q",   "i", "1st deriv of q");

  /* sometimes we save the old vars before the ell. solve */
  AddVar("BNSdata_Psiold",    "",  "old Psi");
  AddVar("BNSdata_Bold",      "I", "old B");
  AddVar("BNSdata_alphaPold", "",  "old alphaP");
  AddVar("BNSdata_Sigmaold",  "",  "old Sigma");
  AddVar("BNSdata_qold",      "",  "old q");

  AddVar("BNSdata_A", "", "store value of A in box4/5");
  AddVar("BNSdata_B", "", "store value of B in box4/5");
  AddVar("BNSdata_phi", "", "store value of phi in box4/5");

  AddVar("BNSdata_temp1", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp2", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp3", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp4", "", "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("BNSdata_m01",   "0.141202", "rest mass of NS1");
  AddPar("BNSdata_m02",   "0.141202", "rest mass of NS2");
  AddPar("BNSdata_iterate_m0", "no", "whether we iterate rest masses [no,yes]");
  if(Getv("BNSdata_iterate_m0", "yes"))
  {
    AddPar("BNSdata_init_m01", "-1", "initial rest mass of NS1, if positive "
           "BNSdata_m01 will be changed during the iterations until we reach "
           "BNSdata_desired_m01");
    AddPar("BNSdata_init_m02", "-1", "initial rest mass of NS2, if positive "
           "BNSdata_m02 will be changed during the iterations until we reach "
           "BNSdata_desired_m02");
    AddPar("BNSdata_m0change", "0.1", "amount by which change BNSdata_m01/2 "
           "during iterations");
    AddPar("BNSdata_desired_m01", "-1", "desired rest mass1 if we iterate it");
    AddPar("BNSdata_desired_m02", "-1", "desired rest mass2 if we iterate it");
  }
  AddPar("BNSdata_qmax1", "0", "max q of NS1 along x-axis");
  AddPar("BNSdata_qmax2", "0", "max q of NS2 along x-axis");
  AddPar("BNSdata_xmax1", "0", "pos. of max q of NS1 along x-axis");
  AddPar("BNSdata_xmax2", "0", "pos. of max q of NS2 along x-axis");
  AddPar("BNSdata_Omega", "estimate",  "orbital angular velocity");
  AddPar("BNSdata_b",     "1",  "separation parameter (distance~2b)");
  AddPar("BNSdata_n",     "1",  "polytropic index n, Gamma = 1 + 1/n");
  AddPar("BNSdata_kappa", "1",  "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("BNSdata_x_CM",  "estimate",  "center of mass in x-direction");
  AddPar("BNSdata_C1",    "-1", "C1 in q = (C1/F-1)/(n+1) "
         "[needs to be adjusted so that m01 stays the constant]");
  AddPar("BNSdata_C2",    "-1", "C2 in q = (C2/F-1)/(n+1) "
         "[needs to be adjusted so that m02 stays the constant]");
  AddPar("BNSdata_guess", "TOV", "init. guess for Psi, alphaP and "
         "q [TOV,TOVaverage,TOVproduct] for shift [TaniguchiShift]");
  if(Getd("BNSdata_m02")==0.0)
  {
    AddPar("BNSdata_yshift1", "0", "shift NS1 in y-direction for testing");
    AddPar("BNSdata_adjustdomain01", "yes", "if we adjust domainshapes "
           "after shift [yes,no]");
  }
  AddPar("BNSdata_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("BNSdata_itmax", "10", "max. number of iterations in BNSdata_solve");
  AddPar("BNSdata_tol",   "1e-6", "tolerance for BNSdata_solve");
  AddPar("BNSdata_Newton_tolFac", "0.1", "tol for Newton is "
         "Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac)");
  AddPar("BNSdata_esw",   "1", "ell. solve weight: after ell. solve new "
         "and old values are averaged according to: "
         "new = esw*new + (1-esw)*old");
  AddPar("BNSdata_esw1",  "1", "second weight esw1 is used if better");
  AddPar("BNSdata_allow_esw1_first_at", "-1", "first iteration when esw=esw1 "
         "will be tried and allowed if better [#,-1]. -1 means never");
  AddPar("BNSdata_adjust", "nothing", "what we adjust (apart from C1/2) "
         "after ell. solve. E.g. \"keep_xout keep_one_xout\" adjusts Omega "
         "and x_CM to keep either xout1 or xout2 in place. "
         "[nothing,WT_L2_method,keep_xout [keep_one_xout,always],"
         "keep_xmax [keep_one_xmax,always,reset_xmax_if_problem]]");
  AddPar("BNSdata_adjust_first_at", "0", 
         "first iteration when we use BNSdata_adjust. -1 means never");
  AddPar("BNSdata_adjust_mintol", "1e-10", "always use tol>=mintol in adjust, "
         "e.g. when BNSdata_adjust=keep_xout");
  AddPar("BNSdata_dOmega_fac", "0.1", "dOmega = Omega*dOmega_Fac");
  AddPar("BNSdata_dx_CM_fac",  "0.1", "dx_CM = b*dx_CM_fac");
  AddPar("BNSdata_EllSolver_method", "BNS_Eqn_Iterator",
         "how we solve for Psi,B^i,alphaP,Sigma "
         "[allatonce, BNS_Eqn_Iterator, sequence1, sequence2]");
  AddPar("BNSdata_linSolver", "UMFPACK", 
         "linear solver used [LAPACK,templates_GMRES,bicgstab,UMFPACK]");
  AddPar("BNSdata_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("BNSdata_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_guess", "test", "initial guess [test]");
         
  AddPar("BNSdata_grid", "4ABphi_2xyz",
         "what grid we use [SphericalDF, AnsorgNS, 4ABphi_2xyz]");
  AddPar("BNSdata_box0_Amax", "0.85", "max A we use in box0 [0...1]");
  AddPar("BNSdata_box3_Amax", "0.85", "max A we use in box3 [0...1]");
  AddPar("BNSdata_2xyz_n", "6", "if 4ABphi_2xyz: use n*n*n points in box4/5");
  AddPar("BNSdata_regularization", "none",
         "options for 4ABphi_2xyz and AnsorgNS " 
         "[regularity_on_axis,regularity_on_axis_at_center]");
	     	   	   	 
  return 0;
}
