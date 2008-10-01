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
//  AddFun(ANALYZE, BNSdata_analyze, "compute error");

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
  AddVar("BNSdata_qold","",  "old q");

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
  AddPar("BNSdata_qmax1", "0", "max q of NS1 along x-axis");
  AddPar("BNSdata_qmax2", "0", "max q of NS2 along x-axis");
  AddPar("BNSdata_Omega", "0",  "orbital angular velocity");
  AddPar("BNSdata_b",     "1",  "separation parameter (distance~2b)");
  AddPar("BNSdata_n",     "1",  "polytropic index n, Gamma = 1 + 1/n");
  AddPar("BNSdata_kappa", "1",  "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("BNSdata_x_CM",  "0",  "center of mass in x-direction");
  AddPar("BNSdata_C1",    "-1", "C1 in q = (C1/F-1)/(n+1) "
         "[needs to be adjusted so that m01 stays the constant]");
  AddPar("BNSdata_C2",    "-1", "C2 in q = (C2/F-1)/(n+1) "
         "[needs to be adjusted so that m02 stays the constant]");
  AddPar("BNSdata_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("BNSdata_itmax", "10", "maximal number of Newton iterations");
  AddPar("BNSdata_tol",   "1e-6","tolerance for Newton step");
  AddPar("BNSdata_EllSolver_method", "sequential",
         "how we solve for Psi,B^i,alphaP,Sigma [allatonce, sequential]");
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
  AddPar("BNSdata_regularization", "none",
         "options for 4ABphi_2xyz and AnsorgNS " 
         "[regularity_on_axis,regularity_on_axis_at_center]");
	     	   	   	 
  return 0;
}
