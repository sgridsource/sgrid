/* sgrid_BNSdata.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "BNSdata.h"


int sgrid_BNSdata() 
{
  if (!Getv("physics", "BNSdata")) return 0;
  printf("Adding BNSdata\n");

  /* functions */
  AddFun(PRE_INITIALDATA, BNSdata_startup, "initialize BNSdata");
  AddFun(INITIALDATA, BNSdata_solve, "solve BNSdata Eq.");
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

  AddVar("BNSdata_temp1", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp2", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp3", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp4", "", "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("BNSdata_n",     "1", "polytropic index n, Gamma = 1 + 1/n");
  AddPar("BNSdata_kappa", "1", "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("BNSdata_Omega", "0", "orbital angular velocity");
  AddPar("BNSdata_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("BNSdata_itmax", "10", "maximal number of Newton iterations");
  AddPar("BNSdata_tol",   "1e-6","tolerance for Newton step");
  AddPar("BNSdata_linSolver", "UMFPACK", 
         "linear solver used [LAPACK,templates_GMRES,bicgstab,UMFPACK]");
  AddPar("BNSdata_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("BNSdata_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_guess", "test", "initial guess [test]");
         
  AddPar("BNSdata_grid", "SphericalDF",
         "what grid we use [SphericalDF, AnsorgNS, 4ABphi_2xyz]");
  AddPar("BNSdata_regularization", "none",
         "options for 4ABphi_2xyz and AnsorgNS " 
         "[regularity_on_axis,regularity_on_axis_at_center]");
	     	   	   	 
  return 0;
}
