/* sgrid_Poisson.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "Poisson.h"


int sgrid_Poisson() 
{
  if (!Getv("physics", "Poisson")) return 0;
  printf("Adding Poisson\n");

  /* functions */
  AddFun(PRE_INITIALDATA, Poisson_startup, "initialize Poisson");
  AddFun(INITIALDATA, Poisson_solve, "solve Poisson Eq.");
  AddFun(ANALYZE, Poisson_analyze, "compute error");

  /* variables */
  AddVar("Poisson_Psi",     "",
         "function that satisfies Poisson Eq: Laplace Psi = rh1");
  AddVar("Poisson_Psi",     "i",    "1st deriv of Psi");
  AddVar("Poisson_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("Poisson_rh1",     "",     "Source on RHS1");
  AddVar("Poisson_Err_Psi", "",     "Error in Psi");

  AddVar("Poisson_Chi",     "",
         "another function that satisfies Poisson: Laplace Chi = rh2");
  AddVar("Poisson_Chi",     "i",    "1st deriv of Chi");
  AddVar("Poisson_Chi",     "(ij)", "2nd deriv of Chi");
  AddVar("Poisson_rh2",     "",     "Source on RHS2");
  AddVar("Poisson_Err_Chi", "",     "Error in Chi");

  AddVar("Poisson_temp1", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson_temp2", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson_temp3", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson_temp4", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson_temp5", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson_temp6", "",   "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("Poisson_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("Poisson_itmax", "10", "maximal number of Newton iterations");
  AddPar("Poisson_tol",   "1e-6","tolerance for Newton step or multigrid");
//  AddPar("Poisson_linSolver", "bicgstab", 
//         "linear solver used [bicgstab, HYPRE, gaussseidel]");
  AddPar("Poisson_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("Poisson_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("Poisson_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");

  AddPar("Poisson_grid", "SphericalDF",
         "what grid we use [SphericalDF, AnsorgNS, 4ABphi_2xyz]");
  AddPar("Poisson_regularization", "none",
         "options for 4ABphi_2xyz and AnsorgNS " 
         "[regularity_on_axis,regularity_on_axis_at_center]");
	     	   	   	 
  return 0;
}
