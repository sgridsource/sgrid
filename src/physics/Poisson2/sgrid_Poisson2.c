/* sgrid_Poisson2.c */
/* Wolfgang Tichy, June 2016 */

#include "sgrid.h"
#include "Poisson2.h"


int sgrid_Poisson2() 
{
  if (!Getv("physics", "Poisson2")) return 0;
  printf("Adding Poisson2\n");

  /* functions */
  AddFun(PRE_INITIALDATA, Poisson2_startup, "initialize Poisson2");
  AddFun(INITIALDATA, Poisson2_solve, "solve Poisson2 Eq.");
  AddFun(ANALYZE, Poisson2_analyze, "compute error");

  /* variables */
  AddVar("Poisson2_Psi",     "",
         "function that satisfies Poisson2 Eq: Laplace Psi = rh1");
  AddVar("Poisson2_Psi",     "i",    "1st deriv of Psi");
  AddVar("Poisson2_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("Poisson2_rh1",     "",     "Source on RHS1");
  AddVar("Poisson2_Err_Psi", "",     "Error in Psi");

  AddVar("Poisson2_Chi",     "",
         "another function that satisfies Poisson2: Laplace Chi = rh2");
  AddVar("Poisson2_Chi",     "i",    "1st deriv of Chi");
  AddVar("Poisson2_Chi",     "(ij)", "2nd deriv of Chi");
  AddVar("Poisson2_rh2",     "",     "Source on RHS2");
  AddVar("Poisson2_Err_Chi", "",     "Error in Chi");

  AddVar("Poisson2_temp1", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson2_temp2", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson2_temp3", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson2_temp4", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson2_temp5", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson2_temp6", "",   "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("Poisson2_itmax", "10", "maximal number of Newton iterations");
  AddPar("Poisson2_tol",   "1e-6","tolerance for Newton step or multigrid");
  AddPar("Poisson2_linSolver", "bicgstab", "linear solver used "
         "[bicgstab,UMFPACK,templates_GMRES_with_Jacobi_precon]");
  AddPar("Poisson2_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("Poisson2_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("Poisson2_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("Poisson2_linSolver_Precon", "I", 
         "Preconditioner used [I,fd_UMFPACK]");
         
  AddPar("Poisson2_grid", "SphericalDF",
         "what grid we use [SphericalDF]");
	     	   	   	 
  return 0;
}
