/* sgrid_Poisson3.c */
/* Wolfgang Tichy, Oct. 2017 */

#include "sgrid.h"
#include "Poisson3.h"


int sgrid_Poisson3() 
{
  if (!Getv("physics", "Poisson3")) return 0;
  printf("Adding Poisson3\n");

  /* functions */
  AddFun(PRE_COORDINATES, Poisson3_initboxes, "initialize boxes we use");
  AddFun(PRE_INITIALDATA, Poisson3_startup, "initialize Poisson3");
  AddFun(INITIALDATA, Poisson3_solve, "solve Poisson3 Eq.");
  AddFun(ANALYZE, Poisson3_analyze, "compute error");

  /* variables */
  AddVar("Poisson3_Psi",     "",
         "function that satisfies Poisson3 Eq: Laplace Psi = rh1");
  AddVar("Poisson3_Psi",     "i",    "1st deriv of Psi");
  AddVar("Poisson3_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("Poisson3_rh1",     "",     "Source on RHS1");
  AddVar("Poisson3_Err_Psi", "",     "Error in Psi");

  AddVar("Poisson3_Chi",     "",
         "another function that satisfies Poisson3: Laplace Chi = rh2");
  AddVar("Poisson3_Chi",     "i",    "1st deriv of Chi");
  AddVar("Poisson3_Chi",     "(ij)", "2nd deriv of Chi");
  AddVar("Poisson3_rh2",     "",     "Source on RHS2");
  AddVar("Poisson3_Err_Chi", "",     "Error in Chi");

  AddVar("Poisson3_temp1", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson3_temp2", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson3_temp3", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson3_temp4", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson3_temp5", "", 	"temporary variable(e.g. to store derivs)");
  AddVar("Poisson3_temp6", "",   "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("Poisson3_itmax", "10", "maximal number of Newton iterations");
  AddPar("Poisson3_tol",   "1e-6","tolerance for Newton step or multigrid");
  AddPar("Poisson3_linSolver", "bicgstab", "linear solver used "
         "[bicgstab,UMFPACK,templates_GMRES_with_Jacobi_precon]");
  AddPar("Poisson3_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("Poisson3_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("Poisson3_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("Poisson3_linSolver_Precon", "I", 
         "Preconditioner used [I,fd_UMFPACK]");
  AddPar("Poisson3_grid", "", "what grid we use [2starcubes]");
	     	   	   	 
  return 0;
}
