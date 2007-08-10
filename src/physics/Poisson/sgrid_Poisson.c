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
  AddVar("Poisson_Psi",     "",     "function that satisfies Poisson Eq.");
  AddVar("Poisson_Psi",     "i",    "1st deriv of Psi");
  AddVar("Poisson_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("Poisson_Err_Psi", "",     "Error in Psi");
  AddVar("Poisson_Chi",     "",     "another function that satisfies Poisson");
  AddVar("Poisson_Chi",     "i",    "1st deriv of Chi");
  AddVar("Poisson_Chi",     "(ij)", "2nd deriv of Chi");
  AddVar("Poisson_Err_Chi", "",     "Error in Chi");
   
  /* parameters */
  AddPar("Poisson_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("Poisson_itmax", "10", "maximal number of Newton iterations");
  AddPar("Poisson_tol",   "1e-6","tolerance for Newton step or multigrid");
//  AddPar("Poisson_linSolver", "bicgstab", 
//         "linear solver used [bicgstab, HYPRE, gaussseidel]");
  AddPar("Poisson_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("Poisson_linSolver_tolFac","0.1", 
         "tol for linSolver is tol * linSolver_tolFac");
	     	   	   	 
  return 0;
}
