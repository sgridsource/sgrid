/* sgrid_evolve.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "evolve.h"




int sgrid_evolve() 
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(EVOLVE, evolve, "evolve with method of lines like icn, rk, dlf");

  /* variables */

  /* parameters */
  AddPar("evolution_method",   "icn", "evolution method [icn,rk,euler]");
  AddPar("evolution_method_rk", "rk4", 
	 "Runge-Kutta method [euler,midpoint,rk3,rk3a,rk3b,rk3c,rk4]");
  AddPar("evolve_persist",     "yes", "whether additional memory persists");
  AddPar("evolve_euler_debug", "no",  
	 "obtain rhs in variable after one timestep");
  AddPar("evolve_compute_change", "", 
	 "list of variables for which change is to be computed");
	     	   	   	 
  /* test */
  if (Getv("physics", "evolve_test")) {
    AddFun(INITIALDATA, evolve_test_startup, "initialize test");
    AddFun(ANALYZE, evolve_test_analyze, "compute test results");
    AddVar("evolve_test_u", "", "test function");
    AddVar("evolve_test_error", "", "error function");
    AddPar("evolution_method_order", "0", "expected order of convergence");
  }
  return 0;
}
