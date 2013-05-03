/* sgrid_NumberChecker.c */
/* Wolfgang Tichy 9/2008 */

#include "sgrid.h"
#include "NumberChecker.h"



int sgrid_NumberChecker(void)
{
  printf("Adding NumberChecker\n");

  /* pars to check for NANs */
  AddPar("NumberChecker_exitifNAN", "no", 
	 "exit if NAN is found [varname,no] (not a list)");
  AddPar("NumberChecker_INF", "1e300", 
	 "consider a number INF/NAN if it's above this [any large number]");

  /* add functions only if NumberChecker is on */
  if (!Getv("physics", "NumberChecker")) return 0;
  /* functions */
  // AddFun(ANALYZE, computeADMconstraints, "compute ADM constraints");

  /* 3+1 field variables */
  // AddVar("g",      "(ij)", "metric");

  if (!Getv("NumberChecker_exitifNAN", "no")) 
    AddFun(POST_EVOLVE, NumberChecker_ExitIfNAN, "exit if NAN");

  return 0;
}
