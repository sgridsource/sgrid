/* sgrid_testwave.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "testwave.h"


int sgrid_testwave() 
{
  printf("Adding testwave\n");

  /* functions */
  AddFun(INITIALDATA, testwave_startup, "initialize test");
  AddFun(ANALYZE, testwave_analyze, "compute test results");
  AddFun(POST_EVOLVE, testwave_filter, "filter u and k");

  /* variables */
  AddVar("testwave_u", "",  "wave function");
  AddVar("testwave_k", "",  "time deriv of wave function");
  AddVar("testwave_u11", "", "derivs of wave function");
  AddVar("testwave_u22", "", "derivs of wave function");
  AddVar("testwave_u33", "", "derivs of wave function");
  AddVar("testwave_dum", "", "dummy var");
  
  /* parameters */
  AddPar("testwave_persist", "yes", "whether additional memory persists");
	     	   	   	 
  return 0;
}
