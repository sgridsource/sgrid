/* sgrid_testwave.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "testwave.h"


int sgrid_testwave() 
{
  printf("Adding testwave\n");

  /* functions */
  AddFun(INITIALDATA, testwave_startup, "initialize test");
  // AddFun(ANALYZE, testwave_analyze, "compute test results");
  AddFun(POST_EVOLVE, testwave_filter, "filter u and k");

  /* variables */
  AddVar("testwave_u", "",  "wave function");
  AddVar("testwave_k", "",  "time deriv of wave function");
  AddVar("testwave_u1", "", "derivs of wave function");
  AddVar("testwave_u2", "", "derivs of wave function");
  AddVar("testwave_u3", "", "derivs of wave function");
  AddVar("testwave_u11", "", "2nd derivs of wave function");
  AddVar("testwave_u22", "", "2nd derivs of wave function");
  AddVar("testwave_u33", "", "2nd derivs of wave function");
  AddVar("testwave_dum", "", "dummy var");
  
  /* parameters */
  AddPar("testwave_persist", "yes", "whether additional memory persists");
	     	   	   	 
  return 0;
}
