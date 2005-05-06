/* sgrid_ConvTest.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "ConvTest.h"


int sgrid_ConvTest() 
{
  printf("Adding ConvTest\n");

  /* functions */
  AddFun(INITIALDATA, ConvTest_startup, "initialize test");
  AddFun(ANALYZE, ConvTest_analyze, "compute test results");
  AddFun(POST_EVOLVE, ConvTest_filter, "filter u and k");

  /* variables */
  AddVar("ConvTest_u", "",  "wave function");
  AddVar("ConvTest_u1", "", "derivs of wave function");
  AddVar("ConvTest_u2", "", "derivs of wave function");
  AddVar("ConvTest_u3", "", "derivs of wave function");
  AddVar("ConvTest_err", "", "relative error");
  
  /* parameters */
  AddPar("ConvTest_persist", "yes", "whether additional memory persists");
	     	   	   	 
  return 0;
}
