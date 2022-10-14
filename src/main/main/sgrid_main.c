/* sgrid_main.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "main.h"



int sgrid_main(void) 
{
  printf("Adding main\n");

  /* functions */

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");
  AddPar("dt"  , "0.125", "time step dt");
  
  AddPar("iterations", "0", "number of grid iterations");
  AddPar("finaltime", "0", "iterate until grid reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  AddPar("errorexit", "exit", "how we exit in case of error [exit,abort]");
  AddPar("verbose", "yes", " [yes,no]");
  return 0;
}
