/* sgrid_GridIterators.c */
/* Wolfgang Tichy 8/2008 */

#include "sgrid.h"
#include "GridIterators.h"



int sgrid_GridIterators() 
{
  printf("Adding GridIterators\n");

  /* global parameter */
  AddPar("GridIterators_verbose", "yes", "talk about it");
  AddPar("GridIterators_GMRES_restart", "max", 
         "Restart parameter for GMRES[some #,max");
  AddPar("GridIterators_setABStozero_below", "0", 
         "some numbers will be set equal to 0 if their magnitude "
         "is below this value [any positive number]");

  /* check whether there is more to do */
  if (!Getv("physics", "GridIterators")) return 0;

  /* functions */

  /* variables */
  //AddVar("u", "", "scalar field for Poisson equation");

/* temporary variables, should be inside if statement to avoid
     cluttering the variable name space 
  AddVar("GridIterators_p",  "", "temporary variable");
  AddVar("GridIterators_ph", "", "temporary variable");
  AddVar("GridIterators_rt", "", "temporary variable");
  AddVar("GridIterators_s",  "", "temporary variable");
  AddVar("GridIterators_sh", "", "temporary variable");
  AddVar("GridIterators_t",  "", "temporary variable");
  AddVar("GridIterators_v",  "", "temporary variable");
*/
  return 0;
}
