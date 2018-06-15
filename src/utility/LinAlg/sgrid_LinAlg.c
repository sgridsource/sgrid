/* sgrid_LinAlg.c */
/* Wolfgang Tich 8/2007 */

#include "sgrid.h"
#include "LinAlg.h"


int sgrid_LinAlg(void) 
{
  /* if (!Getv("physics", "LinAlg")) return; */
  printf("Adding LinAlg\n");

  /* functions */

  /* variables */

  /* parameters */
  AddPar("LinAlg_setABStozero_below", "0", "some numbers will be set to 0 "
         "if their magnitude is below this value [any positive number]");

  return 0;
}
