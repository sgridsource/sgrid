/* sgrid_ParManipulator.c */
/* Bernd Bruegmann 01/00 */

#include "sgrid.h"
#include "ParManipulator.h"



int sgrid_ParManipulator() 
{
  printf("Adding ParManipulator\n");

  /* functions */
  AddFun(PRE_GRID, ParManipulator_setPars, "set some pars automatically");
         
  /* parameters */
  /* add some pars that can be used in each box, e.g.
     set "box2_n1 = n1" to use n1 */
  AddPar("n1", "1", "number of points in X-direction");
  AddPar("n2", "1", "number of points in Y-direction");
  AddPar("n3", "1", "number of points in Z-direction");
  AddPar("min1", "-1", "lower boundary in X-direction");
  AddPar("max1", "+1", "upper boundary in X-direction");
  AddPar("min2", "-1", "lower boundary in Y-direction");
  AddPar("max2", "+1", "upper boundary in Y-direction");
  AddPar("min3", "-1", "lower boundary in Z-direction");
  AddPar("max3", "+1", "upper boundary in Z-direction");
  AddPar("basis1", "ChebExtrema", "basis functions in X-direction");
  AddPar("basis2", "ChebExtrema", "basis functions in Y-direction");
  AddPar("basis3", "ChebExtrema", "basis functions in Z-direction");

  return 0;
}
