/* sgrid_Coordinates.c */
/* Wolfgang Tichy 2003 */

#include "sgrid.h"
#include "Coordinates.h"


int sgrid_Coordinates(void) 
{
  int nboxes = Geti("nboxes");
  int b;
  
  printf("Adding Coordinates\n");

  /* functions */
  AddFun(PRE_INITIALDATA, init_CoordTransform_And_Derivs, 
         "initialize coords and coord transforms");

  /* variables */
  AddConstantVar("x", "", "cartesian x coordinate");
  AddConstantVar("y", "", "cartesian y coordinate");
  AddConstantVar("z", "", "cartesian z coordinate");

  /* parameters */
  for(b=0; b<nboxes; b++)
  {
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    AddPar(str, "Cartesian", 
           "coordinates used in box [Cartesian, Polar]");
  }

  AddPar("compactSphericalDF_r0", "-1", "radius r at xi=0");

  return 0;
}
