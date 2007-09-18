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
  AddVar("x", "", "cartesian x coordinate");
  AddVar("y", "", "cartesian y coordinate");
  AddVar("z", "", "cartesian z coordinate");
  AddVar("dXd", "i", "coordinate derivative dX/dx^i");
  AddVar("dYd", "i", "coordinate derivative dY/dx^i");
  AddVar("dZd", "i", "coordinate derivative dZ/dx^i");
  AddVar("ddXdd", "(ij)", "2nd coordinate derivative d^2 X/dx^i dx^j");
  AddVar("ddYdd", "(ij)", "2nd coordinate derivative d^2 Y/dx^i dx^j");
  AddVar("ddZdd", "(ij)", "2nd coordinate derivative d^2 Z/dx^i dx^j");

  /* parameters */
  for(b=0; b<nboxes; b++)
  {
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    AddPar(str, "Cartesian", 
           "coordinates used in box [Cartesian, Polar, ...]");
  }

  AddPar("CoordinateTransforms_stored", "yes",
         "whether we store Coordinate Transforms in dXdx,... ddXddxx,...");
  AddPar("Coordinates_newtTOLF", "1e-10", "newton tolerence");
  AddPar("Coordinates_newtMAXITS", "100000", "max. newton iterations");
  AddPar("compactSphericalDF_r0", "-1", "radius r at xi=0");
  AddPar("tan_stretch_s", "0", "how much we stretch [0,Xmax]");
  AddPar("Coordinates_Spherical3_c",  "5",
         "constant c in Spherical3 coord trafo");

  return 0;
}
