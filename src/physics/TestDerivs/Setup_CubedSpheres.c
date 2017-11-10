/* Setup_CubedSpheres.c */
/* Wolfgang Tichy 2017 */

#include "sgrid.h"
#include "TestDerivs.h"



int Setup_CubedSpheres_initboxes(tGrid *grid)
{
  if(Getv("TestDerivs_grid","CubedSpheres"))
  {
    double xc[4];

    printf("Initializing Cubed Spheres:\n");
    xc[2] = xc[3] = 0.0;
    xc[1] = 1.0;

    /* test cubed spheres */
    arrange_1box12CubSph_into_full_cube(grid, 0, xc, 0.2,0.4,0.6);
  }
  return 0;
}
