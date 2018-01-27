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
    switch(grid->nboxes)
    {
      case 13:
        arrange_1box12CubSph_into_full_cube(grid, 0, xc, 0.2,0.4,0.6);
        break;
      case 26:
        two_full_cubes_touching_at_x0(grid, 0, xc[1], 0.2,0.4, 0.3,0.5);
        break;
      case 32:
        sphere_around_two_full_cubes_touching_at_x0(grid, 0, xc[1],
                                                    0.2,0.4, 0.3,0.5, 4.0);
        break;
      case 38:
        two_spheres_around_two_full_cubes(grid, 0, xc[1],
                                          0.2,0.4, 0.3,0.5, 4.0,40.0);
        break;
      default:
        errorexit("nboxes should be 13, 26, 32, or 38");
    }
  }
  parameterio_write_current_pars(grid);
  return 0;
}
