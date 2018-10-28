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

    /* do we need sigma var */
    if(Getv("TestDerivs_grid","DeformCubSph"))
    {
      int b;
      for(b=1; b<min2(13, grid->nboxes); b++)
      {
        tBox *box = grid->box[b];
        enablevar_inbox(box, Ind("Coordinates_CubedSphere_sigma01_def"));
      }
    }

    /* test cubed spheres */
    switch(grid->nboxes)
    {
      case 13:
        arrange_1box12CubSph_into_full_cube(grid, 0, xc, 0.2,0.4,0.6);
        break;
      case 26:
        two_full_cubes_touching_at_x0(grid, 0, xc[1], 0.2,0.4, 0.25,0.5);
        break;
      case 32:
        sphere_around_two_full_cubes_touching_at_x0(grid, 0, xc[1],
                                                    0.2,0.4, 0.25,0.5, 4.0);
        break;
      case 38:
        two_spheres_around_two_full_cubes(grid, 0, xc[1],
                                          0.2,0.4, 0.25,0.5, 4.0,40.0);
        break;
      default:
        errorexit("nboxes should be 13, 26, 32, or 38");
    }

    /* deform sigma01_def, and switch on sigma func and its derivs */
    if(Getv("TestDerivs_grid","DeformCubSph"))
    {
      int b;
      for(b=1; b<min2(13, grid->nboxes); b++)
      {
        tBox *box = grid->box[b];
        int n1=box->n1;
        int n2=box->n2;
        int n3=box->n3;
        int i,j,k, p;
        double *Y = box->v[Ind("Y")];
        double *Z = box->v[Ind("Z")];
        double *sig;
        int si;

        /* select surface */
        if(b<7) si=1;
        else    si=0;
        sig = box->v[box->CI->iFS[si]];
        p = si*(n1-1);

        /* deform sigma01_def */
        forplane1(i,j,k, n1,n2,n3, p)
        {
          double th,ph;
          int ijk = Index(i,j,k);

          ThetaPhi_of_AB_CubSph(box, Y[ijk],Z[ijk], &th,&ph);

          sig[ijk] *= (1. + 0.1*sqrt(15./(8*PI))*sin(th)*cos(th)*cos(ph)
                          + 0.1*sqrt(15./(8*PI))*sin(th)*sin(th)*sin(2*ph));
        }
        /* set flags to sigma func */
        if(Getv("TestDerivs_grid","useF")) box->CI->useF = 1;
      }
    } /* end DeformCubSph */
  }

  parameterio_write_current_pars(grid);
  return 0;
}
