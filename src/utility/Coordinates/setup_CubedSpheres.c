/* setup_CubedSpheres.c */
/* Wolfgang Tichy, Nov 2017 */
/* create various arrangements of cubed spheres */

#include "sgrid.h"
#include "Coordinates.h"


/* convert 6 boxes starting with b0 to some kind of cubed spheres */
/* call this after all boxes exist already, so at POST_GRID */
/* type can be "PyramidFrustum", "innerCubedSphere", "outerCubedSphere",
   "CubedShell"
   It returns the index of the box right after the last converted box */
int convert_6boxes_to_CubedSphere(tGrid *grid, int b0, int type,
                                  double *xc, double din, double dout)
{
  int i;
  char par[1000];
  char val[1000];
  char name[1000];

  if(b0+6 > grid->nboxes)
    errorexit("convert_6boxes_to_CubedSphere: nboxes is too small");

  switch(type)
  {
    case PyramidFrustum:
      snprintf(name, 999, "%s", "PyramidFrustum");  break;
    
    case innerCubedSphere:
      snprintf(name, 999, "%s", "innerCubedSphere");  break;
  
    case outerCubedSphere:
      snprintf(name, 999, "%s", "outerCubedSphere");  break;
  
    case CubedShell:
      snprintf(name, 999, "%s", "CubedShell");  break;

    default:
      errorexit("convert_6boxes_to_CubedSphere: unkonwn type");
  }

  /* set box pars */
  for(i=0; i<6; i++)
  {
    tBox *box = grid->box[b0+i];
    int d;
    /* set coord type */
    snprintf(par, 999, "box%d_Coordinates", b0+i);
    snprintf(val, 999, "%s", name);
    Sets(par, val);
    /* set basis and min/max in each direction */
    for(d=1; d<=3; d++)
    {
      snprintf(par, 999, "box%d_basis%d", b0+i, d);
      Sets(par, "ChebExtrema");

      snprintf(par, 999, "box%d_min%d", b0+i, d);
      if(d==1) Sets(par, "0");
      else     Sets(par, "-1");

      snprintf(par, 999, "box%d_max%d", b0+i, d);
      Sets(par, "1");
    }

    /* set center */
    for(d=1; d<=3; d++)  box->CI->xc[d] = xc[d];

    /* set inner and outer distance from center */
    box->CI->s[0] = din;
    box->CI->s[1] = dout;

    /* set domain index and type */
    box->CI->dom = i;
    box->CI->type= type;

    /* erase all bface info in this box */
    free_all_bfaces(box);

    /* FIXME: maybe add varlist for inner or outer distance from center??? */
  }

  /* convert boxes to what the pars now say */
  set_BoxStructures_fromPars(grid, 0);

  return b0+6; /* return box index right after last added box */
}

/* convert 1 box to a cube centered at xc[i],
   returns the index of the box right after the last converted box */
int convert_1box_to_cube(tGrid *grid, int b0, double *xc, double dout)
{
  int i;
  char par[1000];
  char val[1000];
  tBox *box = grid->box[b0+i];
  int d;

  if(b0+1 > grid->nboxes)
    errorexit("convert_1box_to_cube: nboxes is too small");

  /* set coord type */
  snprintf(par, 999, "box%d_Coordinates", b0);
  Sets(par, "Cartesian");
  /* set basis and min/max in each direction */
  for(d=1; d<=3; d++)
  {
    snprintf(par, 999, "box%d_basis%d", b0, d);
    Sets(par, "ChebExtrema");

    snprintf(par, 999, "box%d_min%d", b0, d);
    snprintf(val, 999, "%.16g", xc[d]-dout);
    Sets(par, val);

    snprintf(par, 999, "box%d_max%d", b0, d);
    snprintf(val, 999, "%.16g", xc[d]+dout);
    Sets(par, val);
  }

  /* erase all bface info in this box */
  free_all_bfaces(box);

  /* convert boxes to what the pars now say */
  set_BoxStructures_fromPars(grid, 0);

  return b0+1; /* return box index right after last added box */
}


/* make this:
    __________
   |\    9   /|
   | \  __  / |
   |  \/3 \/  |
   |6 /\__/\  |
   | |0|__|1| |
   |  \/2 \/ 7|
   |  /\__/\  |
   | /   8  \ |
   |/________\|

   returns the index of the box right after the last converted box */
int arrange_12CubSph_into_empty_cube(tGrid *grid, int b0, double *xc,
                                     double din, double dmid, double dout)
{
  int bl=b0;
  bl = convert_6boxes_to_CubedSphere(grid, bl, outerCubedSphere,
                                     xc, din,dmid);
  bl = convert_6boxes_to_CubedSphere(grid, bl, innerCubedSphere,
                                     xc, dmid,dout);
  return bl;
}

/* same as arrange_12CubSph_into_empty_cube, but add one cube at the center */
int arrange_1box12CubSph_into_full_cube(tGrid *grid, int b0, double *xc,
                                        double din, double dmid, double dout)
{
  int bl=b0;
  bl = convert_1box_to_cube(grid, bl, xc, din);
  bl = arrange_12CubSph_into_empty_cube(grid, bl, xc, din,dmid,dout);
  return bl;
}


/* make two big touching cubes like this
    __________ __________
   |\        /|\        /|
   | \  __  / | \  __  / |
   |  \/  \/  |  \/  \/  |
   |  /\__/\  |  /\__/\  |
   | | |__| | | | |__| | |
   |  \/  \/  |  \/  \/  |
   |  /\__/\  |  /\__/\  |
   | /      \ | /      \ |
   |/________\|/________\|
*/
int two_full_cubes_touching_at_x0(tGrid *grid, int b0, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2)
{
  int bl=b0;
  double xc[4];

  /* put centers on x-axis */
  xc[2] = xc[3] = 0.0;

  /* full cube 1 is centered at xc[1]=dc and has width 2dc */
  xc[1] = dc;
  bl = arrange_1box12CubSph_into_full_cube(grid, bl, xc, din1,dmid1, dc);

  /* full cube 2 is centered at xc[1]=-dc and has width 2dc */
  xc[1] = -dc;
  bl = arrange_1box12CubSph_into_full_cube(grid, bl, xc, din2,dmid2, dc);

  /* set bface info for the 1 face at x=0 where the two outer cubes
     are now touching */
  //....
  // do all faces have common points??? Assume this next:              
  /* all faces are touching and have the same points */
  //set_touching_bfaces_of_boxes_with_same_facepoints(grid, b0, bl-b0);

  /* do not set bface info yet, because this is likely called before
     init_CoordTransform_And_Derivs */

  return bl;
}