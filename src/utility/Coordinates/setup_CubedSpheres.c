/* setup_CubedSpheres.c */
/* Wolfgang Tichy, Nov 2017 */
/* create various arrangements of cubed spheres */

#include "sgrid.h"
#include "Coordinates.h"




/* init box->CI->iSurf[si] and its derivs from values in box->CI->iFS[si] */
void init_CubedSphere_from_CI_iFS(tBox *box, int si)
{
  int iFS = box->CI->iFS[si];

  /* do nothing if box->CI->iFS[si] is not set */
  if(iFS<1) return;
  
  if(box->CI->useF)
    FSurf_CubSph_init_from_CI_iFS(box, si);
  else /* copy from box->CI->iFS[si] */
  {
    int isigma    = box->CI->iSurf[si];
    int isigma_dA = box->CI->idSurfdX[si][2];
    int isigma_dB = box->CI->idSurfdX[si][3];

    copy_CubedSphere_sigma01_inplane(box, si, iFS, isigma);
    /* and set derivs */
    compute_CubedSphere_dsigma01(box, isigma, isigma_dA, isigma_dB);
  }
}

/* copy a sigma */
void copy_CubedSphere_sigma01_inplane(tBox *box, int si,
                                      int isigma_src, int isigma_dest)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int p = (n1-1) * (si==1);   /* set plane index p */
  int i,j,k;
  double *s = box->v[isigma_src];
  double *d = box->v[isigma_dest];

  forplane1(i,j,k, n1,n2,n3, p)
  {
    int ijk = Index(i,j,k);
    d[ijk] = s[ijk];
  }
}

/* set var sigma01 to a const in first or last plane of dir1,
   si selects plane: si=0: first, si=1: last */
void set_const_CubedSphere_sigma01_inplane(tBox *box, int isigma, int si,
                                           double sigma)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int p = (n1-1) * (si==1);   /* set plane index p */
  int i,j,k;

  forplane1(i,j,k, n1,n2,n3, p)
  {
    int ijk = Index(i,j,k);
    box->v[isigma][ijk] = sigma;
  }
}

/* compute var dsigma01_dA/B everywhere.
   But this could also be done in first or last plane of dir1,
   where si selects plane: si=0: first, si=1: last */
void compute_CubedSphere_dsigma01(tBox *box, int isigma,
                                  int isigma_dA, int isigma_dB)
{
  double *sigma    = box->v[isigma];
  double *dsigmadA = box->v[isigma_dA];
  double *dsigmadB = box->v[isigma_dB];
  /* 
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int p = (n1-1) * (si==1); // set plane index p
  int i,j,k;
  for now we ignore si and compute the deriv everywhere */
  spec_Deriv1(box, 2, sigma, dsigmadA);
  spec_Deriv1(box, 3, sigma, dsigmadB);
}

/* this deforms sigma for testing purposes */
void deform_CubedSphere_sigma01(tBox *box, int isigma, double cA, double cB)
{
  int iA = Ind("Y");
  int iB = Ind("Z");
  int i;
  double s = (box->CI->dom - 2.5)/2.5;
  forallpoints(box,i)
  {
    double A = box->v[iA][i];
    double B = box->v[iB][i];
    double fac = (1.0 - cA*(A-1)*(A+1)*s) * (1.0 - cB*(B-1)*(B+1)*s);
    box->v[isigma][i] *= fac;
  }
}

/* find Amax,Amin, Bmax,Bmin in a domain using distances from center */
void set_AB_min_max_from_Din(int dom, double *Din,
                             double *Amin, double *Amax,
                             double *Bmin, double *Bmax)
{
  switch(dom)
  {
    case 0:
      *Amin = -Din[3]/Din[dom];
      *Amax =  Din[2]/Din[dom];
      *Bmin = -Din[5]/Din[dom];
      *Bmax =  Din[4]/Din[dom];
      break;

    case 1:
      *Amin = -Din[2]/Din[dom];
      *Amax =  Din[3]/Din[dom];
      *Bmin = -Din[4]/Din[dom];
      *Bmax =  Din[5]/Din[dom];
      break;

    case 2:
      *Amin = -Din[1]/Din[dom];
      *Amax =  Din[0]/Din[dom];
      *Bmin = -Din[5]/Din[dom];
      *Bmax =  Din[4]/Din[dom];
      break;

    case 3:
      *Amin = -Din[0]/Din[dom];
      *Amax =  Din[1]/Din[dom];
      *Bmin = -Din[4]/Din[dom];
      *Bmax =  Din[5]/Din[dom];
      break;

    case 4:
      *Amin = -Din[3]/Din[dom];
      *Amax =  Din[2]/Din[dom];
      *Bmin = -Din[1]/Din[dom];
      *Bmax =  Din[0]/Din[dom];
      break;

    case 5:
      *Amin = -Din[2]/Din[dom];
      *Amax =  Din[3]/Din[dom];
      *Bmin = -Din[0]/Din[dom];
      *Bmax =  Din[1]/Din[dom];
      break;

   default:
      errorexit("set_AB_min_max_from_D: unknown dom");
  }
}

/* convert 6 boxes starting with b0 to some kind of cubed spheres */
/* call this after all boxes exist already, so at POST_GRID */
/* type can be "PyramidFrustum", "innerCubedSphere", "outerCubedSphere",
   "CubedShell"
   xc[1..3] = (x,y,z) of coord center for the 6 cubed spheres
   Din[0...5] inner distance from center for cubed sph. domain 0-5
   Dout[0...5] outer distance from center for cubed sph. domain 0-5 */
/* It returns the index of the box right after the last converted box */
int convert_6boxes_to_CubedSphere(tGrid *grid, int b0, int type, int stretch,
                                  double *xc, double *Din, double *Dout)
{
  int isigma    = Ind("Coordinates_CubedSphere_sigma01");
  int isigma_dA = Ind("Coordinates_CubedSphere_dsigma01_dA");
  int isigma_dB = Ind("Coordinates_CubedSphere_dsigma01_dB");
  int isigdef   = Ind("Coordinates_CubedSphere_sigma01_def");
  int isigma_co = Ind("Coordinates_CubedSphere_sigma01_co");
  int i;
  char par[1000];
  char val[1000];
  char name[1000];

  if(b0+6 > grid->nboxes)
    printf("convert_6boxes_to_CubedSphere: nboxes is too small\n");

  switch(type)
  {
    case PyramidFrustum:
    case innerCubedSphere:
    case outerCubedSphere:
    case CubedShell:
      snprintf(name, 999, "%s", "CubedSphere");  break;

    default:
      errorexit("convert_6boxes_to_CubedSphere: unknown type");
  }

  /* do we stretch a CubedShell? */
  if(stretch)
    snprintf(name, 999, "%s", "stretchedCubedSphere");

  /* set box pars */
  for(i=0; i<6; i++)
  {
    tBox *box;
    int d;
    double Amin,Amax, Bmin,Bmax;

    /* break i-loop if we do not have enough boxes */
    if(b0+i >= grid->nboxes) break;
    box = grid->box[b0+i];

    /* set coord type */
    snprintf(par, 999, "box%d_Coordinates", b0+i);
    snprintf(val, 999, "%s", name);
    Sets(par, val);
    /* set basis and min/max in each direction */
    set_AB_min_max_from_Din(i, Din, &Amin,&Amax, &Bmin,&Bmax);
    for(d=1; d<=3; d++)
    {
      snprintf(par, 999, "box%d_basis%d", b0+i, d);
      Sets(par, "ChebExtrema");

      snprintf(par, 999, "box%d_min%d", b0+i, d);
      if(d==1)      Sets(par, "0");
      else if(d==2) Setd(par, Amin);
      else          Setd(par, Bmin);

      snprintf(par, 999, "box%d_max%d", b0+i, d);
      if(d==1)      Sets(par, "1");
      else if(d==2) Setd(par, Amax);
      else          Setd(par, Bmax);
    }

    /* set center */
    for(d=1; d<=3; d++)  box->CI->xc[d] = xc[d];

    /* set inner and outer distance from center */
    box->CI->s[0] = Din[i];
    box->CI->s[1] = Dout[i];

    /* set domain index and type */
    box->CI->dom = i;
    box->CI->type= type;

    /* set sigma vars and iSurf, idSurfdX for them */
    if( (box->v[isigdef]) &&  (stretch==0) )
    {
      switch(type)
      {
        case innerCubedSphere:
          /* compute sigma on first plane in dir1 from box->CI->s[0] */
          set_const_CubedSphere_sigma01_inplane(box, isigdef,0, box->CI->s[0]);
          /* now set coord. info structure */
          box->CI->iFS[0] = isigdef;
          box->CI->iSurf[0] = isigma;
          box->CI->idSurfdX[0][2] = isigma_dA;
          box->CI->idSurfdX[0][3] = isigma_dB;
          break;
      
        case outerCubedSphere:
          /* compute sigma on last plane in dir1 from box->CI->s[1] */
          set_const_CubedSphere_sigma01_inplane(box, isigdef,1, box->CI->s[1]);
          /* now set coord. info structure */
          box->CI->iFS[1] = isigdef;
          box->CI->iSurf[1] = isigma;
          box->CI->idSurfdX[1][2] = isigma_dA;
          box->CI->idSurfdX[1][3] = isigma_dB;
          break;
      
        case CubedShell:
          /* compute sigma on first plane in dir1 from box->CI->s[0] */
          set_const_CubedSphere_sigma01_inplane(box, isigdef,0, box->CI->s[0]);
          /* compute sigma on last plane in dir1 from box->CI->s[1] */
          set_const_CubedSphere_sigma01_inplane(box, isigdef,1, box->CI->s[1]);
          /* now set coord. info structure */
          box->CI->iFS[0] = isigdef;
          box->CI->iSurf[0] = isigma;
          box->CI->idSurfdX[0][2] = isigma_dA;
          box->CI->idSurfdX[0][3] = isigma_dB;
          box->CI->iFS[1] = isigdef;
          box->CI->iSurf[1] = isigma;
          box->CI->idSurfdX[1][2] = isigma_dA;
          box->CI->idSurfdX[1][3] = isigma_dB;
          break;

        default:
          errorexit("convert_6boxes_to_CubedSphere: not sure what to do...");
      }
      /* compute sigma derivs */
      if(box->v[isigdef]!=NULL)
      {
        // /* deform sigma for testing */
        // if(1) deform_CubedSphere_sigma01(box, isigdef, 0.2, -0.1);
        /* enable isigma and its derivs, and set them  */
        enablevar_inbox(box, isigma);
        enablevar_inbox(box, isigma_dA);
        enablevar_inbox(box, isigma_dB);
        init_CubedSphere_from_CI_iFS(box, 0);
        init_CubedSphere_from_CI_iFS(box, 1);
      }
    }

    /* erase all bface info in this box */
    free_all_bfaces(box);
  }

  /* convert boxes to what the pars now say */
  set_BoxStructures_fromPars(grid, 0);

  return b0+i; /* return box index right after last added box */
}


/* disable vars in box->CI->iSurf box->CI->idSurfdX and set all to zero */
void disable_and_reset_CI_iSurf_vars(tBox *box)
{
  int i, j, vi;
  for(i=0; i<6; i++)
  {
    vi = box->CI->iFS[i];
    if(vi>0) disablevar_inbox(box, vi);
    box->CI->iFS[i] = 0;

    vi = box->CI->iSurf[i];
    if(vi>0) disablevar_inbox(box, vi);
    box->CI->iSurf[i] = 0;

    for(j=0; j<4; j++)
    {
      vi = box->CI->idSurfdX[i][j];
      if(vi>0) disablevar_inbox(box, vi);
      box->CI->idSurfdX[i][j] = 0;
    }
  }
}

/* disable Coordinates_CubedSphere_sigma01 and its derivs in a box */
void disable_Coordinates_CubedSphere_sigma01(tBox *box)
{
  int isigma     = Ind("Coordinates_CubedSphere_sigma01");
  int isigma_dA  = Ind("Coordinates_CubedSphere_dsigma01_dA");
  int isigma_dB  = Ind("Coordinates_CubedSphere_dsigma01_dB");
  int isigma_def = Ind("Coordinates_CubedSphere_sigma01_def");
  int isigma_co  = Ind("Coordinates_CubedSphere_sigma01_co");
  int i, vi;
  for(i=0; i<6; i++)
  {
    vi = box->CI->iSurf[i];
    if(vi == isigma) disable_and_reset_CI_iSurf_vars(box);
  }
  disablevar_inbox(box, isigma);
  disablevar_inbox(box, isigma_dA);
  disablevar_inbox(box, isigma_dB);
  disablevar_inbox(box, isigma_def);
  disablevar_inbox(box, isigma_co);
}

/* convert 1 box to a cube centered at xc[i],
   returns the index of the box right after the last converted box */
int convert_1box_to_cube(tGrid *grid, int b0, double *xc, double dout)
{
  char par[1000];
  char val[1000];
  tBox *box = grid->box[b0];
  int d;

  if(b0+1 > grid->nboxes)
  {
    printf("convert_1box_to_cube: nboxes is too small\n");
    return b0;
  }

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
  int i;
  double Din[6];
  double Dmid[6];
  double Dout[6];

  /* set distances from center to din, dmid, dout for all domains */
  for(i=0; i<6; i++)
  {
    Din[i] = din;
    Dmid[i]= dmid;
    Dout[i]= dout;
  }
  /* convert the 12 boxes */
  bl = convert_6boxes_to_CubedSphere(grid, bl, outerCubedSphere,0,
                                     xc, Din,Dmid);
  bl = convert_6boxes_to_CubedSphere(grid, bl, innerCubedSphere,0,
                                     xc, Dmid,Dout);
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


/* surround two big touching cubes with cubed spheres that have
   A,B that are in an extended range
         _______
      __/       \__ 
     /             \   
    /__     3     __\  
   /   -- _____ --   \     e.g. dom0/1 have A = [-0.8, 0.8] 
  |   0  |  |  |   1  |         dom2/3 have A = [-1.25,1.25]
  |      |__|__|      |
   \ __--       --__ /      r0 is radius of outer sphere
    \       2       /       2*dc is sidelength of one cube 
     \__         __/
        \_______/ 
*/
int sphere_around_two_full_cubes_touching_at_x0(tGrid *grid, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0)
{
  int bl=b0;
  double xc[4], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes */
  bl = two_full_cubes_touching_at_x0(grid, b0, dc, din1,dmid1, din2,dmid2);

  /* set distances to make 6 more cubed spheres around these 2 full cubes */
  for(i=0; i<6; i++)
  {
    if(i<2) Din[i] = 2.0*dc;
    else    Din[i] = dc;
    Dout[i] = r0;
  }  
  xc[1] = xc[2] = xc[3] = 0.0;
  bl = convert_6boxes_to_CubedSphere(grid, bl, outerCubedSphere,0,
                                     xc, Din,Dout);
  return bl;
}


/* put 6 stretchedCubedShell's around the sphere from
   sphere_around_two_full_cubes_touching_at_x0
                   ___________
             _____/           \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/       \__  -/   ...
    |           \/             \/   
   /            /__           __\  
  |            /   -- _____ --   \    r0 is radius of inner sphere (sphere0)
  |           |      |  |  |      |   r1 is radius of outer sphere (sphere1)
  |           |      |__|__|      |
  |            \ __--       --__ / 
   \            \               /  
    |          _/\__         __/
     |       _/     \_______/ 
      \    _/  
       -_ /
         \      ...
*/
int two_spheres_around_two_full_cubes(tGrid *grid, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1)
{
  int bl=b0;
  double xc[4], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes and sphere0 around them */
  bl = sphere_around_two_full_cubes_touching_at_x0(grid, b0, dc,
                                                   din1,dmid1, din2,dmid2, r0);
  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }  
  xc[1] = xc[2] = xc[3] = 0.0;
  bl = convert_6boxes_to_CubedSphere(grid, bl, CubedShell,1, xc, Din,Dout);
  return bl;
}
