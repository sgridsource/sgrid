/* coordtrans_CubedSphere.c */
/* Wolfgang Tichy, Nov 2017 */
/* transformations for cubed spheres which really are sphered cubes */

#include "sgrid.h"
#include "Coordinates.h"

/*
We only show domains 0-3. Domains 4 and 5 are in z-dir and are not shown here:

    1 excised cube    1 excised sphere   1 excised cube       spheroidal
    on the inside     on the inside      on the inside        inside and 
    cubical outside   cubical outside    spheroidal outside   outside
y                                              ____
 ^                                           _/    \_ 
 |    ________          ________            /        \            __
 |   |\  3   /|        |\  3   /|          / \  3   / \          /3 \
 |   | \ __ / |        | \ __ / |         /   \ __ /   \        /\  /\
 |   |0 |  | 1|        |0 /  \ 1|        |  0  |  |  1  |      |0 () 1|
 |   |  |__|  |        |  \__/  |        |     |__|     |       \/2 \/
 |   | / 2  \ |        | / 2  \ |         \   /    \   /         \__/
 |   |/______\|        |/______\|          \ /  2   \ /
 |                                          \_      _/  domain size
 |                                            \____/    is arbitrary
-|-----------------------------------------------------------------------> x
  PyramidFrustum    innerCubedSphere      outerCubedSphere     CubedShell
*/

/* Type of cubed sphere or rather sphered cube coord transform 
   PyramidFrustum:    both inner & outer surfaces are flat
   innerCubedSphere:  inner surface is curved, but outer surface is flat
   outerCubedSphere:  outer surface is curved, but inner surface is flat
   CubedShell:        both inner & outer surfaces are curved */


/* compute a sphered cube coord trafo. The type is enumerated above. */
/* Note: lam \in [0,1], (A,B) \in [-1,1] */
void xyz_of_lamAB_CubSph(tBox *box, int ind, double lam, double A, double B,
                         double *x, double *y, double *z)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int type = box->CI->type;
  int dir, p;
  double pm, dx,dy,dz, xc,yc,zc, a0,a1, sigma0,sigma1;

  /* get direction, pm, and center */
  dir = domain/2 + 1;
  p  = 2*(domain%2) - 1; /* p  = +/- 1 */
  pm = p;                /* pm = +/- 1.0 */
  xc = box->CI->xc[1];
  yc = box->CI->xc[2];
  zc = box->CI->xc[3];

  /* check type of trafo */
  if(type==PyramidFrustum)
  {
    /* this gives a pyramid frustum */
    a0 = pm * box->CI->s[0];
    a1 = pm * box->CI->s[1];
  }
  else if(type==innerCubedSphere)
  {
    /* this gives a cubed sphere piece where the inner surface is curved
       and the outer surface is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    a0 = pm * sigma0/sqrt(1.0 + A*A + B*B);
    a1 = pm * box->CI->s[1];
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a1 = pm * sigma1/sqrt(1.0 + A*A + B*B);
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a0 = pm * sigma0/sqrt(1.0 + A*A + B*B);
    a1 = pm * sigma1/sqrt(1.0 + A*A + B*B);
  }
  else errorexit("unknown type");

  /* compute coord trafo for each domain */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc) */
    dx = (a1-a0)*lam + a0;
    *x = xc + dx;
    *y = yc + A*dx;
    *z = zc + B*dx;
  }
  else if(dir==2)
  {
    dy = (a1-a0)*lam + a0;
    *y = yc + dy;
    *x = xc + A*dy;
    *z = zc + B*dy;
  }
  else /* dir==3 */
  {
    dz = (a1-a0)*lam + a0;
    *z = zc + dz;
    *y = yc + A*dz;
    *x = xc + B*dz;
  }
}


void dlamAB_dxyz_SphCube(tBox *box, int ind, int domain,
                         double lam, double A, double B, 
                         double *x, double *y, double *z,
                         double *dlamdx, double *dlamdy, double *dlamdz,
                         double *dAdx,   double *dAdy,   double *dAdz,
                         double *dBdx,   double *dBdy,   double *dBdz)
{
}



#define CALL_xyz_of_lamAB_CubSph \
  tBox *box = (tBox *) aux;\
  double x,y,z;\
  xyz_of_lamAB_CubSph(box, ind, lam,A,B, &x,&y,&z)
#define RET_x CALL_xyz_of_lamAB_CubSph; return x
#define RET_y CALL_xyz_of_lamAB_CubSph; return y
#define RET_z CALL_xyz_of_lamAB_CubSph; return z

/* Coord. trafos of for Cubed Spheres */
double x_of_CubedSphere(void *aux, int ind, double lam, double A, double B)
{ RET_x; }
double y_of_CubedSphere(void *aux, int ind, double lam, double A, double B)
{ RET_y; }
double z_of_CubedSphere(void *aux, int ind, double lam, double A, double B)
{ RET_z; }
