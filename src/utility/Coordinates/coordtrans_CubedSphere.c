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
  double pm, rx,ry,rz, xc,yc,zc, a0,a1, sigma0,sigma1;
  double sqrt_1_A2_B2 = sqrt(1.0 + A*A + B*B);

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
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0], func sigma_CubSph below */
    a0 = pm * sigma0/sqrt_1_A2_B2;
    a1 = pm * box->CI->s[1];
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a1 = pm * sigma1/sqrt_1_A2_B2;
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a0 = pm * sigma0/sqrt_1_A2_B2;
    a1 = pm * sigma1/sqrt_1_A2_B2;
  }
  else errorexit("unknown type");

  /* compute coord trafo for each domain */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc) */
    rx = (a1-a0)*lam + a0;
    *x = xc + rx;
    *y = yc + A*rx;
    *z = zc + B*rx;
  }
  else if(dir==2)
  {
    ry = (a1-a0)*lam + a0;
    *y = yc + ry;
    *x = xc + A*ry;
    *z = zc + B*ry;
  }
  else /* dir==3 */
  {
    rz = (a1-a0)*lam + a0;
    *z = zc + rz;
    *y = yc + A*rz;
    *x = xc + B*rz;
  }
}

/* compute inverse sphered cube coord trafo. */
/* Note: lam \in [0,1], (A,B) \in [-1,1] */
void lamAB_of_xyz_CubSph(tBox *box, int ind, double x, double y, double z,
                         double *lam, double *A, double *B)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int type = box->CI->type;
  int dir, p;
  double pm, rx,ry,rz,rc, xc,yc,zc, a0,a1, sigma0,sigma1;
  double oosqrt_1_A2_B2;

  /* get direction, pm, and center */
  dir = domain/2 + 1;
  p  = 2*(domain%2) - 1; /* p  = +/- 1 */
  pm = p;                /* pm = +/- 1.0 */
  xc = box->CI->xc[1];
  yc = box->CI->xc[2];
  zc = box->CI->xc[3];

  /* get d=istance from center, and sqrt_1_A2_B2 */
  rx = x-xc;
  ry = y-yc;
  rz = z-zc;
  rc = sqrt(rx*rx + ry*ry + rz*rz);
  if(dir==1)       oosqrt_1_A2_B2 = fabs(rx)/rc;
  else if(dir==2)  oosqrt_1_A2_B2 = fabs(ry)/rc;
  else             oosqrt_1_A2_B2 = fabs(rz)/rc;

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
    a0 = pm * sigma0*oosqrt_1_A2_B2;
    a1 = pm * box->CI->s[1];
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a1 = pm * sigma1*oosqrt_1_A2_B2;
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a0 = pm * sigma0*oosqrt_1_A2_B2;
    a1 = pm * sigma1*oosqrt_1_A2_B2;
  }
  else errorexit("unknown type");

  /* compute coord trafo for each domain */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc)
       (1 + A^2 + B^2) rx^2 = rx^2 + ry^2 + rz^2
       1/sqrt(1 + A^2 + B^2) = |rx|/sqrt(rx^2 + ry^2 + rz^2) = |rx|/rc */
    *lam = (rx-a0)/(a1-a0);
    *A   = ry/rx;
    *B   = rz/rx;
  }
  else if(dir==2)
  {
    *lam = (ry-a0)/(a1-a0); 
    *A   = rx/ry;
    *B   = rz/ry;
  }
  else /* dir==3 */
  {
    *lam = (rz-a0)/(a1-a0);
    *A   = ry/rz;
    *B   = rx/rz;
  }
}

/* compute derivs of inverse trafo */
void dlamAB_dxyz_CubSph(tBox *box, int ind, double lam, double A, double B, 
                        double *x, double *y, double *z,
                        double *dlamdx, double *dlamdy, double *dlamdz,
                        double *dAdx,   double *dAdy,   double *dAdz,
                        double *dBdx,   double *dBdy,   double *dBdz)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int type = box->CI->type;
  int dir, p;
  double pm, rx,ry,rz, rx2,ry2,rz2, rc2, xc,yc,zc, a0,a1, sigma0,sigma1;
  double sqrt_1_A2_B2 = sqrt(1.0 + A*A + B*B);
  double mx,my,mz, da0_dx,da0_dy,da0_dz, da1_dx,da1_dy,da1_dz;
  double mx0,my0,mz0, mx1,my1,mz1;
  double dsigma_dA_osigma0, dsigma_dB_osigma0;
  double dsigma_dA_osigma1, dsigma_dB_osigma1;

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
    dsigma_dA_osigma0 = 0.0;
    dsigma_dB_osigma0 = 0.0;
    dsigma_dA_osigma1 = 0.0;
    dsigma_dB_osigma1 = 0.0;
    a0 = pm * box->CI->s[0];
    a1 = pm * box->CI->s[1];
    da0_dx = 0.0;
    da0_dy = 0.0;
    da0_dz = 0.0;
    da1_dx = 0.0;
    da1_dy = 0.0;
    da1_dz = 0.0;
  }
  else if(type==innerCubedSphere)
  {
    /* this gives a cubed sphere piece where the inner surface is curved
       and the outer surface is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    dsigma_dA_osigma0 = 0.0;
    dsigma_dB_osigma0 = 0.0;
    a0 = pm * sigma0/sqrt_1_A2_B2;
    dsigma_dA_osigma1 = 0.0;
    dsigma_dB_osigma1 = 0.0;
    a1 = pm * box->CI->s[1];
    da0_dx = a0;  /* mark as non-zero */
    da0_dy = a0;  /* mark as non-zero */
    da0_dz = a0;  /* mark as non-zero */
    da1_dx = 0.0;
    da1_dy = 0.0;
    da1_dz = 0.0;
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    dsigma_dA_osigma0 = 0.0;
    dsigma_dB_osigma0 = 0.0;
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    dsigma_dA_osigma1 = 0.0;
    dsigma_dB_osigma1 = 0.0;
    a1 = pm * sigma1/sqrt_1_A2_B2;
    da0_dx = 0.0;
    da0_dy = 0.0;
    da0_dz = 0.0;
    da1_dx = a1;  /* mark as non-zero */
    da1_dy = a1;  /* mark as non-zero */
    da1_dz = a1;  /* mark as non-zero */
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    sigma0 = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    dsigma_dA_osigma0 = 0.0;
    dsigma_dB_osigma0 = 0.0;
    sigma1 = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    dsigma_dA_osigma1 = 0.0;
    dsigma_dB_osigma1 = 0.0;
    a0 = pm * sigma0/sqrt_1_A2_B2;
    a1 = pm * sigma1/sqrt_1_A2_B2;
    da0_dx = a0; /* mark as non-zero */
    da0_dy = a0; /* mark as non-zero */
    da0_dz = a0; /* mark as non-zero */
    da1_dx = a1; /* mark as non-zero */
    da1_dy = a1; /* mark as non-zero */
    da1_dz = a1; /* mark as non-zero */
  }
  else errorexit("unknown type");

  /* compute coord trafo for each domain */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc)
       (1 + A^2 + B^2) rx^2 = rx^2 + ry^2 + rz^2
       1/sqrt(1 + A^2 + B^2) = |rx|/sqrt(rx^2 + ry^2 + rz^2) = |rx|/rc */
    rx = (a1-a0)*lam + a0;
    *x = xc + rx;
    *y = yc + A*rx;
    *z = zc + B*rx;
    ry = *y-yc;
    rz = *z-zc;
    rx2 = rx*rx;
    ry2 = ry*ry;
    rz2 = rz*rz;
    rc2= rx2 + ry2 + rz2;
    /* da0_dx = (1/rx - rx/rc^2)*a0 + ((dsigma0/dx)/sigma0)*a0 */
    /* da0_dy =       - ry/rc^2 *a0 + ((dsigma0/dy)/sigma0)*a0 */
    mx = (1.0/rx - rx/rc2);
    my = (       - ry/rc2);
    mz = (       - rz/rc2);
    /* dsigma0/dx = -(ry/rx2)*dsigma0_dA - (rz/rx2)*dsigma0_dB */
    /* dsigma0/dy =    (1/rx)*dsigma0_dA + 0 */
    mx0 = -dsigma_dA_osigma0*(ry/rx2) - dsigma_dB_osigma0*(rz/rx2);
    my0 =  dsigma_dA_osigma0/rx;
    mz0 =                               dsigma_dB_osigma0/rx;
    mx1 = -dsigma_dA_osigma1*(ry/rx2) - dsigma_dB_osigma1*(rz/rx2);
    my1 =  dsigma_dA_osigma1/rx;
    mz1 =                               dsigma_dB_osigma1/rx;
    da0_dx *= (mx + mx0);
    da0_dy *= (my + my0);
    da0_dz *= (mz + mz0);
    da1_dx *= (mx + mx1);
    da1_dy *= (my + my1);
    da1_dz *= (mz + mz1);
    /* lam = (rx-a0)/(a1-a0);
       A   = ry/rx;
       B   = rz/rx;
       lam = (rx-a0(x,y,z))/(a1(x,y,z)-a0(x,y,z)) */
    *dlamdx = (1.0 - da0_dx -(da1_dx-da0_dx)*lam)/(a1-a0);
    *dlamdy = (    - da0_dy -(da1_dy-da0_dy)*lam)/(a1-a0);
    *dlamdz = (    - da0_dz -(da1_dz-da0_dz)*lam)/(a1-a0);
    *dAdx = -ry/rx2;
    *dAdy = 1.0/rx;
    *dAdz = 0.0;
    *dBdx = -rz/rx2;
    *dBdy = 0.0;
    *dBdz = *dAdy;
  }
  else if(dir==2)
  {
    ry = (a1-a0)*lam + a0;
    *y = yc + ry;
    *x = xc + A*ry;
    *z = zc + B*ry;
    rx = *x-xc;
    rz = *z-zc;
    rx2 = rx*rx;
    ry2 = ry*ry;
    rz2 = rz*rz;
    rc2= rx2 + ry2 + rz2;
    mx = (       - rx/rc2);
    my = (1.0/ry - ry/rc2);
    mz = (       - rz/rc2);
    mx0 =  dsigma_dA_osigma0/ry;
    my0 = -dsigma_dA_osigma0*(rx/ry2) - dsigma_dB_osigma0*(rz/ry2);
    mz0 =                               dsigma_dB_osigma0/ry;
    mx1 =  dsigma_dA_osigma1/ry;
    my1 = -dsigma_dA_osigma1*(rx/ry2) - dsigma_dB_osigma1*(rz/ry2);
    mz1 =                               dsigma_dB_osigma1/ry;
    da0_dx *= (mx + mx0);
    da0_dy *= (my + my0);
    da0_dz *= (mz + mz0);
    da1_dx *= (mx + mx1);
    da1_dy *= (my + my1);
    da1_dz *= (mz + mz1);
    /* lam = (ry-a0)/(a1-a0); 
       A   = rx/ry;
       B   = rz/ry; */
    *dlamdx = (    - da0_dx -(da1_dx-da0_dx)*lam)/(a1-a0);
    *dlamdy = (1.0 - da0_dy -(da1_dy-da0_dy)*lam)/(a1-a0);
    *dlamdz = (    - da0_dz -(da1_dz-da0_dz)*lam)/(a1-a0);
    *dAdx = 1.0/ry;
    *dAdy = -rx/ry2;
    *dAdz = 0.0;
    *dBdx = 0.0;
    *dBdy = -rz/ry2;
    *dBdz = *dAdx;
  }
  else /* dir==3 */
  {
    rz = (a1-a0)*lam + a0;
    *z = zc + rz;
    *y = yc + A*rz;
    *x = xc + B*rz;
    rx = *x-xc;
    ry = *y-yc;
    rx2 = rx*rx;
    ry2 = ry*ry;
    rz2 = rz*rz;
    rc2= rx2 + ry2 + rz2;
    mx = (       - rx/rc2);
    my = (       - ry/rc2);
    mz = (1.0/rz - rz/rc2);
    mx0 =                               dsigma_dB_osigma0/rz;
    my0 =  dsigma_dA_osigma0/rz;
    mz0 = -dsigma_dA_osigma0*(ry/rz2) - dsigma_dB_osigma0*(rx/rz2);
    mx1 =                               dsigma_dB_osigma1/rz;
    my1 =  dsigma_dA_osigma1/rz;
    mz1 = -dsigma_dA_osigma1*(ry/rz2) - dsigma_dB_osigma1*(rx/rz2);
    da0_dx *= (mx + mx0);
    da0_dy *= (my + my0);
    da0_dz *= (mz + mz0);
    da1_dx *= (mx + mx1);
    da1_dy *= (my + my1);
    da1_dz *= (mz + mz1);
    /* lam = (rz-a0)/(a1-a0);
       A   = ry/rz;
       B   = rx/rz; */
    *dlamdx = (    - da0_dx -(da1_dx-da0_dx)*lam)/(a1-a0);
    *dlamdy = (    - da0_dy -(da1_dy-da0_dy)*lam)/(a1-a0);
    *dlamdz = (1.0 - da0_dz -(da1_dz-da0_dz)*lam)/(a1-a0);
    *dAdx = 0.0;
    *dAdy = 1.0/rz;
    *dAdz = -ry/rz2;
    *dBdx = *dAdy;
    *dBdy = 0.0;
    *dBdz = -rx/rz2;
  }
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


#define CALL_dlamAB_dxyz_CubSph \
  tBox *box = (tBox *) aux;\
  double x,y,z;\
  double dlamdx,dlamdy,dlamdz, dAdx,dAdy,dAdz, dBdx,dBdy,dBdz;\
  dlamAB_dxyz_CubSph(box, ind, lam,A,B, &x,&y,&z,\
                     &dlamdx, &dlamdy, &dlamdz,\
                     &dAdx,   &dAdy,   &dAdz,\
                     &dBdx,   &dBdy,   &dBdz)
#define RET_dlam_dx CALL_dlamAB_dxyz_CubSph; return dlamdx
#define RET_dlam_dy CALL_dlamAB_dxyz_CubSph; return dlamdy
#define RET_dlam_dz CALL_dlamAB_dxyz_CubSph; return dlamdz
#define RET_dA_dx CALL_dlamAB_dxyz_CubSph; return dAdx
#define RET_dA_dy CALL_dlamAB_dxyz_CubSph; return dAdy
#define RET_dA_dz CALL_dlamAB_dxyz_CubSph; return dAdz
#define RET_dB_dx CALL_dlamAB_dxyz_CubSph; return dBdx
#define RET_dB_dy CALL_dlamAB_dxyz_CubSph; return dBdy
#define RET_dB_dz CALL_dlamAB_dxyz_CubSph; return dBdz

/* Derivs for Cubed Spheres */
double dlam_dx_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dlam_dx; }
double dlam_dy_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dlam_dy; }
double dlam_dz_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dlam_dz; }
double dA_dx_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dA_dx; }
double dA_dy_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dA_dy; }
double dA_dz_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dA_dz; }
double dB_dx_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dB_dx; }
double dB_dy_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dB_dy; }
double dB_dz_CubedSphere(void *aux, int ind, double lam, double A, double B)
{  RET_dB_dz; }




/* Functions to get sigma0/1 for cubed spheres from box->CI->iSurf[0/1].
   box->CI->iSurf[0/1] is supposed to contain the index of the var where
   we store sigma.
   Here si is 0/1, and ind is the index of the point where we need sigma,
   and A,B are the coordinates where we need sigma */
double sigma_CubSph(tBox *box, int si, int ind, double A, double B)
{
  int isig = box->CI->iSurf[si]; /* get index of var with sigma */

  /* if var index does not seem set, use sigma = box->CI->s[si] */
  if(isig<1) return box->CI->s[si];

  /* if we are on a grid point (ind>=0), use value at this point */
  if(ind>=0) return box->v[isig][ind];
  else /* we need to interpolate */
  {
    return 0;// FIXME: interpolate_isig_to_point_given by A,B(box, isig, B,phi, 0);
  }
}
/* dsigma/dA and dsigma/dB need to be computed as well. */
