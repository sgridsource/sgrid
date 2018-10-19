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
int xyz_of_lamAB_CubSph(tBox *box, int ind, double lam, double A, double B,
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
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    a0 = pm * sigma0/sqrt_1_A2_B2;
    a1 = pm * box->CI->s[1];
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
    a1 = pm * sigma1/sqrt_1_A2_B2;
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
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
  return 0;
}

/* compute inverse sphered cube coord trafo. */
/* Note: lam \in [0,1], (A,B) \in [-1,1] */
int lamAB_of_xyz_CubSph(tBox *box, int ind, double x, double y, double z,
                        double *lam, double *A, double *B)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int type = box->CI->type;
  int dir, p;
  double pm, rx,ry,rz,rc, xc,yc,zc, a0,a1, sigma0,sigma1;
  double oosqrt_1_A2_B2;
  int stat=0;

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

  /* get A,B for each domain */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc)
       (1 + A^2 + B^2) rx^2 = rx^2 + ry^2 + rz^2
       1/sqrt(1 + A^2 + B^2) = |rx|/sqrt(rx^2 + ry^2 + rz^2) = |rx|/rc */
    *A   = ry/rx;
    *B   = rz/rx;
  }
  else if(dir==2)
  {
    *A   = rx/ry;
    *B   = rz/ry;
  }
  else /* dir==3 */
  {
    *A   = ry/rz;
    *B   = rx/rz;
  }

  /* check if A is out of range */
  if( *A < box->bbox[2] || box->bbox[3] < *A )
  {
    /* put A within range if we are just a bit out */
    if(dequal(*A,box->bbox[2]))      *A = box->bbox[2];
    else if(dequal(*A,box->bbox[3])) *A = box->bbox[3];
    else stat=-1;
  }
  /* check if B is out of range */
  if( *B < box->bbox[4] || box->bbox[5] < *B )
  {
    /* put B within range if we are just a bit out */
    if(dequal(*B,box->bbox[4]))      *B = box->bbox[4];
    else if(dequal(*B,box->bbox[5])) *B = box->bbox[5];
    else stat=-1;
  }

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
    sigma0 = CubedSphere_sigma(box, 0, ind, *A,*B);
    a0 = pm * sigma0*oosqrt_1_A2_B2;
    a1 = pm * box->CI->s[1];
    if(ind<0) stat*=2; /* we want signal that now we interpolated sigma */
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = pm * box->CI->s[0];
    sigma1 = CubedSphere_sigma(box, 1, ind, *A,*B);
    a1 = pm * sigma1*oosqrt_1_A2_B2;
    if(ind<0) stat*=2; /* we want signal that now we interpolated sigma */
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, *A,*B);
    sigma1 = CubedSphere_sigma(box, 1, ind, *A,*B);
    a0 = pm * sigma0*oosqrt_1_A2_B2;
    a1 = pm * sigma1*oosqrt_1_A2_B2;
    if(ind<0) stat*=2; /* we want signal that now we interpolated sigma */
  }
  else errorexit("unknown type");

  /* complete coord trafo for each domain and get lam */
  if(dir==1)
  { /* lam = (x-xc-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc)
       (1 + A^2 + B^2) rx^2 = rx^2 + ry^2 + rz^2
       1/sqrt(1 + A^2 + B^2) = |rx|/sqrt(rx^2 + ry^2 + rz^2) = |rx|/rc */
    *lam = (rx-a0)/(a1-a0);
  }
  else if(dir==2)
  {
    *lam = (ry-a0)/(a1-a0); 
  }
  else /* dir==3 */
  {
    *lam = (rz-a0)/(a1-a0);
  }
//if(!finite(*lam))
//{
//printf("box->b=%d ind=%d x=%g y=%g z=%g  lam=%g A=%g B=%g\n",
//box->b,ind, x,y,z, *lam,*A,*B);
//printf("sigma0=%g sigma1=%g a0=%g a1=%g\n", sigma0,sigma1, a0,a1);
//}
  return stat;
}

/* compute derivs of inverse trafo */
int dlamAB_dxyz_CubSph(tBox *box, int ind, double lam, double A, double B,
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
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    dsigma_dA_osigma0 = CubedSphere_dsigma_dA(box, 0, ind, A,B)/sigma0;
    dsigma_dB_osigma0 = CubedSphere_dsigma_dB(box, 0, ind, A,B)/sigma0;
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
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
    dsigma_dA_osigma1 = CubedSphere_dsigma_dA(box, 1, ind, A,B)/sigma1;
    dsigma_dB_osigma1 = CubedSphere_dsigma_dB(box, 1, ind, A,B)/sigma1;
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
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    dsigma_dA_osigma0 = CubedSphere_dsigma_dA(box, 0, ind, A,B)/sigma0;
    dsigma_dB_osigma0 = CubedSphere_dsigma_dB(box, 0, ind, A,B)/sigma0;
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
    dsigma_dA_osigma1 = CubedSphere_dsigma_dA(box, 1, ind, A,B)/sigma1;
    dsigma_dB_osigma1 = CubedSphere_dsigma_dB(box, 1, ind, A,B)/sigma1;
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
  return 0;
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



/* interpolate ibox->v[isig] */
/* The use of Temp1 is not thread safe!
   If tsafe=1 it is thread safe, because we allocate and the free
   memory for each call. This is slower though. */
double interpolate_isig_in_plane1_at_p(tBox *box, int isig,
                                       int p, double A, double B, int tsafe)
{
  double interp;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double *c;

  if(tsafe)  c = malloc(n1*n2*n3 * sizeof(double));
  else       c = box->v[Ind("Temp1")];

  /* interpolate in plane1 at index p */
  spec_Coeffs_inplaneN(box, 1,p, box->v[isig], c);
  interp = spec_interpolate_inplaneN(box, 1,p, c, A,B);
//if(!finite(interp))
//{
//printf("interp=%g box->b=%d p=%d A=%.18g B=%.18g\n", interp, box->b, p, A,B);
////errorexit("interp is not finite");
//}
  if(tsafe) free(c);
  return interp;
}

/* Functions to get sigma0/1 for cubed spheres at any point A,B.
   If box->CI->FSurf[si] contains a function that returns sigma it uses it,
   if not we interpolate using the vars in CI->iSurf[si].
   Here si is 0/1. */
double CubedSphere_sigma_AB(tBox *box, int si, double A, double B)
{
  tCoordInfo *CI = box->CI;
  int n1, p, isig;

  if(CI->useF && CI->FSurf[si] != NULL)
    return       CI->FSurf[si](box, si, A,B);

  /* if useF=0 interpolate box vars */
  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */
  isig = CI->iSurf[si]; /* get index of var with sigma */
  return interpolate_isig_in_plane1_at_p(box, isig, p, A,B, 1);
}
/* dsigma/dA and dsigma/dB are derivs of sigma. */
double CubedSphere_dsigma_dA_AB(tBox *box, int si, double A, double B)
{
  tCoordInfo *CI = box->CI;
  int n1, p, isig;

  if(CI->useF && CI->dFSurfdX[si][2] != NULL)
    return       CI->dFSurfdX[si][2](box, si, A,B);

  /* if useF=0 interpolate box vars */
  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */
  isig = CI->idSurfdX[si][2]; /* get index of var with dsigma/dA */
  return interpolate_isig_in_plane1_at_p(box, isig, p, A,B, 1);
}
double CubedSphere_dsigma_dB_AB(tBox *box, int si, double A, double B)
{
  tCoordInfo *CI = box->CI;
  int n1, p, isig;

  if(CI->useF && CI->dFSurfdX[si][3] != NULL)
    return       CI->dFSurfdX[si][3](box, si, A,B);

  /* if useF=0 interpolate box vars */
  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */
  isig = CI->idSurfdX[si][3]; /* get index of var with dsigma/dB */
  return interpolate_isig_in_plane1_at_p(box, isig, p, A,B, 1);
}


/* Functions to get sigma0/1 for cubed spheres from box->CI->iSurf[0/1].
   box->CI->iSurf[0/1] is supposed to contain the index of the var where
   we store sigma.
   Here si is 0/1, and ind is the index of the point where we need sigma,
   and A,B are the coordinates where we need sigma */
double CubedSphere_sigma(tBox *box, int si, int ind, double A, double B)
{
  int isig = box->CI->iSurf[si]; /* get index of var with sigma */
  int p, n1;

  /* if var index does not seem set, use sigma = box->CI->s[si] */
  if(isig<1) return box->CI->s[si];

  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */

  /* if we are on a grid point (ind>=0), use value at this j,k, but i=p */
  if(ind>=0)
  {
    int n2= box->n2;
    int k = kOfInd_n1n2(ind,n1,n2);
    int j = jOfInd_n1n2_k(ind,n1,n2,k);
    int ijk = Index(p,j,k);
    return box->v[isig][ijk];
  }
  else /* we need to interpolate */
    return CubedSphere_sigma_AB(box, si, A,B);
}
/* dsigma/dA and dsigma/dB are derivs of sigma. */
double CubedSphere_dsigma_dA(tBox *box, int si, int ind, double A, double B)
{
  int isig = box->CI->idSurfdX[si][2]; /* get index of var with dsigma/dA */
  int p, n1;

  /* if var index does not seem set, use dsigma = 0 */
  if(isig<1) return 0.0;

  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */

  /* if we are on a grid point (ind>=0), use value at this j,k, but i=p */
  if(ind>=0)
  {
    int n2= box->n2;
    int k = kOfInd_n1n2(ind,n1,n2);
    int j = jOfInd_n1n2_k(ind,n1,n2,k);
    int ijk = Index(p,j,k);
    return box->v[isig][ijk];
  }
  else /* we need to interpolate */
    return CubedSphere_dsigma_dA_AB(box, si, A,B);
}
double CubedSphere_dsigma_dB(tBox *box, int si, int ind, double A, double B)
{
  int isig = box->CI->idSurfdX[si][3]; /* get index of var with dsigma/dB */
  int p, n1;

  /* if var index does not seem set, use dsigma = 0 */
  if(isig<1) return 0.0;

  n1 = box->n1;
  p = (n1-1) * (si==1);   /* set plane index p */

  /* if we are on a grid point (ind>=0), use value at this j,k, but i=p */
  if(ind>=0)
  {
    int n2= box->n2;
    int k = kOfInd_n1n2(ind,n1,n2);
    int j = jOfInd_n1n2_k(ind,n1,n2,k);
    int ijk = Index(p,j,k);
    return box->v[isig][ijk];
  }
  else /* we need to interpolate */
    return CubedSphere_dsigma_dB_AB(box, si, A,B);
}


/************************************************************************/
/* some functions to relate r and lam                                   */
/************************************************************************/
/*
This is for a trafo from (lam,A,B) to (r,A,B):
Recall in e.g. dom1: 
x = (a1-a0)*lam + a0,  y = A*x,  z = B*x
==> r^2 = (1 + A^2 + B^2)*x^2 = (1 + A^2 + B^2) * ((a1-a0)*lam + a0)
def: sig1 = a1 * sqrt(1 + A^2 + B^2),  sig0 = a0 * sqrt(1 + A^2 + B^2),  
Thus:  r = (sig1-sig0)*lam + sig0
Def: L = sig1-sig0,   so: r = L lam + sig0
*/
double r_of_lam_sig0sig1(double lam, double sig0, double sig1)
{
  double L = sig1-sig0;
  return L*lam + sig0;
}
double lam_of_r_sig0sig1(double r, double sig0, double sig1)
{
  double L = sig1-sig0;
  return (r-sig0)/L;
}
/* Derivatives */
double dr_dlam_of_lam_sig0sig1(double lam, double sig0, double sig1)
{
  double L = sig1-sig0;
  return L;
}
double dlam_dr_of_lam_sig0sig1(double lam, double sig0, double sig1)
{
  double L = sig1-sig0;
  return 1./L;
}
/* r and dr/dlam under Coordtrafo (lam,A,B) -> (r,A,B) */
int r_dr_dlam_of_lamAB_CubSph(tBox *box, int ind, double lam,
                              double A, double B, double *r, double *drdlam)
{
  int type = box->CI->type;
  double sig0,sig1;
  double sqrt_1_A2_B2 = sqrt(1.0 + A*A + B*B);
  double a0, a1;

  /* check type of trafo */
  if(type==PyramidFrustum)
  {
    /* this gives a pyramid frustum */
    a0 = box->CI->s[0];
    a1 = box->CI->s[1];
    sig0 = a0 * sqrt_1_A2_B2;
    sig1 = a1 * sqrt_1_A2_B2;
  }
  else if(type==innerCubedSphere)
  {
    /* this gives a cubed sphere piece where the inner surface is curved
       and the outer surface is flat */
    sig0 = CubedSphere_sigma(box, 0, ind, A,B);
    a1 = box->CI->s[1];
    sig1 = a1 * sqrt_1_A2_B2;
  }
  else if(type==outerCubedSphere)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    a0 = box->CI->s[0];
    sig0 = a0 * sqrt_1_A2_B2;
    sig1 = CubedSphere_sigma(box, 1, ind, A,B);
  }
  else if(type==CubedShell)
  {
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sig0 = CubedSphere_sigma(box, 0, ind, A,B);
    sig1 = CubedSphere_sigma(box, 1, ind, A,B);
  }
  else errorexit("unknown type");

  /* set r and dr/dlam */
  *r      = r_of_lam_sig0sig1(lam, sig0, sig1);
  *drdlam = dr_dlam_of_lam_sig0sig1(lam, sig0, sig1);
  return 0;
}

/************************************************************************/
/* some functions to relate Theta,Phi and A,B                           */
/************************************************************************/
/* Arg_plus(X,Y) gives arctan(Y/X), it returns a value in (-PI,PI] */

/* get Theta, Phi from A,B in one box */
int ThetaPhi_of_AB_CubSph(tBox *box, double A, double B,
                          double *Theta, double *Phi)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int dir, p;
  double pm;

  /* get direction, pm */
  dir = domain/2 + 1;
  p  = 2*(domain%2) - 1; /* p  = +/- 1 */
  pm = p;                /* pm = +/- 1.0 */
  //z1 = (1. - pm)*0.5;  /* z1 = 0 or 1 */

  if(dir==1)
  { /* A = ry/rx;   B = rz/rx;
       rx = r*cos(Phi)*sin(Theta)
       ry = r*sin(Phi)*sin(Theta)
       rz = r*cos(Theta)
       ry/rx = tan(Phi);                rz/rx = 1/(tan(Theta)*cos(Phi));
       rx/rz = cos(Phi)*tan(Theta)      rz/ry = 1/(tan(Theta)*sin(Phi));
       Phi = atan(ry/rx);
       A = ry/rx = tan(Phi);
       B = rz/rx = 1/(tan(Theta)*cos(Phi));  => tan(Theta) = 1/(B*cos(Phi)) */
    *Phi = Arg_plus(pm, pm*A);
    *Theta = Arg_plus(B*cos(*Phi), 1.);
  }
  else if(dir==2)
  { /* A = rx/ry;   B = rz/ry; */
    *Phi = Arg_plus(pm*A, pm);
    *Theta = Arg_plus(B*sin(*Phi), 1.);
  }
  else /* dir==3 */
  { /* A = ry/rz;   B = rx/rz;   A/B = ry/rx;
       A = sin(Phi)*tan(Theta)
       B = cos(Phi)*tan(Theta)  => A/B = tan(Phi) => Phi = atan(A/B)= Arg(B,A)
       s^2 = A^2 + B^2 = tan(Theta)^2
        => Theta = atan(sqrt(A^2 + B^2)) = Arg(1, sqrt(A^2 + B^2))   */
    *Phi = Arg_plus(pm*B, pm*A);
    *Theta = Arg_plus(pm, sqrt(A*A + B*B));
    /* Derivs of Phi, Theta using dArgdx and dArgdy:
       dPhi/dA = B/sqrt(A^2 + B^2)   dPhi/dB = -A/sqrt(A^2 + B^2) 
       dPhi/dA = cos(Phi)            dPhi/dB = -sin(Phi)
       Theta =  Arg(1, s) => dTheta/ds = 1/(1+s^2)
       dTheta/dA = dTheta/ds ds/dA
       dTheta/dA = (1/(1+s^2)) A/sqrt(A^2 + B^2)
       dTheta/dB = (1/(1+s^2)) B/sqrt(A^2 + B^2)
       dTheta/dA = (1/(1+s^2)) sin(Phi)
       dTheta/dB = (1/(1+s^2)) cos(Phi)
       I.e. derivs are not continous across north pole!!! */
  }

  return 0;
}
/* get Theta, Phi, dTheta/dA, dTheta/dB, dPhi/dA, dPhi/dB, 
   from A,B in one box */
int ThetaPhi_dThetaPhidAB_of_AB_CubSph(tBox *box, double A, double B,
                                       double *Theta,    double *Phi,
                                       double *dThetadA, double *dThetadB,
                                       double *dPhidA,   double *dPhidB)
{
  int domain = box->CI->dom;  /* get domain and type info from box */
  int dir, p;
  double pm;

  /* get direction, pm */
  dir = domain/2 + 1;
  p  = 2*(domain%2) - 1; /* p  = +/- 1 */
  pm = p;                /* pm = +/- 1.0 */
  //z1 = (1. - pm)*0.5;  /* z1 = 0 or 1 */

  if(dir==1)
  { /* A = ry/rx;   B = rz/rx;
       rx = r*cos(Phi)*sin(Theta)
       ry = r*sin(Phi)*sin(Theta)
       rz = r*cos(Theta)
       ry/rx = tan(Phi);                rz/rx = 1/(tan(Theta)*cos(Phi));
       rx/rz = cos(Phi)*tan(Theta)      rz/ry = 1/(tan(Theta)*sin(Phi));
       Phi = atan(ry/rx);
       A = ry/rx = tan(Phi);
       B = rz/rx = 1/(tan(Theta)*cos(Phi));  => tan(Theta) = 1/(B*cos(Phi)) */
    *Phi = Arg_plus(pm, pm*A);
    *Theta = Arg_plus(B*cos(*Phi), 1.);
    /* note:  1+tan^2(phi) = 1/cos^2(phi) ==> cos(phi) = 1/sqrt(1+tan^2(phi)) */

    /* Derivs */
    *dPhidA   = dArgdy(pm, pm*A)*pm;
    *dPhidB   = 0.;
    *dThetadA = dArgdx(B*cos(*Phi), 1.) * B*(-sin(*Phi)) * (*dPhidA);
    *dThetadB = dArgdx(B*cos(*Phi), 1.) *
                ( cos(*Phi) + B*(-sin(*Phi)) * (*dPhidB) );
  }
  else if(dir==2)
  { /* A = rx/ry;   B = rz/ry; */
    *Phi = Arg_plus(pm*A, pm);
    *Theta = Arg_plus(B*sin(*Phi), 1.);

    /* Derivs */
    *dPhidA   = dArgdx(pm*A, pm)*pm;
    *dPhidB   = 0.;
    *dThetadA = dArgdx(B*sin(*Phi), 1.) * B*(cos(*Phi)) * (*dPhidA);
    *dThetadB = dArgdx(B*sin(*Phi), 1.) *
                ( sin(*Phi) + B*(cos(*Phi)) * (*dPhidB) );
  }
  else /* dir==3 */
  { /* A = ry/rz;   B = rx/rz;   A/B = ry/rx;
       A = sin(Phi)*tan(Theta)
       B = cos(Phi)*tan(Theta)  => A/B = tan(Phi) => Phi = atan(A/B)= Arg(B,A)
       s^2 = A^2 + B^2 = tan(Theta)^2
        => Theta = atan(sqrt(A^2 + B^2)) = Arg(1, sqrt(A^2 + B^2))   */
    double sqrtA2B2 = sqrt(A*A + B*B);
    double dsqrtA2B2_dA = A/sqrtA2B2; 
    double dsqrtA2B2_dB = B/sqrtA2B2; 

    *Phi = Arg_plus(pm*B, pm*A);
    *Theta = Arg_plus(pm, sqrtA2B2);
    /* Derivs of Phi, Theta using dArgdx and dArgdy:
       dPhi/dA = B/sqrt(A^2 + B^2)   dPhi/dB = -A/sqrt(A^2 + B^2) 
       dPhi/dA = cos(Phi)            dPhi/dB = -sin(Phi)
       Theta =  Arg(1, s) => dTheta/ds = 1/(1+s^2)
       dTheta/dA = dTheta/ds ds/dA
       dTheta/dA = (1/(1+s^2)) A/sqrt(A^2 + B^2)
       dTheta/dB = (1/(1+s^2)) B/sqrt(A^2 + B^2)
       dTheta/dA = (1/(1+s^2)) sin(Phi)
       dTheta/dB = (1/(1+s^2)) cos(Phi)
       I.e. derivs are not continous across north pole!!! */

    /* Derivs */
    *dPhidA = dArgdy(pm*B, pm*A)*pm;
    *dPhidB = dArgdx(pm*B, pm*A)*pm;
    *dThetadA = dArgdy(pm, sqrtA2B2) * dsqrtA2B2_dA;
    *dThetadB = dArgdy(pm, sqrtA2B2) * dsqrtA2B2_dB;
  }

  return 0;
}



/************************************************************************/
/* some functions to stretch a cubed sphere to a large radius sigma0    */
/************************************************************************/
/*
Recall in e.g. dom1: r = (sig1-sig0)*lam + sig0
Def: L = sig1-sig0,   so: r = L lam + sig0
We now replace lam by rho where
rho = (sig1/L)*( 1.0 - sig0/(L*lam + sig0) )
lam = (sig0/L)*( sig1/(sig1 - L*rho) - 1.0 )
Then
rho = (sig1/L)*( 1 - sig0/r )
*/
double rho_of_lam_sig0sig1(double lam, double sig0, double sig1)
{
  double L = sig1-sig0;
  return (sig1/L)*( 1.0 - sig0/(L*lam + sig0) );
}
double lam_of_rho_sig0sig1(double rho, double sig0, double sig1)
{
  double L = sig1-sig0;
  return (sig0/L)*( sig1/(sig1 - L*rho) - 1.0 );
}

/* Derivatives */
double drho_dlam_of_rho_sig0sig1(double rho, double sig0, double sig1)
{
  double L = sig1-sig0;
  double d = sig1 - L*rho;
  return d*d/(sig0*sig1);
}

/* Coord. trafos of for stretched Cubed Spheres */
/* Coordtrafo (rho,A,B) -> (x,y,z) */
int xyz_of_rhoAB_CubSph(tBox *box, int ind, double rho, double A, double B,
                        double *x, double *y, double *z)
{
  int type = box->CI->type;
  double sigma0,sigma1, lam;

  /* check type of trafo */
  if(type==CubedShell)
  {
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
  }
  else errorexit("xyz_of_rhoAB_CubSph works only with CubedShell");

  /* get lam, and then set x,y,z */
  lam = lam_of_rho_sig0sig1(rho, sigma0,sigma1);
  return xyz_of_lamAB_CubSph(box, ind, lam,A,B, x,y,z);
}

/* inverse (x,y,z) -> (rho,A,B) */
int rhoAB_of_xyz_CubSph(tBox *box, int ind, double x, double y, double z,
                        double *rho, double *A, double *B)
{
  int type = box->CI->type;
  int stat;
  double sigma0,sigma1, lam;

  /* get lam,A,B from x,y,z */
  stat = lamAB_of_xyz_CubSph(box, ind, x,y,z, &lam,A,B);

  /* check type of trafo */
  if(type==CubedShell)
  {
    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, *A,*B);
    sigma1 = CubedSphere_sigma(box, 1, ind, *A,*B);
  }
  else errorexit("works only for CubedShell");

  /* get rho from lambda */
  *rho = rho_of_lam_sig0sig1(lam, sigma0,sigma1);
  return stat;
}

/* compute derivs of inverse trafo */
int drhoAB_dxyz_CubSph(tBox *box, int ind, double rho, double A, double B,
                       double *x, double *y, double *z,
                       double *drhodx, double *drhody, double *drhodz,
                       double *dAdx,   double *dAdy,   double *dAdz,
                       double *dBdx,   double *dBdy,   double *dBdz)
{
  int type = box->CI->type;
  double sigma0,sigma1, lam, dlamdx,dlamdy,dlamdz, drho_dlam;

  /* check type of trafo */
  if(type==CubedShell)
  {
    int isig0 = box->CI->iSurf[0]; /* get index of var with sigma */
    int isig1 = box->CI->iSurf[1]; /* get index of var with sigma */
    if(isig0>0 || isig1>0)
      errorexit("drhoAB_dxyz_CubSph works only with constant sigmas");

    /* this gives a cubed sphere where both inner outer surfaces are curved */
    sigma0 = CubedSphere_sigma(box, 0, ind, A,B);
    sigma1 = CubedSphere_sigma(box, 1, ind, A,B);
  }
  else errorexit("works only for CubedShell");

  /* get lam and derivs of lam */
  lam = lam_of_rho_sig0sig1(rho, sigma0,sigma1);
  dlamAB_dxyz_CubSph(box, ind, lam,A,B, x,y,z, &dlamdx, &dlamdy, &dlamdz,
                                                  dAdx,    dAdy,    dAdz,
                                                  dBdx,    dBdy,    dBdz);
  /* get drho_dlam */
  drho_dlam = drho_dlam_of_rho_sig0sig1(rho, sigma0,sigma1);
  
  /* now set drhodx, drhody, drhodz */
  *drhodx = drho_dlam * dlamdx;
  *drhody = drho_dlam * dlamdy;
  *drhodz = drho_dlam * dlamdz;

  /* if sigmas are not const we need to add: */
  /* *drhodx += drho_dA * (*dAdx) + drho_dB * (*dBdx);
     *drhody += drho_dA * (*dAdy) + drho_dB * (*dBdy);
     *drhodz += drho_dA * (*dAdz) + drho_dB * (*dBdz); */
  /* Note: drho_dA depends on dsigma_{0/1}/dA */
  return 0;
}


#define CALL_xyz_of_rhoAB_CubSph \
  tBox *box = (tBox *) aux;\
  double x,y,z;\
  xyz_of_rhoAB_CubSph(box, ind, rho,A,B, &x,&y,&z)
#define RET2_x CALL_xyz_of_rhoAB_CubSph; return x
#define RET2_y CALL_xyz_of_rhoAB_CubSph; return y
#define RET2_z CALL_xyz_of_rhoAB_CubSph; return z

/* Coord. trafos of for stretched Cubed Spheres */
double x_of_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{ RET2_x; }
double y_of_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{ RET2_y; }
double z_of_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{ RET2_z; }


#define CALL_drhoAB_dxyz_CubSph \
  tBox *box = (tBox *) aux;\
  double x,y,z;\
  double drhodx,drhody,drhodz, dAdx,dAdy,dAdz, dBdx,dBdy,dBdz;\
  drhoAB_dxyz_CubSph(box, ind, rho,A,B, &x,&y,&z,\
                     &drhodx, &drhody, &drhodz,\
                     &dAdx,   &dAdy,   &dAdz,\
                     &dBdx,   &dBdy,   &dBdz)
#define RET2_drho_dx CALL_drhoAB_dxyz_CubSph; return drhodx
#define RET2_drho_dy CALL_drhoAB_dxyz_CubSph; return drhody
#define RET2_drho_dz CALL_drhoAB_dxyz_CubSph; return drhodz
#define RET2_dA_dx CALL_drhoAB_dxyz_CubSph; return dAdx
#define RET2_dA_dy CALL_drhoAB_dxyz_CubSph; return dAdy
#define RET2_dA_dz CALL_drhoAB_dxyz_CubSph; return dAdz
#define RET2_dB_dx CALL_drhoAB_dxyz_CubSph; return dBdx
#define RET2_dB_dy CALL_drhoAB_dxyz_CubSph; return dBdy
#define RET2_dB_dz CALL_drhoAB_dxyz_CubSph; return dBdz

/* Derivs for stretched Cubed Spheres */
double drho_dx_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_drho_dx; }
double drho_dy_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_drho_dy; }
double drho_dz_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_drho_dz; }
double dA_dx_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dA_dx; }
double dA_dy_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dA_dy; }
double dA_dz_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dA_dz; }
double dB_dx_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dB_dx; }
double dB_dy_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dB_dy; }
double dB_dz_sCubedSphere(void *aux, int ind, double rho, double A, double B)
{  RET2_dB_dz; }
