/* coordtrans_CubSph.c */
/* Wolfgang Tichy, Nov 2017 */
/* transformations for cubed spheres which really are sphered cubes */

#include "sgrid.h"

/*
We only show domains 0-3. Domains 4 and 5 are in z-dir and are not shown here:

    1 excised cube    1 excised sphere   1 excised cube
    in the middle     in the middle      in the middle
    cubical outside   cubical outside    spheroidal outside
y                                              ____
 ^                                           _/    \_ 
 |    ________          ________            /        \
 |   |\  3   /|        |\  3   /|          / \  3   / \
 |   | \ __ / |        | \ __ / |         /   \ __ /   \
 |   |0 |  | 1|        |0 /  \ 1|        |  0  |  |  1  |
 |   |  |__|  |        |  \__/  |        |     |__|     |
 |   | / 2  \ |        | / 2  \ |         \   /    \   /
 |   |/______\|        |/______\|          \ /  2   \ /
 |                                          \_      _/
 |                                            \____/ 
-|------------------------------------------------------------> x
       FLAT              CURV0                 CURV1
*/

/* Type of cubed sphere or rather sphered cube coord transform */
enum
{
  FLAT,     /* both inner & outer surfaces are flat */  
  CURV0,    /* inner surface is curved, but outer surface is flat */
  CURV1,    /* outer surface is curved, but inner surface is flat */
};


/* compute a sphered cube coord trafo. The type is enumerated above. */
/* Note: lam \in [0,1], (A,B) \in [-1,1] */
void xyz_of_lamAB_SphCube(tBox *box, int ind, int domain, int type,
                          double lam, double A, double B,
                          double *x, double *y, double *z)
{
  int dir, p;
  double pm, xc,yc,zc, a0,a1, sigma;

  /* get direction, pm, and center */
  dir = domain/2 + 1;
  p  = 2*(domain%2) - 1; /* p  = +/- 1 */
  pm = p;                /* pm = +/- 1.0 */
  xc = box->CI->xc[1];
  yc = box->CI->xc[2];
  zc = box->CI->xc[3];

  /* check type of trafo */
  if(type==FLAT)
  {
    /* this gives a pyramid frustum */
    a0 = pm * box->CI->s[0];
    a1 = pm * box->CI->s[1];
  }
  else if(type==CURV0)
  {
    /* this gives a cubed sphere piece where the inner surface is curved
       and the outer surface is flat */
    sigma = box->CI->s[0];  /* or get it from box->CI->iSurf[0] */
    a0 = pm * sigma/sqrt(1.0 + A*A + B*B);
    a1 = pm * box->CI->s[1];
  }
  else if(type==CURV1)
  {
    /* this gives a cubed sphere piece where the outer surface is curved
       and the inner one is flat */
    sigma = box->CI->s[1];  /* or get it from box->CI->iSurf[1] */
    a0 = pm * box->CI->s[0];
    a1 = pm * sigma/sqrt(1.0 + A*A + B*B);
  }
  else errorexit("unknown type");

  /* compute coord trafo for each domain */
  if(dir==1)
  { /* lam = (x-a0)/(a1-a0),  A = (y-yc)/(x-xc),  B = (z-zc)/(x-xc) */
    *x = (a1-a0)*lam + a0;
    *y = yc + A * (*x-xc);
    *z = zc + B * (*x-xc);
  }
  else if(dir==2)
  {
    *y = (a1-a0)*lam + a0;
    *x = xc + A * (*y-yc);
    *z = zc + B * (*y-yc);
  }
  else /* dir==3 */
  {
    *z = (a1-a0)*lam + a0;
    *y = yc + A * (*z-zc);
    *x = xc + B * (*z-zc);
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



#define CALL_xyz_of_lamAB_SphCube(dom, typ) \
  tBox *box = (tBox *) aux;\
  double x,y,z;\
  xyz_of_lamAB_SphCube(box, ind, dom, typ, lam,A,B, &x,&y,&z);
#define RETURN_x(dom,typ) CALL_xyz_of_lamAB_SphCube(dom,typ); return x
#define RETURN_y(dom,typ) CALL_xyz_of_lamAB_SphCube(dom,typ); return y
#define RETURN_z(dom,typ) CALL_xyz_of_lamAB_SphCube(dom,typ); return z


/* Coord. trafos of type=FLAT, i.e. a pyramid frustum */
double   x_of_PyraFrust0(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(0,FLAT); }
double   y_of_PyraFrust0(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(0,FLAT); }
double   z_of_PyraFrust0(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(0,FLAT); }

double   x_of_PyraFrust1(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(1,FLAT); }
double   y_of_PyraFrust1(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(1,FLAT); }
double   z_of_PyraFrust1(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(1,FLAT); }

double   x_of_PyraFrust2(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(2,FLAT); }
double   y_of_PyraFrust2(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(2,FLAT); }
double   z_of_PyraFrust2(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(2,FLAT); }

double   x_of_PyraFrust3(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(3,FLAT); }
double   y_of_PyraFrust3(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(3,FLAT); }
double   z_of_PyraFrust3(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(3,FLAT); }

double   x_of_PyraFrust4(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(4,FLAT); }
double   y_of_PyraFrust4(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(4,FLAT); }
double   z_of_PyraFrust4(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(4,FLAT); }

double   x_of_PyraFrust5(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(5,FLAT); }
double   y_of_PyraFrust5(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(5,FLAT); }
double   z_of_PyraFrust5(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(5,FLAT); }


/* Coord. trafos of type=CURV0, i.e. a cubed sphere where inner surface is curved */
double   x_of_innerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(0,CURV0); }
double   y_of_innerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(0,CURV0); }
double   z_of_innerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(0,CURV0); }

double   x_of_innerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(1,CURV0); }
double   y_of_innerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(1,CURV0); }
double   z_of_innerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(1,CURV0); }

double   x_of_innerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(2,CURV0); }
double   y_of_innerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(2,CURV0); }
double   z_of_innerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(2,CURV0); }

double   x_of_innerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(3,CURV0); }
double   y_of_innerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(3,CURV0); }
double   z_of_innerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(3,CURV0); }

double   x_of_innerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(4,CURV0); }
double   y_of_innerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(4,CURV0); }
double   z_of_innerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(4,CURV0); }

double   x_of_innerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(5,CURV0); }
double   y_of_innerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(5,CURV0); }
double   z_of_innerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(5,CURV0); }

/* Coord. trafos of type=CURV1, i.e. a cubed sphere where outer surface is curved */
double   x_of_outerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(0,CURV1); }
double   y_of_outerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(0,CURV1); }
double   z_of_outerCubSph0(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(0,CURV1); }

double   x_of_outerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(1,CURV1); }
double   y_of_outerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(1,CURV1); }
double   z_of_outerCubSph1(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(1,CURV1); }

double   x_of_outerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(2,CURV1); }
double   y_of_outerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(2,CURV1); }
double   z_of_outerCubSph2(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(2,CURV1); }

double   x_of_outerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(3,CURV1); }
double   y_of_outerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(3,CURV1); }
double   z_of_outerCubSph3(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(3,CURV1); }

double   x_of_outerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(4,CURV1); }
double   y_of_outerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(4,CURV1); }
double   z_of_outerCubSph4(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(4,CURV1); }

double   x_of_outerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_x(5,CURV1); }
double   y_of_outerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_y(5,CURV1); }
double   z_of_outerCubSph5(void *aux, int ind, double lam, double A, double B)
{ RETURN_z(5,CURV1); }
