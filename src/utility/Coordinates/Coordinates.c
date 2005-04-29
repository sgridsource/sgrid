/* Coordinates.c */
/* Wolfgang Tichy 4/2005 */

#include "sgrid.h"
#include "Coordinates.h"


/* initialize the coord transforms */
int init_CoordTransform_And_Derivs(tGrid *grid)
{
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    if( Getv(str, "Polar") )
    {
      box->x_of_X[1] = x_ofPolar;
      box->x_of_X[2] = y_ofPolar;
      box->x_of_X[3] = z_ofPolar;

      box->dX_dx[1][1] = drho_dx;
      box->dX_dx[1][2] = drho_dy;
      box->dX_dx[1][3] = zero_of_xyz;
      box->dX_dx[2][1] = dphi_dx;
      box->dX_dx[2][2] = dphi_dy;
      box->dX_dx[2][3] = zero_of_xyz;
      box->dX_dx[3][1] = zero_of_xyz;
      box->dX_dx[3][2] = zero_of_xyz;
      box->dX_dx[3][3] = one_of_xyz;

      //box->Sing_d_dx[1] = set_d_dx_at_rhoEQzero;
      //box->Sing_d_dx[2] = set_d_dy_at_rhoEQzero;
      /*
      box->dX_dxdx[1][1][1] = drho_dxdx;
      box->dX_dxdx[1][1][2] = drho_dxdy;
      box->dX_dxdx[1][1][3] = zero_of_xyz;
      box->dX_dxdx[1][2][2] = drho_dydy;
      box->dX_dxdx[1][2][3] = zero_of_xyz;
      box->dX_dxdx[1][3][3] = zero_of_xyz;

      box->dX_dxdx[2][1][1] = dphi_dxdx;
      box->dX_dxdx[2][1][2] = dphi_dxdy;
      box->dX_dxdx[2][1][3] = zero_of_xyz;
      box->dX_dxdx[2][2][2] = dphi_dydy;
      box->dX_dxdx[2][2][3] = zero_of_xyz;
      box->dX_dxdx[2][3][3] = zero_of_xyz;

      box->dX_dxdx[3][1][1] = zero_of_xyz;
      box->dX_dxdx[3][1][2] = zero_of_xyz;
      box->dX_dxdx[3][1][3] = zero_of_xyz;
      box->dX_dxdx[3][2][2] = zero_of_xyz;
      box->dX_dxdx[3][2][3] = zero_of_xyz;
      box->dX_dxdx[3][3][3] = zero_of_xyz;
      */
    }
    if( Getv(str, "PolarCE") )
    {
      box->x_of_X[1] = x_ofPolarCE;
      box->x_of_X[2] = y_ofPolarCE;
      box->x_of_X[3] = z_equals_Z;

      box->dX_dx[1][1] = drho_dx;
      box->dX_dx[1][2] = drho_dy;
      box->dX_dx[1][3] = zero_of_xyz;
      box->dX_dx[2][1] = dYPolarCE_dx;
      box->dX_dx[2][2] = dYPolarCE_dy;
      box->dX_dx[2][3] = zero_of_xyz;
      box->dX_dx[3][1] = zero_of_xyz;
      box->dX_dx[3][2] = zero_of_xyz;
      box->dX_dx[3][3] = one_of_xyz;

      box->Sing_d_dx[1] = NULL;
      box->Sing_d_dx[2] = NULL;
    }

    /* compute cartesian coordinates x,y,z from X,Y,Z */
    if( box->x_of_X[1] != NULL )
    {
      int ind;
      enablevar_inbox(box, var_x);
      enablevar_inbox(box, var_y);
      enablevar_inbox(box, var_z);

      forallpoints(box, ind)
      {
        double X = box->v[var_X][ind];
        double Y = box->v[var_Y][ind];
        double Z = box->v[var_Z][ind];
        box->v[var_x][ind] = box->x_of_X[1]((void *) box, X,Y,Z);
        box->v[var_y][ind] = box->x_of_X[2]((void *) box, X,Y,Z);
        box->v[var_z][ind] = box->x_of_X[3]((void *) box, X,Y,Z);
      }
    }
  
  }
  return 0;
}


/* Some trivial functions */
double zero_of_xyz(void *aux, double X, double Y, double Z)
{
  return 0.0;
}
double one_of_xyz(void *aux, double X, double Y, double Z)
{
  return 1.0;
}
double x_equals_X(void *aux, double X, double Y, double Z)
{
  return X;
}
double y_equals_Y(void *aux, double X, double Y, double Z)
{
  return Y;
}
double z_equals_Z(void *aux, double X, double Y, double Z)
{
  return Z;
}

/* start: Polar coordinates: */

/* Coord. trafos */
double x_ofPolar(void *aux, double rho, double phi, double Z)
{
  return rho*cos(phi);
}
double y_ofPolar(void *aux, double rho, double phi, double Z)
{
  return rho*sin(phi);
}
double z_ofPolar(void *aux, double rho, double phi, double Z)
{
  return Z;
}

/* first coord. derivs */
double drho_dx(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  
  if(rho2>0.0) return x/sqrt(rho2);
  else         return 1.0; /* result if we go along y=0 line */
}
double drho_dy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  if(rho2>0.0) return y/sqrt(rho2);
  else         return 1.0; /* result if we go along x=0 line */
}

double dphi_dx(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  if(rho2>0.0) return -y/(rho2);
  else         return 0.0; /* result if we go along y=0 line */
}
double dphi_dy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  if(rho2>0.0) return x/(rho2);
  else         return 0.0; /* result if we go along x=0 line */
}

/* functions to treat cartesian derivs at singular points */
void set_d_dx_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3)
{
  tBox *box = (tBox *) bo;
  double *vx = (double *) va;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int j,k;

  /* take deriv d/dx at phi=0 <=> y=0 and use it everywhere */
  for(k=0; k<n3; k++)
    for(j=1; j<n2; j++)
      vx[Index(0,j,k)] = vx[Index(0,0,k)];
}
void set_d_dy_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3)
{
  tBox *box = (tBox *) bo;
  double *vy = (double *) va;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int j,k;

  if(n2 % 4)
    errorexit("set_d_dy_at_rhoEQzero: box->n2 has to be divisible by 4.");

  /* take deriv d/dy at phi=pi/2 <=> x=0 and use it everywhere */
  for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)
      vy[Index(0,j,k)] = vy[Index(0, n2/4, k)];
}

/* second coord. derivs currently not needed */
/*
double drho_dxdx(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return y*y/pow(rho2, 1.5);
}
double drho_dxdy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return -x*y/pow(rho2, 1.5);
}
double drho_dydy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return x*x/pow(rho2, 1.5);
}

double dphi_dxdx(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return 2.0*x*y/( rho2*rho2 );
}
double dphi_dxdy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return (y*y - x*x)/( rho2*rho2 );
}
double dphi_dydy(void *aux, double x, double y, double z)
{
  double rho2 = x*x + y*y;
  return -2.0*x*y/( rho2*rho2 );
}
*/
/* end: Polar coordinates: */


/* start: PolarCE coordinates: */

/* Coord. trafos */
double x_ofPolarCE(void *aux, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );

  return rho*cos(phi);
}
double y_ofPolarCE(void *aux, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );

  return rho*sin(phi);
}

/* first coord. derivs */
/* NOTE: Y = PI*( 1.0 - cos( ((N-1.0)/(2.0*N)) * phi ) )             */
/* dY/dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi ) */
double dYPolarCE_dx(void *aux, double x, double y, double z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = atan(y/x);
  double dY_dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi );
  
  return dY_dphi * dphi_dx(aux, x, y, z);
}
double dYPolarCE_dy(void *aux, double x, double y, double z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = atan(y/x);
  double dY_dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi );
  
  return dY_dphi * dphi_dy(aux, x, y, z);
}
/* functions to treat cartesian derivs at singular points are currently not
   implemented */
/* second coord. derivs are currently not needed */

/* end: PolarCE coordinates: */
