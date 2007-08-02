/* Coordinates.c */
/* Wolfgang Tichy 4/2005 */

#include "sgrid.h"
#include "Coordinates.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     (log((double) (x)))
#define Sin(x)     (sin((double) (x)))
#define Cos(x)     (cos((double) (x)))
#define Csc(x)     (1.0/sin((double) (x)))
#define Tan(x)     (tan((double) (x)))
#define Cot(x)     (1.0/tan((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))


/* initialize the coord transforms */
int init_CoordTransform_And_Derivs(tGrid *grid)
{
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int b, ind;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    if( Getv(str, "Polar") )
    {
      printf("Coordinates: initializing Polar coordinates...\n");
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
      box->ddX_dxdx[1][1][1] = drho_dxdx;
      box->ddX_dxdx[1][1][2] = drho_dxdy;
      box->ddX_dxdx[1][1][3] = zero_of_xyz;
      box->ddX_dxdx[1][2][2] = drho_dydy;
      box->ddX_dxdx[1][2][3] = zero_of_xyz;
      box->ddX_dxdx[1][3][3] = zero_of_xyz;

      box->ddX_dxdx[2][1][1] = dphi_dxdx;
      box->ddX_dxdx[2][1][2] = dphi_dxdy;
      box->ddX_dxdx[2][1][3] = zero_of_xyz;
      box->ddX_dxdx[2][2][2] = dphi_dydy;
      box->ddX_dxdx[2][2][3] = zero_of_xyz;
      box->ddX_dxdx[2][3][3] = zero_of_xyz;

      box->ddX_dxdx[3][1][1] = zero_of_xyz;
      box->ddX_dxdx[3][1][2] = zero_of_xyz;
      box->ddX_dxdx[3][1][3] = zero_of_xyz;
      box->ddX_dxdx[3][2][2] = zero_of_xyz;
      box->ddX_dxdx[3][2][3] = zero_of_xyz;
      box->ddX_dxdx[3][3][3] = zero_of_xyz;
      */
    }
    if( Getv(str, "PolarCE") )
    {
      printf("Coordinates: initializing PolarCE coordinates...\n");
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
    if( Getv(str, "SphericalDF") )
    {
      printf("Coordinates: initializing SphericalDF coordinates...\n");
      box->x_of_X[1] = x_ofSphericalDF;
      box->x_of_X[2] = y_ofSphericalDF;
      box->x_of_X[3] = z_ofSphericalDF;

      box->dX_dx[1][1] = dr_dx;
      box->dX_dx[1][2] = dr_dy;
      box->dX_dx[1][3] = dr_dz;
      box->dX_dx[2][1] = dthm_dx;
      box->dX_dx[2][2] = dthm_dy;
      box->dX_dx[2][3] = dthm_dz;
      box->dX_dx[3][1] = dphiSphericalDF_dx;
      box->dX_dx[3][2] = dphiSphericalDF_dy;
      box->dX_dx[3][3] = zero_of_xyz;

      box->ddX_dxdx[1][1][1] = ddr_SphericalDF_dxdx;
      box->ddX_dxdx[1][1][2] = ddr_SphericalDF_dxdy;
      box->ddX_dxdx[1][1][3] = ddr_SphericalDF_dxdz;
      box->ddX_dxdx[1][2][2] = ddr_SphericalDF_dydy;
      box->ddX_dxdx[1][2][3] = ddr_SphericalDF_dydz;
      box->ddX_dxdx[1][3][3] = ddr_SphericalDF_dzdz;

      box->ddX_dxdx[2][1][1] = ddthm_SphericalDF_dxdx;
      box->ddX_dxdx[2][1][2] = ddthm_SphericalDF_dxdy;
      box->ddX_dxdx[2][1][3] = ddthm_SphericalDF_dxdz;
      box->ddX_dxdx[2][2][2] = ddthm_SphericalDF_dydy;
      box->ddX_dxdx[2][2][3] = ddthm_SphericalDF_dydz;
      box->ddX_dxdx[2][3][3] = ddthm_SphericalDF_dzdz;

      box->ddX_dxdx[3][1][1] = ddphi_SphericalDF_dxdx;
      box->ddX_dxdx[3][1][2] = ddphi_SphericalDF_dxdy;
      box->ddX_dxdx[3][1][3] = ddphi_SphericalDF_dxdz;
      box->ddX_dxdx[3][2][2] = ddphi_SphericalDF_dydy;
      box->ddX_dxdx[3][2][3] = ddphi_SphericalDF_dydz;
      box->ddX_dxdx[3][3][3] = ddphi_SphericalDF_dzdz;
    }
    if( Getv(str, "compactSphericalDF") )
    {
      printf("Coordinates: initializing compactSphericalDF coordinates...\n");
      box->x_of_X[1] = x_ofcompactSphericalDF;
      box->x_of_X[2] = y_ofcompactSphericalDF;
      box->x_of_X[3] = z_ofcompactSphericalDF;

      box->dX_dx[1][1] = dxi_dx;
      box->dX_dx[1][2] = dxi_dy;
      box->dX_dx[1][3] = dxi_dz;
      box->dX_dx[2][1] = dthmcompactSphericalDF_dx;
      box->dX_dx[2][2] = dthmcompactSphericalDF_dy;
      box->dX_dx[2][3] = dthmcompactSphericalDF_dz;
      box->dX_dx[3][1] = dphicompactSphericalDF_dx;
      box->dX_dx[3][2] = dphicompactSphericalDF_dy;
      box->dX_dx[3][3] = zero_of_xyz;

      box->ddX_dxdx[1][1][1] = ddxi_compactSphericalDF_dxdx;
      box->ddX_dxdx[1][1][2] = ddxi_compactSphericalDF_dxdy;
      box->ddX_dxdx[1][1][3] = ddxi_compactSphericalDF_dxdz;
      box->ddX_dxdx[1][2][2] = ddxi_compactSphericalDF_dydy;
      box->ddX_dxdx[1][2][3] = ddxi_compactSphericalDF_dydz;
      box->ddX_dxdx[1][3][3] = ddxi_compactSphericalDF_dzdz;

      box->ddX_dxdx[2][1][1] = ddthm_compactSphericalDF_dxdx;
      box->ddX_dxdx[2][1][2] = ddthm_compactSphericalDF_dxdy;
      box->ddX_dxdx[2][1][3] = ddthm_compactSphericalDF_dxdz;
      box->ddX_dxdx[2][2][2] = ddthm_compactSphericalDF_dydy;
      box->ddX_dxdx[2][2][3] = ddthm_compactSphericalDF_dydz;
      box->ddX_dxdx[2][3][3] = ddthm_compactSphericalDF_dzdz;

      box->ddX_dxdx[3][1][1] = ddphi_compactSphericalDF_dxdx;
      box->ddX_dxdx[3][1][2] = ddphi_compactSphericalDF_dxdy;
      box->ddX_dxdx[3][1][3] = ddphi_compactSphericalDF_dxdz;
      box->ddX_dxdx[3][2][2] = ddphi_compactSphericalDF_dydy;
      box->ddX_dxdx[3][2][3] = ddphi_compactSphericalDF_dydz;
      box->ddX_dxdx[3][3][3] = ddphi_compactSphericalDF_dzdz;
    }
    if( Getv(str, "tan_stretch") )
    {
      printf("Coordinates: initializing tan_stretch coordinates...\n");
      box->x_of_X[1] = x_of_tan_stretch;
      box->x_of_X[2] = y_of_tan_stretch;
      box->x_of_X[3] = z_of_tan_stretch;

      box->dX_dx[1][1] = dxs_dx_tan_stretch;
      box->dX_dx[1][2] = zero_of_xyz;
      box->dX_dx[1][3] = zero_of_xyz;
      box->dX_dx[2][1] = zero_of_xyz;
      box->dX_dx[2][2] = dys_dy_tan_stretch;
      box->dX_dx[2][3] = zero_of_xyz;
      box->dX_dx[3][1] = zero_of_xyz;
      box->dX_dx[3][2] = zero_of_xyz;
      box->dX_dx[3][3] = dzs_dz_tan_stretch;
    }
    if( Getv(str, "AnsorgNS1") )
    {
      printf("Coordinates: initializing AnsorgNS1 coordinates...\n");
      box->x_of_X[1] = x_of_AnsorgNS1;
      box->x_of_X[2] = y_of_AnsorgNS1;
      box->x_of_X[3] = z_of_AnsorgNS1;

      box->dX_dx[1][1] = dA_dx_AnsorgNS1;
      box->dX_dx[1][2] = dA_dy_AnsorgNS1;
      box->dX_dx[1][3] = dA_dz_AnsorgNS1;
      box->dX_dx[2][1] = dB_dx_AnsorgNS1;
      box->dX_dx[2][2] = dB_dy_AnsorgNS1;
      box->dX_dx[2][3] = dB_dz_AnsorgNS1;
      box->dX_dx[3][1] = dphi_dx_AnsorgNS1;
      box->dX_dx[3][2] = dphi_dy_AnsorgNS1;
      box->dX_dx[3][3] = dphi_dz_AnsorgNS1;
    }

    /* compute cartesian coordinates x,y,z from X,Y,Z */
    enablevar_inbox(box, var_x);
    enablevar_inbox(box, var_y);
    enablevar_inbox(box, var_z);

    if( box->x_of_X[1] != NULL )
      forallpoints(box, ind)
      {
        double X = box->v[var_X][ind];
        double Y = box->v[var_Y][ind];
        double Z = box->v[var_Z][ind];
        box->v[var_x][ind] = box->x_of_X[1]((void *) box, ind, X,Y,Z);
        box->v[var_y][ind] = box->x_of_X[2]((void *) box, ind, X,Y,Z);
        box->v[var_z][ind] = box->x_of_X[3]((void *) box, ind, X,Y,Z);
      }
    else
      forallpoints(box, ind)
      {
        double X = box->v[var_X][ind];
        double Y = box->v[var_Y][ind];
        double Z = box->v[var_Z][ind];
        box->v[var_x][ind] = X;
        box->v[var_y][ind] = Y;
        box->v[var_z][ind] = Z;
      }

  }
  return 0;
}


/* Some trivial functions */
double zero_of_xyz(void *aux, int ind, double X, double Y, double Z)
{
  return 0.0;
}
double one_of_xyz(void *aux, int ind, double X, double Y, double Z)
{
  return 1.0;
}
double x_equals_X(void *aux, int ind, double X, double Y, double Z)
{
  return X;
}
double y_equals_Y(void *aux, int ind, double X, double Y, double Z)
{
  return Y;
}
double z_equals_Z(void *aux, int ind, double X, double Y, double Z)
{
  return Z;
}


/* ****************************************************************** */
/* start: Polar coordinates:                                          */

/* Coord. trafos */
double x_ofPolar(void *aux, int ind, double rho, double phi, double Z)
{
  return rho*cos(phi);
}
double y_ofPolar(void *aux, int ind, double rho, double phi, double Z)
{
  return rho*sin(phi);
}
double z_ofPolar(void *aux, int ind, double rho, double phi, double Z)
{
  return Z;
}

/* first coord. derivs */
double drho_dx(void *aux, int ind, double rho, double phi, double Z)
{
  return cos(phi);
}
double drho_dy(void *aux, int ind, double rho, double phi, double Z)
{
  return sin(phi);
}

double dphi_dx(void *aux, int ind, double rho, double phi, double Z)
{
  if(rho>0.0) return -sin(phi)/rho;
  else        return 0.0; /* result if we go along y=0 line */
}
double dphi_dy(void *aux, int ind, double rho, double phi, double Z)
{
  if(rho>0.0) return cos(phi)/rho;
  else        return 0.0; /* result if we go along x=0 line */
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
double drho_dxdx(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return y*y/pow(rho2, 1.5);
}
double drho_dxdy(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return -x*y/pow(rho2, 1.5);
}
double drho_dydy(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return x*x/pow(rho2, 1.5);
}

double dphi_dxdx(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return 2.0*x*y/( rho2*rho2 );
}
double dphi_dxdy(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return (y*y - x*x)/( rho2*rho2 );
}
double dphi_dydy(void *aux, int ind, double rho, double phi, double Z)
{
  double rho2 = x*x + y*y;
  return -2.0*x*y/( rho2*rho2 );
}
*/
/* end: Polar coordinates: */


/* ****************************************************************** */
/* start: PolarCE coordinates:                                        */

/* Coord. trafos */
double x_ofPolarCE(void *aux, int ind, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );

  return rho*cos(phi);
}
double y_ofPolarCE(void *aux, int ind, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );

  return rho*sin(phi);
}

/* first coord. derivs */
/* NOTE: Y = PI*( 1.0 - cos( ((N-1.0)/(2.0*N)) * phi ) )             */
/* dY/dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi )   */
double dYPolarCE_dx(void *aux, int ind, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi, dY_dphi;

  phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );
  dY_dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi );
  
  return dY_dphi * dphi_dx(aux,ind, rho, phi, Z);
}
double dYPolarCE_dy(void *aux, int ind, double rho, double Y, double Z)
{
  tBox *box = (tBox *) aux;
  int N = box->n2 - 1;
  double phi, dY_dphi;

  phi = ((2.0*N)/(N-1.0)) * acos( 1.0 - Y/PI );
  dY_dphi = PI*((N-1.0)/(2.0*N)) * sin( ((N-1.0)/(2.0*N)) * phi );
  
  return dY_dphi * dphi_dy(aux,ind, rho, phi, Z);
}
/* functions to treat cartesian derivs at singular points are currently not
   implemented */
/* second coord. derivs are currently not needed */

/* end: PolarCE coordinates: */


/* ****************************************************************** */
/* start: SphericalDF coordinates:                                      */

/* Coord. trafos */
double x_ofSphericalDF(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return r*cos(phi)*sin(theta);
}
double y_ofSphericalDF(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return r*sin(phi)*sin(theta);
}
double z_ofSphericalDF(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return r*cos(theta);
}

/* first coord. derivs */
double dr_dx(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return cos(phi)*sin(theta);
}
double dr_dy(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return sin(phi)*sin(theta);
}
double dr_dz(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return cos(theta);
}

double dthm_dx(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  if(r>0.0) return cos(theta)*cos(phi)/r;
  else      return 0.0; /* result if we go along y=0 line */
}
double dthm_dy(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  if(r>0.0) return cos(theta)*sin(phi)/r;
  else      return 0.0; /* result if we go along x=0 line */
}
double dthm_dz(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  if(r>0.0) return -sin(theta)/r;
  else      return 0.0; /* result if we go along z=0 line */
}

double dphiSphericalDF_dx(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return -sin(phi)/(r*sin(theta));
}
double dphiSphericalDF_dy(void *aux, int ind, double r, double thm, double phi)
{
  tBox *box = (tBox *) aux;
  int N = box->n2;
  double theta = thm + PI/((1+N%2)*N);

  return cos(phi)/(r*sin(theta));
}

/* second coord. derivs */
double ddr_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (Power(Cos(theta),2) + Power(Sin(phi),2)*Power(Sin(theta),2))/r;
}

double ddr_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(phi)*Sin(phi)*Power(Sin(theta),2))/r;
}

double ddr_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(phi)*Cos(theta)*Sin(theta))/r;
}

double ddr_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (Power(Cos(theta),2) + Power(Cos(phi),2)*Power(Sin(theta),2))/r;
}

double ddr_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(theta)*Sin(phi)*Sin(theta))/r;
}

double ddr_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return Power(Sin(theta),2)/r;
}

double ddthm_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-(Cos(2.*phi) - Power(Cos(phi),2)*Cos(2.*theta))*Cot(theta))/Power(r,2);
}

double ddthm_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (0.5*(-2. + Cos(2.*theta))*Cot(theta)*Sin(2.*phi))/Power(r,2);
}

double ddthm_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(phi)*Cos(2.*theta))/Power(r,2);
}

double ddthm_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (Cot(theta)*(Cos(2.*phi) + Cos(2.*theta)*Power(Sin(phi),2)))/Power(r,2);
}

double ddthm_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(2.*theta)*Sin(phi))/Power(r,2);
}

double ddthm_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return Sin(2.*theta)/Power(r,2);
}

double ddphi_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (Power(Csc(theta),2)*Sin(2.*phi))/Power(r,2);
}

double ddphi_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Cos(2.*phi)*Power(Csc(theta),2))/Power(r,2);
}

double ddphi_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi)
{
return 0.;
}

double ddphi_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);
return (-Power(Csc(theta),2)*Sin(2.*phi))/Power(r,2);
}

double ddphi_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi)
{
return 0.;
}

double ddphi_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi)
{
return 0.;
}

/* end: SphericalDF coordinates: */


/* *********************************************************** */
/* start: compactSphericalDF coordinates:                      */
/* xi = (2/PI) * arctan( r - r0)  <=>  r = r0 + tan(xi * PI/2) */
/* xi in [0,1] then r in [r0,infty] */
double r_of_xi(double xi)
{
  static double r0=-1.0;

//  if(xi==1.0) return  
  if(r0 < 0.0)
  {
    r0=Getd("compactSphericalDF_r0");
    printf(" Coordinates: r_of_xi: compactSphericalDF_r0=%f\n", r0);
  }
  return tan(xi * PI/2.0) + r0;
}
double dxi_dr(double xi)
{
  static double r0=-1.0;
  double r;

  if(xi==1.0) return 0.0;

  if(r0 < 0.0)
  {
    r0=Getd("compactSphericalDF_r0");
    printf(" Coordinates: dxi_dr: compactSphericalDF_r0=%f\n", r0);
  }
  
  r = tan(xi * PI/2.0) + r0;
  return (2.0/PI) / ( 1.0 + (r-r0)*(r-r0) );
}
/* xi = (2/PI) * arctan( r - r0)  <=>  r = r0 + tan(xi * PI/2) */
/* xi in [0,1] then r in [r0,infty] */

/* Coord. trafos */
double x_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi)
{
  return x_ofSphericalDF(aux,ind, r_of_xi(xi), thm, phi);
}
double y_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi)
{
  return y_ofSphericalDF(aux,ind, r_of_xi(xi), thm, phi);
}
double z_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi)
{
  return z_ofSphericalDF(aux,ind, r_of_xi(xi), thm, phi);
}

/* first coord. derivs */
/* xi = (2.0/PI) * arctan( r - r0) */
double dxi_dx(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return  dxi_dr(xi) * dr_dx(aux,ind, r_of_xi(xi), thm, phi);
}
double dxi_dy(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return  dxi_dr(xi) * dr_dy(aux,ind, r_of_xi(xi), thm, phi);
}
double dxi_dz(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return  dxi_dr(xi) * dr_dz(aux,ind, r_of_xi(xi), thm, phi);
}

double dthmcompactSphericalDF_dx(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return dthm_dx(aux,ind, r_of_xi(xi), thm, phi);
}
double dthmcompactSphericalDF_dy(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return dthm_dy(aux,ind, r_of_xi(xi), thm, phi);
}
double dthmcompactSphericalDF_dz(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return dthm_dz(aux,ind, r_of_xi(xi), thm, phi);
}

double dphicompactSphericalDF_dx(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return dphiSphericalDF_dx(aux,ind, r_of_xi(xi), thm, phi);
}
double dphicompactSphericalDF_dy(void *aux, int ind, double xi, double thm, double phi)
{
  if(xi==1.0) return 0.0;
  return dphiSphericalDF_dy(aux,ind, r_of_xi(xi), thm, phi);
}

/* second coord. derivs */
double ddxi_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (2.*Power(Cos(0.5*PI*xi),4)*(Power(Cos(theta),4)*Power(r0 + Tan(0.5*PI*xi),2) + Power(Cos(theta),2)*(1. + Power(r0,2) - 2.*r0*(r0 + Tan(0.5*PI*xi)) - Power(Cos(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) + 2.*Power(Sin(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2)) + Power(Sin(theta),2)*(2.*r0*Power(Cos(phi),2)*(r0 + Tan(0.5*PI*xi)) - 2.*Power(Cos(phi),4)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) + Power(Sin(phi),4)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) - Power(Sin(phi),2)*(-1. + Power(r0,2) + 2.*r0*Tan(0.5*PI*xi) + Power(Cos(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2)))))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddxi_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Power(Cos(0.5*PI*xi),2)*Sin(2.*phi)*Power(Sin(theta),2)*(-2. + Cos(PI*xi) - r0*Sin(PI*xi)))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddxi_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Cos(phi)*Power(Cos(0.5*PI*xi),2)*Sin(2.*theta)*(-2. + Cos(PI*xi) - r0*Sin(PI*xi)))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddxi_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (2.*Power(Cos(0.5*PI*xi),4)*(Power(Cos(theta),4)*Power(r0 + Tan(0.5*PI*xi),2) + Power(Cos(theta),2)*(1. + Power(r0,2) - 2.*r0*(r0 + Tan(0.5*PI*xi)) + 2.*Power(Cos(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) - Power(Sin(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2)) + Power(Sin(theta),2)*(2.*r0*Power(Sin(phi),2)*(r0 + Tan(0.5*PI*xi)) + Power(Cos(phi),4)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) - 2.*Power(Sin(phi),4)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2) - Power(Cos(phi),2)*(-1. + Power(r0,2) + 2.*r0*Tan(0.5*PI*xi) + Power(Sin(phi),2)*Power(Sin(theta),2)*Power(r0 + Tan(0.5*PI*xi),2)))))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddxi_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Power(Cos(0.5*PI*xi),2)*Sin(phi)*Sin(2.*theta)*(-2. + Cos(PI*xi) - r0*Sin(PI*xi)))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddxi_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (-2.*Power(Cos(0.5*PI*xi),2)*(Cos(2.*theta) + Power(Cos(theta),2)*(-Cos(PI*xi) + r0*Sin(PI*xi))))/(PI*(r0 + Tan(0.5*PI*xi)));
}

double ddthm_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return -(((Cos(2.*phi) - Power(Cos(phi),2)*Cos(2.*theta))*Power(Cos(0.5*PI*xi),2)*Cot(theta))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2));
}

double ddthm_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (0.5*(-2. + Cos(2.*theta))*Power(Cos(0.5*PI*xi),2)*Cot(theta)*Sin(2.*phi))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2);
}

double ddthm_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return -((Cos(phi)*Cos(2.*theta)*Power(Cos(0.5*PI*xi),2))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2));
}

double ddthm_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Power(Cos(0.5*PI*xi),2)*Cot(theta)*(Cos(2.*phi) + Cos(2.*theta)*Power(Sin(phi),2)))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2);
}

double ddthm_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return -((Cos(2.*theta)*Power(Cos(0.5*PI*xi),2)*Sin(phi))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2));
}

double ddthm_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Power(Cos(0.5*PI*xi),2)*Sin(2.*theta))/Power(r0*Cos(0.5*PI*xi) + Sin(0.5*PI*xi),2);
}

double ddphi_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (Power(Csc(theta),2)*Sin(2.*phi))/Power(r0 + Tan(0.5*PI*xi),2);
}

double ddphi_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return -((Cos(2.*phi)*Power(Csc(theta),2))/Power(r0 + Tan(0.5*PI*xi),2));
}

double ddphi_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi)
{
return 0.0;
}

double ddphi_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi)
{
tBox *box = (tBox *) aux;  int N = box->n2;  double theta = thm + PI/((1+N%2)*N);  static double r0=-1.0;
   if(r0 < 0.0)  {    r0=Getd("compactSphericalDF_r0");    printf(" Coordinates: compactSphericalDF_r0=%f\n", r0);  }
return (-2.*Cos(phi)*Power(Csc(theta),2)*Sin(phi))/Power(r0 + Tan(0.5*PI*xi),2);
}

double ddphi_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi)
{
return 0.0;
}

double ddphi_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi)
{
return 0.0;
}
/* end: compactSphericalDF coordinates: */


/* *********************************************** */
/* start: tan_stretch coordinates:                 */
/* x  = tan(xs * s * PI/2.0) * 2.0/(PI*s)          */
/* xs = arctan( x * s * PI/2.0) * 2.0/(PI*s)       */
/* s=1:  xs in [-1,1] then x in [-infty,infty]     */
/* s<<1: xs in [-1,1] then x ~ xs                  */
double x_of_xs(double xs)
{
  static double s=-1.0;

  if(s < 0.0)
  {
    s=Getd("tan_stretch_s");
    printf(" Coordinates: x_of_xs: tan_stretch_s=%f\n", s);
  }
  if(s==0) return xs;

  return tan(xs * s * PI/2.0) * 2.0/(PI*s);
}
double dxs_dx(double xs)
{
  static double s=-1.0;
  double u, arg;

  if(s < 0.0)
  {
    s=Getd("tan_stretch_s");
    printf(" Coordinates: dxs_dx: tan_stretch_s=%f\n", s);
  }
  
  arg = xs * s;
  if( fabs(arg)>=1.0 ) return 0.0;

  u = tan(arg * PI/2.0);
  
  return 1.0 / ( 1.0 + u*u );
}

/* Coord. trafos */
double x_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return x_of_xs(xs);
}
double y_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return x_of_xs(ys);
}
double z_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return x_of_xs(zs);
}

/* first coord. derivs */
/* xi = (2.0/PI) * arctan( r - r0) */
double dxs_dx_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return dxs_dx(xs);
}
double dys_dy_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return dxs_dx(ys);
}
double dzs_dz_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return dxs_dx(zs);
}

/* second coord. derivs are currently not needed */

/* end: _tan_stretch coordinates: */


/* ******************************************************************** */
/* start: AnsorgNS coordinates:						*/
/* 4 domains: 0=inside NS+, 1=outside NS+, 2=outside NS-, 3=inside NS-	*/
/* see gr-qc/0612081 v1							*/
/* A,B,phi are computational coords					*/
/* We need to do several coord. trafos:					*/
/* (A,B,phi) -> (X,R,phi)=C -> (x,rho,phi)=c -> (x,y,z)			*/
/* we may skip over (x,rho,phi)... */


/* compute (x,y,z) from (A,B,ph) and save result to speed up
   repeated calls */
void xyz_of_AnsorgNS(tBox *box, int ind, int domain,
                     double A, double B, double phi, 
                     double *x, double *y, double *z)
{
  static int domainsav=-1;
  static double Asav=-1, Bsav=-1, phisav=-1;
  static double xsav, ysav, zsav;
  double X,R;
  double Rsqr, Xsqr, ooRsqr_p_Xsqr_sqr;
  double b;

  /* check if we have saved values */  
  if(A==Asav && B==Bsav && phi==phisav && domain==domainsav) 
  {
    *x=xsav; *y=ysav; *z=zsav;
    return;
  }
  Asav=A;  Bsav=B;  phisav=phi;  domainsav=domain;
  b = 1; // Getd("BNS_D")*0.5;

  if(domain==0) yo();
  if(domain==1)
  {
    double lep = -1; // Getd("BNS_log_epsp");
    double Ap = sinh(A*lep)/sinh(lep);
    double sigp_Bphi = 1; // change this!
    double sigp_1phi = 1; // change this!
    double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
    
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
  }
  if(domain==2) yo();
  if(domain==3) yo();

/* Begin HACK1 */
//X=A;
//R=B;
/* End HACK1 */

  /* compute x,y,z */
  Rsqr = R*R;
  Xsqr = X*X;
  ooRsqr_p_Xsqr_sqr = 1.0/((Rsqr + Xsqr)*(Rsqr + Xsqr));
  *x = b*(ooRsqr_p_Xsqr_sqr + 1.0)*(Xsqr - Rsqr)*0.5;
  *y = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*cos(phi);
  *z = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*sin(phi);

/* Begin HACK2 */
//*x=X;
//*y=R;
//*z=phi;
/* End HACK2 */

  /* and save x,y,z */
  xsav=*x; ysav=*y; zsav=*z;
}

/* compute d(A,B,ph)/d(x,y,z) and save result to speed up
   repeated calls */
void dABphi_dxyz_AnsorgNS(tBox *box, int ind, int domain,
                          double A, double B, double phi, 
                          double *x, double *y, double *z,
                          double *dAdx,   double *dAdy,   double *dAdz,
                          double *dBdx,   double *dBdy,   double *dBdz,
                          double *dphidx, double *dphidy, double *dphidz)
{
  static int domainsav=-1;
  static double Asav=-1, Bsav=-1, phisav=-1;
  static double xsav, ysav, zsav;
  static double dAdxsav,   dAdysav,   dAdzsav,
                dBdxsav,   dBdysav,   dBdzsav,
                dphidxsav, dphidysav, dphidzsav;

  double X,R;
  double Rsqr, Xsqr, Rsqr_p_Xsqr;
  double b;
  double dABphi_dXRphi[4][4]; /* dABphi_dXRphi[k][l] = dA^k/dX^l */
  double dXRphi_dxyz[4][4];   /* dXRphi_dxyz[l][m]   = dX^l/dx^m */
  int l;

  /* check if we have saved values */  
  if(A==Asav && B==Bsav && phi==phisav && domain==domainsav) 
  {
    *x=xsav; *y=ysav; *z=zsav;
    *dAdx=dAdxsav;     *dAdy=dAdysav;     *dAdz=dAdzsav;
    *dBdx=dBdxsav;     *dBdy=dBdysav;     *dBdz=dBdzsav;
    *dphidx=dphidxsav; *dphidy=dphidysav; *dphidz=dphidzsav;
    return;
  }
  Asav=A;  Bsav=B;  phisav=phi;  domainsav=domain;
  b = 1; // Getd("BNS_D")*0.5;

/* Begin HACK3a */
//A=0.3;
//B=0.35;
//phi=0.4;
/* End HACK3a */

  if(domain==0) yo();
  if(domain==1) /* use Eq. (22) */
  {
    double lep = -1; // Getd("BNS_log_epsp");
    double Ap = sinh(A*lep)/sinh(lep);
    double dApdA = lep*cosh(A*lep)/sinh(lep);

    double sigp_Bphi = 1; // change this!
    double sigp_1phi = 1; // change this!
    double dsigp_dB_Bphi = 0; // change this!
    /* double dsigp_dB_1phi = 0; // change this! */
    double dsigp_dphi_Bphi = 0; // change this!
    double dsigp_dphi_1phi = 0; // change this!

    double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);

    double AbsdCp_dB_Bphi =(0.5/AbsCp_Bphi)*Abssech(0.25*sigp_Bphi, 0.25*PI*B)*
                           Abssech(0.25*sigp_Bphi, 0.25*PI*B)*
                           0.25*sqrt(dsigp_dB_Bphi*dsigp_dB_Bphi + PI*PI);
    double ArgdCp_dB_Bphi =-ArgCp_Bphi+2.0*Argsech(0.25*sigp_Bphi, 0.25*PI*B)+
                            Arg(dsigp_dB_Bphi, PI);
    double AbsdCp_dphi_Bphi =(0.5/AbsCp_Bphi)*
                             Abssech(0.25*sigp_Bphi, 0.25*PI*B)*
                             Abssech(0.25*sigp_Bphi, 0.25*PI*B)*
                             0.25*abs(dsigp_dphi_Bphi);
    double ArgdCp_dphi_Bphi =-ArgCp_Bphi+2.0*Argsech(0.25*sigp_Bphi, 0.25*PI*B);

    /* double AbsdCp_dB_1phi =(0.5/AbsCp_1phi)*Abssech(0.25*sigp_1phi, 0.25*PI*B)*
                           Abssech(0.25*sigp_1phi, 0.25*PI*B)*
                           0.25*sqrt(dsigp_dB_1phi*dsigp_dB_1phi + PI*PI);
       double ArgdCp_dB_1phi =-ArgCp_1phi+2.0*Argsech(0.25*sigp_1phi, 0.25*PI*B)+
                            Arg(dsigp_dB_1phi, PI);  */
    double AbsdCp_dphi_1phi =(0.5/AbsCp_1phi)*
                             Abssech(0.25*sigp_1phi, 0.25*PI*B)*
                             Abssech(0.25*sigp_1phi, 0.25*PI*B)*
                             0.25*abs(dsigp_dphi_1phi);
    double ArgdCp_dphi_1phi =-ArgCp_1phi+2.0*Argsech(0.25*sigp_1phi, 0.25*PI*B);

    double dArgCp_dphi_1phi=(-(sin(2.0*PI*B)*cosh(2.0*sigp_1phi))/
                              (sinh(2.0*sigp_1phi)*sinh(2.0*sigp_1phi)+
                               sin(2.0*PI*B)*sin(2.0*PI*B)) )*dsigp_dphi_1phi;

    double RedCp_dB_Bphi   = AbsdCp_dB_Bphi * cos(ArgdCp_dB_Bphi);
    double ImdCp_dB_Bphi   = AbsdCp_dB_Bphi * sin(ArgdCp_dB_Bphi);
    double RedCp_dphi_Bphi = AbsdCp_dphi_Bphi * cos(ArgdCp_dphi_Bphi);
    double ImdCp_dphi_Bphi = AbsdCp_dphi_Bphi * sin(ArgdCp_dphi_Bphi);
    /* double RedCp_dB_1phi   = AbsdCp_dB_1phi * cos(ArgdCp_dB_1phi);
       double ImdCp_dB_1phi   = AbsdCp_dB_1phi * sin(ArgdCp_dB_1phi); */
    double RedCp_dphi_1phi = AbsdCp_dphi_1phi * cos(ArgdCp_dphi_1phi);
    double ImdCp_dphi_1phi = AbsdCp_dphi_1phi * sin(ArgdCp_dphi_1phi);
    
    double dXdA = -(ReCp_Bphi - B*ReCp_1phi)*dApdA -
                  B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi)*(PIq-ArgCp_1phi)*dApdA;
    double dRdA = -(ImCp_Bphi - B*ImCp_1phi)*dApdA +
                  B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi)*(PIq-ArgCp_1phi)*dApdA;
    double dXdB = (1.0-Ap)*(RedCp_dB_Bphi-ReCp_1phi) +
                  cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    double dRdB = (1.0-Ap)*(ImdCp_dB_Bphi-ImCp_1phi) +
                  sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    double dXdphi=(1.0-Ap)*(RedCp_dphi_Bphi-B*RedCp_dphi_1phi) -
                  B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    double dRdphi=(1.0-Ap)*(ImdCp_dphi_Bphi-B*ImdCp_dphi_1phi) +
                  B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    /* M = {{dXdA, dXdB, dXdphi}, {dRdA, dRdB, dRdphi}, {0,0,1}} 
       nenner = dRdB*dXdA - dRdA*dXdB
      In[4]:= Inverse[M]*nenner
      Out[4]= {{ dRdB, -dXdB,   dRdphi dXdB  - dRdB dXdphi},
               {-dRdA,  dXdA, -(dRdphi dXdA) + dRdA dXdphi},
               {0, 0, nenner}}    */
    double nenner = dRdB*dXdA - dRdA*dXdB;
    double dAdX   = dRdB/nenner;
    double dAdR   =-dXdB/nenner;
    double dAdphi = (dRdphi*dXdB - dRdB*dXdphi)/nenner;
    double dBdX   =-dRdA/nenner;
    double dBdR   = dXdA/nenner;
    double dBdphi = (-(dRdphi*dXdA) + dRdA*dXdphi)/nenner;
    /* dphidX=0; dphidR=0; dphidphi=1; */
    dABphi_dXRphi[1][1] = dAdX;
    dABphi_dXRphi[1][2] = dAdR;
    dABphi_dXRphi[1][3] = dAdphi;
    dABphi_dXRphi[2][1] = dBdX;
    dABphi_dXRphi[2][2] = dBdR;
    dABphi_dXRphi[2][3] = dBdphi;
    dABphi_dXRphi[3][1] = dABphi_dXRphi[3][2] = 0.0;
    dABphi_dXRphi[3][3] = 1.0;

/* Begin HACK3b */
//printf("dXdA=%f dRdA=%f dXdB=%f dRdB=%f dXdphi=%f dRdphi=%f\n",
//        dXdA,dRdA,dXdB,dRdB,dXdphi,dRdphi);
//printf("RedCp_dB_Bphi=%f ImdCp_dB_Bphi=%f\n",        
//        RedCp_dB_Bphi, ImdCp_dB_Bphi);
//printf("AbsdCp_dB_Bphi=%f ArgdCp_dB_Bphi=%f\n",
//        AbsdCp_dB_Bphi, ArgdCp_dB_Bphi);
//printf("AbsCp_Bphi=%f ArgCp_Bphi=%f\n",
//        AbsCp_Bphi, ArgCp_Bphi);
/* End HACK3b */

    /* use Eq. (22) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
  }
  if(domain==2) yo();
  if(domain==3) yo();

/* Begin HACK1 */
//X=A;
//R=B;
//dABphi_dXRphi[1][1] = 1;
//dABphi_dXRphi[1][2] = dABphi_dXRphi[1][3] = 0;
//dABphi_dXRphi[2][1] = dABphi_dXRphi[2][3] = 0;
//dABphi_dXRphi[2][2] = 1;
//dABphi_dXRphi[3][1] = dABphi_dXRphi[3][2] = 0.0;
//dABphi_dXRphi[3][3] = 1.0;
/* End HACK1 */

  /* compute output vars from X,R,phi */
  Rsqr = R*R;
  Xsqr = X*X;
  Rsqr_p_Xsqr = (Rsqr + Xsqr);
  if(Rsqr_p_Xsqr>0) /* if not at infinity */
  {
    double ooRsqr_p_Xsqr_sqr = 1.0/(Rsqr_p_Xsqr*Rsqr_p_Xsqr);
    double ooRsqr_p_Xsqr_cube = ooRsqr_p_Xsqr_sqr/Rsqr_p_Xsqr;

    /* compute derivs */
    double dxDX = (-2*b*X*(-Rsqr + Xsqr))*ooRsqr_p_Xsqr_cube +
                   b*X*(1 + ooRsqr_p_Xsqr_sqr);
    double dxDR = (-2*b*R*(-Rsqr + Xsqr))*ooRsqr_p_Xsqr_cube -
                   b*R*(1 + ooRsqr_p_Xsqr_sqr);
    /* double dxDphi=0; */

    double dyDX=  (-4*b*R*Xsqr*cos(phi))*ooRsqr_p_Xsqr_cube +
                   b*R*(-1 + ooRsqr_p_Xsqr_sqr)*cos(phi);
    double dyDR=  (-4*b*Rsqr*X*cos(phi))*ooRsqr_p_Xsqr_cube +
                   b*X*(-1 + ooRsqr_p_Xsqr_sqr)*cos(phi);
    double dyDphi=-(b*R*X*(-1 + ooRsqr_p_Xsqr_sqr)*sin(phi));

    double dzDX=  (-4*b*R*Xsqr*sin(phi))*ooRsqr_p_Xsqr_cube +
                   b*R*(-1 + ooRsqr_p_Xsqr_sqr)*sin(phi);
    double dzDR=  (-4*b*Rsqr*X*sin(phi))*ooRsqr_p_Xsqr_cube +
                   b*X*(-1 + ooRsqr_p_Xsqr_sqr)*sin(phi);
    double dzDphi=b*R*X*(-1 + ooRsqr_p_Xsqr_sqr)*cos(phi);
    /* M = {{dxDX, dxDR, 0}, {dyDX, dyDR, dyDphi}, {dzDX, dzDR, dzDphi}}
       In[27]:= Inverse[M]
       In[28]:= Simplify[%]
       Out[28]= {{(dyDR dzDphi - dyDphi dzDR) /
                  (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                   dxDR dyDphi dzDX), (dxDR dzDphi) /
                   (-(dxDX dyDR dzDphi) + dxDR dyDX dzDphi + dxDX dyDphi dzDR -
                    dxDR dyDphi dzDX), (dxDR dyDphi) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX)}, {(-(dyDX dzDphi) + dyDphi dzDX) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX), (dxDX dzDphi) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX), (dxDX dyDphi) /
                   (-(dxDX dyDR dzDphi) + dxDR dyDX dzDphi + dxDX dyDphi dzDR -
                    dxDR dyDphi dzDX)}, {(dyDX dzDR - dyDR dzDX) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX), (-(dxDX dzDR) + dxDR dzDX) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX), (dxDX dyDR - dxDR dyDX) /
                   (dxDX dyDR dzDphi - dxDR dyDX dzDphi - dxDX dyDphi dzDR +
                    dxDR dyDphi dzDX)}}    */
    double nenner=dxDX*dyDR*dzDphi -dxDR*dyDX*dzDphi -dxDX*dyDphi*dzDR +
                       dxDR*dyDphi*dzDX;

    dXRphi_dxyz[1][1]=(dyDR*dzDphi - dyDphi*dzDR)/nenner;
    dXRphi_dxyz[1][2]=-dxDR*dzDphi/nenner;
    dXRphi_dxyz[1][3]=dxDR*dyDphi/nenner;
    dXRphi_dxyz[2][1]=(-(dyDX*dzDphi) + dyDphi*dzDX)/nenner;
    dXRphi_dxyz[2][2]=dxDX*dzDphi/nenner;
    dXRphi_dxyz[2][3]=-dxDX*dyDphi/nenner;
    dXRphi_dxyz[3][1]=(dyDX*dzDR - dyDR*dzDX)/nenner;
    dXRphi_dxyz[3][2]=(-(dxDX*dzDR) + dxDR*dzDX)/nenner;
    dXRphi_dxyz[3][3]=(dxDX*dyDR - dxDR*dyDX)/nenner;

    /* compute x,y,z */
    *x = b*(ooRsqr_p_Xsqr_sqr + 1.0)*(Xsqr - Rsqr)*0.5;
    *y = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*cos(phi);
    *z = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*sin(phi);
  }
  else
    printf("we are at infinity!!! Probably all dXRphi_dxyz are zero???\n");

/* Begin HACK2 */
//*x=X;
//*y=R;
//*z=phi;
//dXRphi_dxyz[1][1] = 1;
//dXRphi_dxyz[1][2] = dXRphi_dxyz[1][3] = 0;
//dXRphi_dxyz[2][1] = dXRphi_dxyz[2][3] = 0;
//dXRphi_dxyz[2][2] = 1;
//dXRphi_dxyz[3][1] = dXRphi_dxyz[3][2] = 0.0;
//dXRphi_dxyz[3][3] = 1.0;
/* End HACK2 */

  /* compute dA^k/dx^m */
  *dAdx = *dAdy = *dAdz = 0.0;
  *dBdx = *dBdy = *dBdz = 0.0;
  *dphidx=*dphidy=*dphidz=0.0;
  
  for(l=1; l<=3; l++) *dAdx+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][1];
  for(l=1; l<=3; l++) *dAdy+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][2];
  for(l=1; l<=3; l++) *dAdz+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][3];

  for(l=1; l<=3; l++) *dBdx+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][1];
  for(l=1; l<=3; l++) *dBdy+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][2];
  for(l=1; l<=3; l++) *dBdz+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][3];

  for(l=1; l<=3; l++) *dphidx+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][1];
  for(l=1; l<=3; l++) *dphidy+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][2];
  for(l=1; l<=3; l++) *dphidz+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][3];

  /* and save */
  xsav=*x; ysav=*y; zsav=*z;
  dAdxsav=*dAdx;     dAdysav=*dAdy;     dAdzsav=*dAdz;
  dBdxsav=*dBdx;     dBdysav=*dBdy;     dBdzsav=*dBdz;
  dphidxsav=*dphidx; dphidysav=*dphidy; dphidzsav=*dphidz;
}

/* Coord. trafos for domain 0 */
double x_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 0, A,B,phi, &x,&y,&z);
  return x;
}
double y_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 0, A,B,phi, &x,&y,&z);
  return y;
}
double z_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 0, A,B,phi, &x,&y,&z);
  return z;
}

/* first coord. derivs for domain 0 */
double dA_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdx;                                                                                                          
}
double dA_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdy;                                                                                                          
}
double dA_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdz;                                                                                                          
}
double dB_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdx;                                                                                                          
}
double dB_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdy;                                                                                                          
}
double dB_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdz;                                                                                                          
}
double dphi_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidx;                                                                                                          
}
double dphi_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidy;                                                                                                          
}
double dphi_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 0, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidz;                                                                                                          
}
/* second coord. derivs are currently not implemented */

/* Coord. trafos for domain 1 */
double x_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 1, A,B,phi, &x,&y,&z);
  return x;
}
double y_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 1, A,B,phi, &x,&y,&z);
  return y;
}
double z_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 1, A,B,phi, &x,&y,&z);
  return z;
}

/* first coord. derivs for domain 1 */
double dA_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdx;                                                                                                          
}
double dA_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdy;                                                                                                          
}
double dA_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdz;                                                                                                          
}
double dB_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdx;                                                                                                          
}
double dB_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdy;                                                                                                          
}
double dB_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdz;                                                                                                          
}
double dphi_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidx;                                                                                                          
}
double dphi_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidy;                                                                                                          
}
double dphi_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 1, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidz;                                                                                                          
}
/* second coord. derivs are currently not implemented */

/* Coord. trafos for domain 2 */
double x_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 2, A,B,phi, &x,&y,&z);
  return x;
}
double y_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 2, A,B,phi, &x,&y,&z);
  return y;
}
double z_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 2, A,B,phi, &x,&y,&z);
  return z;
}

/* first coord. derivs for domain 2 */
double dA_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdx;                                                                                                          
}
double dA_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdy;                                                                                                          
}
double dA_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdz;                                                                                                          
}
double dB_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdx;                                                                                                          
}
double dB_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdy;                                                                                                          
}
double dB_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdz;                                                                                                          
}
double dphi_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidx;                                                                                                          
}
double dphi_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidy;                                                                                                          
}
double dphi_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 2, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidz;                                                                                                          
}
/* second coord. derivs are currently not implemented */

/* Coord. trafos for domain 3 */
double x_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 3, A,B,phi, &x,&y,&z);
  return x;
}
double y_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 3, A,B,phi, &x,&y,&z);
  return y;
}
double z_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;

  xyz_of_AnsorgNS(box, ind, 3, A,B,phi, &x,&y,&z);
  return z;
}

/* first coord. derivs for domain 3 */
double dA_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdx;                                                                                                          
}
double dA_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdy;                                                                                                          
}
double dA_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dAdz;                                                                                                          
}
double dB_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdx;                                                                                                          
}
double dB_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdy;                                                                                                          
}
double dB_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dBdz;                                                                                                          
}
double dphi_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidx;                                                                                                          
}
double dphi_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidy;                                                                                                          
}
double dphi_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi)
{
  tBox *box = (tBox *) aux;
  double x,y,z;
  double dAdx,dAdy,dAdz, dBdx,dBdy,dBdz, dphidx,dphidy,dphidz;

  dABphi_dxyz_AnsorgNS(box, ind, 3, A,B,phi,  &x,&y,&z,
                       &dAdx,   &dAdy,   &dAdz, 
                       &dBdx,   &dBdy,   &dBdz,
                       &dphidx, &dphidy, &dphidz);
  return dphidz;                                                                                                          
}
/* second coord. derivs are currently not implemented */

/* end: _AnsorgNS coordinates: */

