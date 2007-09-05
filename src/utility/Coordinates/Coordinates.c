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
#define ArcTan(x)  atan(x)
#define Cosh(x)    cosh(x)
#define Sinh(x)    sinh(x)
#define Csch(x)    (1.0/sinh(x))
#define Sech(x)    (1.0/cosh(x))
#define Coth(x)    (1.0/tanh(x))

/* initialize the coord transforms */
int init_CoordTransform_And_Derivs(tGrid *grid)
{
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int dXd = Ind("dXdx");
  int dYd = Ind("dYdx");
  int dZd = Ind("dZdx");
  int ddXdd = Ind("ddXddxx");
  int ddYdd = Ind("ddYddxx");
  int ddZdd = Ind("ddZddxx");
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    printf("Coordinates: %s = %s\n", str, Gets(str));
    if( Getv(str, "Cartesian") )
    {
      printf("Coordinates: default Cartesian coordinates...\n");
    }
    else if( Getv(str, "Polar") )
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
    else if( Getv(str, "PolarCE") )
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
    else if( Getv(str, "SphericalDF") )
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
    else if( Getv(str, "compactSphericalDF") )
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
    else if( Getv(str, "tan_stretch") )
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

      box->ddX_dxdx[1][1][1] = ddxs_dxdx_tan_stretch;
      box->ddX_dxdx[1][1][2] = zero_of_xyz;
      box->ddX_dxdx[1][1][3] = zero_of_xyz;
      box->ddX_dxdx[1][2][2] = zero_of_xyz;
      box->ddX_dxdx[1][2][3] = zero_of_xyz;
      box->ddX_dxdx[1][3][3] = zero_of_xyz;

      box->ddX_dxdx[2][1][1] = zero_of_xyz;
      box->ddX_dxdx[2][1][2] = zero_of_xyz;
      box->ddX_dxdx[2][1][3] = zero_of_xyz;
      box->ddX_dxdx[2][2][2] = ddys_dydy_tan_stretch;
      box->ddX_dxdx[2][2][3] = zero_of_xyz;
      box->ddX_dxdx[2][3][3] = zero_of_xyz;

      box->ddX_dxdx[3][1][1] = zero_of_xyz;
      box->ddX_dxdx[3][1][2] = zero_of_xyz;
      box->ddX_dxdx[3][1][3] = zero_of_xyz;
      box->ddX_dxdx[3][2][2] = zero_of_xyz;
      box->ddX_dxdx[3][2][3] = zero_of_xyz;
      box->ddX_dxdx[3][3][3] = ddzs_dzdz_tan_stretch;
    }
    else if( Getv(str, "AnsorgNS0") )
    {
      printf("Coordinates: initializing AnsorgNS0 coordinates...\n");
      box->x_of_X[1] = x_of_AnsorgNS0;
      box->x_of_X[2] = y_of_AnsorgNS0;
      box->x_of_X[3] = z_of_AnsorgNS0;

      box->dX_dx[1][1] = dA_dx_AnsorgNS0;
      box->dX_dx[1][2] = dA_dy_AnsorgNS0;
      box->dX_dx[1][3] = dA_dz_AnsorgNS0;
      box->dX_dx[2][1] = dB_dx_AnsorgNS0;
      box->dX_dx[2][2] = dB_dy_AnsorgNS0;
      box->dX_dx[2][3] = dB_dz_AnsorgNS0;
      box->dX_dx[3][1] = dphi_dx_AnsorgNS0;
      box->dX_dx[3][2] = dphi_dy_AnsorgNS0;
      box->dX_dx[3][3] = dphi_dz_AnsorgNS0;
      
      //box->Sing_d_dx[2] = set_d_dy_at_rhoEQzero_AnsorgNS;
      //box->Sing_d_dx[3] = set_d_dz_at_rhoEQzero_AnsorgNS;
    }
    else if( Getv(str, "AnsorgNS1") )
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
      
      //box->Sing_d_dx[2] = set_d_dy_at_rhoEQzero_AnsorgNS;
      //box->Sing_d_dx[3] = set_d_dz_at_rhoEQzero_AnsorgNS;
    }
    else if( Getv(str, "AnsorgNS2") )
    {
      printf("Coordinates: initializing AnsorgNS2 coordinates...\n");
      box->x_of_X[1] = x_of_AnsorgNS2;
      box->x_of_X[2] = y_of_AnsorgNS2;
      box->x_of_X[3] = z_of_AnsorgNS2;

      box->dX_dx[1][1] = dA_dx_AnsorgNS2;
      box->dX_dx[1][2] = dA_dy_AnsorgNS2;
      box->dX_dx[1][3] = dA_dz_AnsorgNS2;
      box->dX_dx[2][1] = dB_dx_AnsorgNS2;
      box->dX_dx[2][2] = dB_dy_AnsorgNS2;
      box->dX_dx[2][3] = dB_dz_AnsorgNS2;
      box->dX_dx[3][1] = dphi_dx_AnsorgNS2;
      box->dX_dx[3][2] = dphi_dy_AnsorgNS2;
      box->dX_dx[3][3] = dphi_dz_AnsorgNS2;
      
      //box->Sing_d_dx[2] = set_d_dy_at_rhoEQzero_AnsorgNS;
      //box->Sing_d_dx[3] = set_d_dz_at_rhoEQzero_AnsorgNS;
    }
    else if( Getv(str, "AnsorgNS3") )
    {
      printf("Coordinates: initializing AnsorgNS3 coordinates...\n");
      box->x_of_X[1] = x_of_AnsorgNS3;
      box->x_of_X[2] = y_of_AnsorgNS3;
      box->x_of_X[3] = z_of_AnsorgNS3;

      box->dX_dx[1][1] = dA_dx_AnsorgNS3;
      box->dX_dx[1][2] = dA_dy_AnsorgNS3;
      box->dX_dx[1][3] = dA_dz_AnsorgNS3;
      box->dX_dx[2][1] = dB_dx_AnsorgNS3;
      box->dX_dx[2][2] = dB_dy_AnsorgNS3;
      box->dX_dx[2][3] = dB_dz_AnsorgNS3;
      box->dX_dx[3][1] = dphi_dx_AnsorgNS3;
      box->dX_dx[3][2] = dphi_dy_AnsorgNS3;
      box->dX_dx[3][3] = dphi_dz_AnsorgNS3;
      
      //box->Sing_d_dx[2] = set_d_dy_at_rhoEQzero_AnsorgNS;
      //box->Sing_d_dx[3] = set_d_dz_at_rhoEQzero_AnsorgNS;
    }
    else
      errorexit("Coordinates: unknown coordinates...");

    /* enable cartesian coordinates x,y,z */
    enablevar_inbox(box, var_x);
    enablevar_inbox(box, var_y);
    enablevar_inbox(box, var_z);

    /* enable storage for trafos? */
    if(Getv("CoordinateTransforms_stored", "yes"))
    {
      if( box->dX_dx[1][1] != NULL )
      {
        enablevar_inbox(box, dXd);
        enablevar_inbox(box, dYd);
        enablevar_inbox(box, dZd);
      }
      if( box->ddX_dxdx[1][1][1] != NULL )
      {
        enablevar_inbox(box, ddXdd);
        enablevar_inbox(box, ddYdd);
        enablevar_inbox(box, ddZdd);
      }
    }
  } /* end for b */

  /* compute cartesian coordinates x,y,z from X,Y,Z */
  compute_xyz_dXYZdxyz_ddXYZddxyz(grid);

  return 0;
}


/* compute coord values and coord trafos and store them in
   x,y,z, dXdx,dXdy,... , ddXddxx,ddXddxy,...               */
int compute_xyz_dXYZdxyz_ddXYZddxyz(tGrid *grid)
{
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int dXd = Ind("dXdx");
  int dYd = Ind("dYdx");
  int dZd = Ind("dZdx");
  int ddXdd = Ind("ddXddxx");
  int ddYdd = Ind("ddYddxx");
  int ddZdd = Ind("ddZddxx");
  int b, ind;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];

    if( box->x_of_X[1] != NULL )
    {
      forallpoints(box, ind)
      {
        int j,k,n;

        /* compute cartesian coordinates x,y,z from X,Y,Z */
        double X = box->v[var_X][ind];
        double Y = box->v[var_Y][ind];
        double Z = box->v[var_Z][ind];
        box->v[var_x][ind] = box->x_of_X[1]((void *) box, ind, X,Y,Z);
        box->v[var_y][ind] = box->x_of_X[2]((void *) box, ind, X,Y,Z);
        box->v[var_z][ind] = box->x_of_X[3]((void *) box, ind, X,Y,Z);
        
        /* compute dXdx, dXdy, ... */
        if( box->v[dXd] != NULL )
          for(j=1; j<=3; j++)
          {
            box->v[dXd + j-1][ind]=box->dX_dx[1][j]((void *) box, ind, X,Y,Z);
            box->v[dYd + j-1][ind]=box->dX_dx[2][j]((void *) box, ind, X,Y,Z);
            box->v[dZd + j-1][ind]=box->dX_dx[3][j]((void *) box, ind, X,Y,Z);
          }

        /* compute ddXddxx, ddXddxy, ... */
        if( box->v[ddXdd] != NULL )
        {
          n=0;
          for(j=1; j<=3; j++)
          for(k=j; k<=3; k++)
          {
            box->v[ddXdd + n][ind] 
              = box->ddX_dxdx[1][j][k]((void *) box, ind, X,Y,Z);
            box->v[ddYdd + n][ind] 
              = box->ddX_dxdx[2][j][k]((void *) box, ind, X,Y,Z);
            box->v[ddZdd + n][ind] 
              = box->ddX_dxdx[3][j][k]((void *) box, ind, X,Y,Z);
            n++;
          }
        } /* end if( box->v[ddXdd] != NULL ) */
      } /* end forallpoints */
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
  double *vx = (double *) v1;
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
  double *vy = (double *) v2;
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
double ddxs_dxdx(double xs)
{
  static double s=-1.0;
  double u, arg, dum;

  if(s < 0.0)
  {
    s=Getd("tan_stretch_s");
    printf(" Coordinates: dxs_dx: tan_stretch_s=%f\n", s);
  }
  
  arg = xs * s;
  if( fabs(arg)>=1.0 ) return 0.0;

  u = tan(arg * PI/2.0);
  
  dum = 1.0 / ( 1.0 + u*u );
  return -PI*s*u*dum*dum;
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

/* second coord. derivs are from an old sgrid version on loki,
   they may have some bugs!!! */
/* second coord. derivs */
double ddxs_dxdx_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return ddxs_dxdx(xs);
}
double ddys_dydy_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return ddxs_dxdx(ys);
}
double ddzs_dzdz_tan_stretch(void *aux, int ind, double xs, double ys, double zs)
{
  return ddxs_dxdx(zs);
}

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
  static int indsav=-1;
  static double Asav=-1, Bsav=-1, phisav=-1;
  static double xsav, ysav, zsav;
  double X,R;
  double Rsqr, Xsqr, Rsqr_p_Xsqr;
  double b, lep;
  double sigp_Bphi, sigp_1phi;
  double Ap;

  /* check if we have saved values */  
  if(ind==indsav && domain==domainsav) 
  {
    *x=xsav; *y=ysav; *z=zsav;
    return;
  }
  Asav=A;  Bsav=B;  phisav=phi;  domainsav=domain;

  /* set some pars */
  b = 1; // Getd("BNS_D")*0.5;
  if(domain==0 || domain==1)
  {
    lep = -1.5; // Getd("BNS_log_epsp");
    Ap = sinh(A*lep)/sinh(lep);
    if(domain==0)
    {
      double rootpower = 1;  // Getd("BNS_rootpower");
      Ap = 1.0 - pow(1.0-A, 1.0/rootpower);
    }
    sigp_Bphi = 1; //0.8 + 0.15*cos(B*2*PI) + 0.15*sin(phi); // 1; // change this!
    sigp_1phi = 1; //0.8 + 0.15 + 0.15*sin(phi); // 1; // change this!
  }
  if(domain==2 || domain==3)
  {
    lep = -1; // Getd("BNS_log_epsm");
    Ap = sinh(A*lep)/sinh(lep);
    if(domain==3)
    {
      double rootpower = 1;  // Getd("BNS_rootpower");
      Ap = 1.0 - pow(1.0-A, 1.0/rootpower);
    }
    sigp_Bphi = -1;// -1.2 + 0.1*cos(B*2*PI) + 0.1*sin(phi); // -1; // change this!
    sigp_1phi = -1;// -1.2 + 0.1 + 0.1*sin(phi); // -1; // change this!
  }

  /* compute coord trafo for each domain */
  if(domain==0) /* use Eq. (24) */
  {
    double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);

    /* use Eq. (24) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos((1.0-Ap)*ArgCp_1phi) + Ap*(1.0-B);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin((1.0-Ap)*ArgCp_1phi);
  }
  if(domain==1 || domain==2) /* use Eq. (22) */
  {
    double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);
    
    /* use Eq. (22) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
  }
  if(domain==3) /* use Eq. (23) */
  {
    double AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    double ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    double ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    double ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    double AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    double ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    double ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    double ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);

    /* use Eq. (23) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIh*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIh*Ap + (1.0-Ap)*ArgCp_1phi) + Ap*(1.0-B);
  }

/* Begin HACK1 */
//X=A;
//R=B;
/* End HACK1 */

  /* compute x,y,z */
  Rsqr = R*R;
  Xsqr = X*X;
  Rsqr_p_Xsqr = (Rsqr + Xsqr);
  if(Rsqr_p_Xsqr>0.0) /* if not at infinity */
  {
    double ooRsqr_p_Xsqr_sqr = 1.0/(Rsqr_p_Xsqr*Rsqr_p_Xsqr);

    *x = b*(ooRsqr_p_Xsqr_sqr + 1.0)*(Xsqr - Rsqr)*0.5;
    *y = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*cos(phi);
    *z = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X*sin(phi);
  }
  else
  {
    *x = *y = *z = 1.0e300;
  }
        
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
  static int indsav=-1;
  static double Asav=-1, Bsav=-1, phisav=-1;
  static double xsav, ysav, zsav;
  static double dAdxsav,   dAdysav,   dAdzsav,
                dBdxsav,   dBdysav,   dBdzsav,
                dphidxsav, dphidysav, dphidzsav;
  double X,R;
  double Rsqr, Xsqr, Rsqr_p_Xsqr;
  double dABphi_dXRphi[4][4]; /* dABphi_dXRphi[k][l] = dA^k/dX^l */
  double dXRphi_dxyz[4][4];   /* dXRphi_dxyz[l][m]   = dX^l/dx^m */
  int l;
  double b, lep;
  double sigp_Bphi, sigp_1phi;
  double dsigp_dphi_Bphi, dsigp_dphi_1phi, dsigp_dB_Bphi; /* dsigp_dB_1phi */
  double Ap;
  double dApdA;

  /* check if we have saved values */  
  if(ind==indsav && domain==domainsav) 
  {
    *x=xsav; *y=ysav; *z=zsav;
    *dAdx=dAdxsav;     *dAdy=dAdysav;     *dAdz=dAdzsav;
    *dBdx=dBdxsav;     *dBdy=dBdysav;     *dBdz=dBdzsav;
    *dphidx=dphidxsav; *dphidy=dphidysav; *dphidz=dphidzsav;
    return;
  }
  Asav=A;  Bsav=B;  phisav=phi;  domainsav=domain;

/* Begin HACK3a */
//A=0.3;
//B=0.35;
//phi=2;
/* End HACK3a */

  /* set some pars */
  b = 1; // Getd("BNS_D")*0.5;
  if(domain==0 || domain==1)
  {
    lep = -1.5; // Getd("BNS_log_epsp");
    Ap = sinh(A*lep)/sinh(lep);
    dApdA = lep*cosh(A*lep)/sinh(lep);
    if(domain==0)
    {
      double rootpower = 1;  // Getd("BNS_rootpower");
      Ap = 1.0 - pow(1.0-A, 1.0/rootpower);
      dApdA = pow(1.0-A, 1.0/rootpower - 1.0)/rootpower;
    }
    sigp_Bphi = 1;// 0.8 + 0.15*cos(B*2*PI) + 0.15*sin(phi); // 1; // change this!
    sigp_1phi = 1;//0.8 + 0.15 + 0.15*sin(phi); // 1; // change this!
    dsigp_dB_Bphi = 0;//-2*PI*0.15*sin(B*2*PI); // 0; // change this!
    /* dsigp_dB_1phi = 0; // change this! */
    dsigp_dphi_Bphi = 0;//0.15*cos(phi); // 0; // change this!
    dsigp_dphi_1phi = 0;//0.15*cos(phi); // 0; // change this!
  }
  if(domain==2 || domain==3)
  {
    lep = -1; // Getd("BNS_log_epsm");
    Ap = sinh(A*lep)/sinh(lep);
    dApdA = lep*cosh(A*lep)/sinh(lep);
    if(domain==3)
    {
      double rootpower = 1;  // Getd("BNS_rootpower");
      Ap = 1.0 - pow(1.0-A, 1.0/rootpower);
      dApdA = pow(1.0-A, 1.0/rootpower - 1.0)/rootpower;
    }
    sigp_Bphi = -1;//-1.2 + 0.1*cos(B*2*PI) + 0.1*sin(phi); // -1; // change this!
    sigp_1phi = -1;//-1.2 + 0.1 + 0.1*sin(phi); // -1; // change this!
    dsigp_dB_Bphi = 0;//-2*PI*0.1*sin(B*2*PI); // 0; // change this!
    /* dsigp_dB_1phi = 0; // change this! */
    dsigp_dphi_Bphi = 0;//0.1*cos(phi); // 0; // change this!
    dsigp_dphi_1phi = 0;//0.1*cos(phi); // 0; // change this!
  }

  /* compute coord trafo for each domain */
  if(domain==0) /* use Eq. (24) */
  {
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
                             0.25*fabs(dsigp_dphi_Bphi);
    double ArgdCp_dphi_Bphi =-ArgCp_Bphi +
                              2.0*Argsech(0.25*sigp_Bphi, 0.25*PI*B) +
                              Arg(dsigp_dphi_Bphi, 0.0);
    /* double AbsdCp_dB_1phi =(0.5/AbsCp_1phi)*Abssech(0.25*sigp_1phi, 0.25*PI)*
                           Abssech(0.25*sigp_1phi, 0.25*PI)*
                           0.25*sqrt(dsigp_dB_1phi*dsigp_dB_1phi + PI*PI);
       double ArgdCp_dB_1phi =-ArgCp_1phi+2.0*Argsech(0.25*sigp_1phi, 0.25*PI)+
                            Arg(dsigp_dB_1phi, PI);  */
    double AbsdCp_dphi_1phi =(0.5/AbsCp_1phi)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             0.25*fabs(dsigp_dphi_1phi);
    double ArgdCp_dphi_1phi =-ArgCp_1phi + 
                              2.0*Argsech(0.25*sigp_1phi, 0.25*PI) +
                              Arg(dsigp_dphi_1phi, 0.0);

    double dArgCp_dphi_1phi=(-(cosh(0.5*sigp_1phi))/
                              (sinh(0.5*sigp_1phi)*sinh(0.5*sigp_1phi)+ 1.0) )*
                              0.25*dsigp_dphi_1phi;

    double RedCp_dB_Bphi   = AbsdCp_dB_Bphi * cos(ArgdCp_dB_Bphi);
    double ImdCp_dB_Bphi   = AbsdCp_dB_Bphi * sin(ArgdCp_dB_Bphi);
    double RedCp_dphi_Bphi = AbsdCp_dphi_Bphi * cos(ArgdCp_dphi_Bphi);
    double ImdCp_dphi_Bphi = AbsdCp_dphi_Bphi * sin(ArgdCp_dphi_Bphi);
    /* double RedCp_dB_1phi   = AbsdCp_dB_1phi * cos(ArgdCp_dB_1phi);
       double ImdCp_dB_1phi   = AbsdCp_dB_1phi * sin(ArgdCp_dB_1phi); */
    double RedCp_dphi_1phi = AbsdCp_dphi_1phi * cos(ArgdCp_dphi_1phi);
    double ImdCp_dphi_1phi = AbsdCp_dphi_1phi * sin(ArgdCp_dphi_1phi);
    
    double dXdA = -(ReCp_Bphi - B*ReCp_1phi)*dApdA -
                  B*sin((1.0-Ap)*ArgCp_1phi)*(-ArgCp_1phi)*dApdA +(1-B)*dApdA;
    double dRdA = -(ImCp_Bphi - B*ImCp_1phi)*dApdA +
                  B*cos((1.0-Ap)*ArgCp_1phi)*(-ArgCp_1phi)*dApdA;
    double dXdB = (1.0-Ap)*(RedCp_dB_Bphi-ReCp_1phi) +
                  cos((1.0-Ap)*ArgCp_1phi) - Ap;
    double dRdB = (1.0-Ap)*(ImdCp_dB_Bphi-ImCp_1phi) +
                  sin((1.0-Ap)*ArgCp_1phi);
    double dXdphi=(1.0-Ap)*(RedCp_dphi_Bphi-B*RedCp_dphi_1phi) -
                  B*sin((1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    double dRdphi=(1.0-Ap)*(ImdCp_dphi_Bphi-B*ImdCp_dphi_1phi) +
                  B*cos((1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    double nenner, dAdX, dAdR, dAdphi, dBdX, dBdR, dBdphi;
    /* M = {{dXdA, dXdB, dXdphi}, {dRdA, dRdB, dRdphi}, {0,0,1}} 
       nenner = dRdB*dXdA - dRdA*dXdB
      In[4]:= Inverse[M]*nenner
      Out[4]= {{ dRdB, -dXdB,   dRdphi dXdB  - dRdB dXdphi},
               {-dRdA,  dXdA, -(dRdphi dXdA) + dRdA dXdphi},
               {0, 0, nenner}}    */
    if(A!=1.0)
    {
      nenner = dRdB*dXdA - dRdA*dXdB;
      dAdX   = dRdB/nenner;
      dAdR   =-dXdB/nenner;
      dAdphi = (dRdphi*dXdB - dRdB*dXdphi)/nenner;
      dBdX   =-dRdA/nenner;
      dBdR   = dXdA/nenner;
      dBdphi = (-(dRdphi*dXdA) + dRdA*dXdphi)/nenner;
    }
    else /* factor dRdB out of nenner and Zaehler */
    {
      double dXdBodRdB = (1 - ReCp_1phi + RedCp_dB_Bphi)/
                         (ArgCp_1phi - ImCp_1phi + ImdCp_dB_Bphi);
      dAdX   = 1.0/(dXdA - dRdA*dXdBodRdB);
      dAdR   =-dXdBodRdB/(dXdA - dRdA*dXdBodRdB);
      dAdphi = (dRdphi*dXdBodRdB - dXdphi)/(dXdA - dRdA*dXdBodRdB);
      dBdX   = dBdR = dBdphi = 0.0; /* allowed since  du/dB=0 at A=1 */
////printf("dXdBodRdB=%f/n", dXdBodRdB);
//printf("dAdX=%f  dAdR=%f  dAdphi=%f\n", dAdX,dAdR,dAdphi);
//if(A==1 && B==0)
//printf("dAdX=%f  dAdR=%f  dAdphi=%f\n", dAdX,dAdR,dAdphi);
    }
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

    /* use Eq. (24) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos((1.0-Ap)*ArgCp_1phi) + Ap*(1.0-B);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin((1.0-Ap)*ArgCp_1phi);
  }

  if(domain==1 || domain==2) /* use Eq. (22) */
  {
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
                             0.25*fabs(dsigp_dphi_Bphi);
    double ArgdCp_dphi_Bphi =-ArgCp_Bphi +
                              2.0*Argsech(0.25*sigp_Bphi, 0.25*PI*B) +
                              Arg(dsigp_dphi_Bphi, 0.0);
    /* double AbsdCp_dB_1phi =(0.5/AbsCp_1phi)*Abssech(0.25*sigp_1phi, 0.25*PI)*
                           Abssech(0.25*sigp_1phi, 0.25*PI)*
                           0.25*sqrt(dsigp_dB_1phi*dsigp_dB_1phi + PI*PI);
       double ArgdCp_dB_1phi =-ArgCp_1phi+2.0*Argsech(0.25*sigp_1phi, 0.25*PI)+
                            Arg(dsigp_dB_1phi, PI);  */
    double AbsdCp_dphi_1phi =(0.5/AbsCp_1phi)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             0.25*fabs(dsigp_dphi_1phi);
    double ArgdCp_dphi_1phi =-ArgCp_1phi +
                              2.0*Argsech(0.25*sigp_1phi, 0.25*PI) +
                              Arg(dsigp_dphi_1phi, 0.0);

    double dArgCp_dphi_1phi=(-(cosh(0.5*sigp_1phi))/
                              (sinh(0.5*sigp_1phi)*sinh(0.5*sigp_1phi)+ 1.0) )*
                              0.25*dsigp_dphi_1phi;

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

/* Begin HACK3c */
//printf("domain=%d\n", domain);
////printf("AbsCp_Bphi=%f ArgCp_Bphi=%f\n",AbsCp_Bphi,ArgCp_Bphi);
//printf("dXdA=%g dRdA=%g dXdB=%g dRdB=%g dXdphi=%g dRdphi=%g\n",
//        dXdA,dRdA,dXdB,dRdB,dXdphi,dRdphi);
//printf("dABphi_dXRphi[1][1]=%f\n", dABphi_dXRphi[1][1]);
//printf("dABphi_dXRphi[1][2]=%f\n", dABphi_dXRphi[1][2]);
//printf("dABphi_dXRphi[1][3]=%f\n", dABphi_dXRphi[1][3]);
//printf("dABphi_dXRphi[2][1]=%f\n", dABphi_dXRphi[2][1]);
//printf("dABphi_dXRphi[2][2]=%f\n", dABphi_dXRphi[2][2]);
//printf("dABphi_dXRphi[2][3]=%f\n", dABphi_dXRphi[2][3]);
//printf("dABphi_dXRphi[3][1]=%f\n", dABphi_dXRphi[3][1]);
//printf("dABphi_dXRphi[3][2]=%f\n", dABphi_dXRphi[3][2]);
//printf("dABphi_dXRphi[3][3]=%f\n", dABphi_dXRphi[3][3]);
//printf("dArgCp_dphi_1phi=%f\n",dArgCp_dphi_1phi);
////printf("RedCp_dB_Bphi=%f ImdCp_dB_Bphi=%f\n",
////        RedCp_dB_Bphi, ImdCp_dB_Bphi);
////printf("AbsdCp_dB_Bphi=%f ArgdCp_dB_Bphi=%f\n",
////        AbsdCp_dB_Bphi, ArgdCp_dB_Bphi);
////printf("AbsCp_Bphi=%f ArgCp_Bphi=%f\n",
////        AbsCp_Bphi, ArgCp_Bphi);
//printf("RedCp_dphi_Bphi=%f ImdCp_dphi_Bphi=%f\n",
//        RedCp_dphi_Bphi, ImdCp_dphi_Bphi);
//printf("AbsdCp_dphi_Bphi=%f ArgdCp_dphi_Bphi=%f\n",
//        AbsdCp_dphi_Bphi, ArgdCp_dphi_Bphi);
////printf("AbsCp_Bphi=%f ArgCp_Bphi=%f\n",
////        AbsCp_Bphi, ArgCp_Bphi);
//printf("RedCp_dphi_1phi=%f ImdCp_dphi_1phi=%f\n",
//        RedCp_dphi_1phi, ImdCp_dphi_1phi);
//printf("AbsdCp_dphi_1phi=%f ArgdCp_dphi_1phi=%f\n",
//        AbsdCp_dphi_1phi, ArgdCp_dphi_1phi);
////printf("AbsCp_1phi=%f ArgCp_1phi=%f\n",
////        AbsCp_1phi, ArgCp_1phi);
/* End HACK3c */

    /* use Eq. (22) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIq*Ap + (1.0-Ap)*ArgCp_1phi);
  }
  if(domain==3) /* use Eq. (23) */
  {
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
                             0.25*fabs(dsigp_dphi_Bphi);
    double ArgdCp_dphi_Bphi =-ArgCp_Bphi +
                              2.0*Argsech(0.25*sigp_Bphi, 0.25*PI*B) +
                              Arg(dsigp_dphi_Bphi, 0.0);
    /* double AbsdCp_dB_1phi =(0.5/AbsCp_1phi)*Abssech(0.25*sigp_1phi, 0.25*PI)*
                           Abssech(0.25*sigp_1phi, 0.25*PI)*
                           0.25*sqrt(dsigp_dB_1phi*dsigp_dB_1phi + PI*PI);
       double ArgdCp_dB_1phi =-ArgCp_1phi+2.0*Argsech(0.25*sigp_1phi, 0.25*PI)+
                            Arg(dsigp_dB_1phi, PI);  */
    double AbsdCp_dphi_1phi =(0.5/AbsCp_1phi)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             Abssech(0.25*sigp_1phi, 0.25*PI)*
                             0.25*fabs(dsigp_dphi_1phi);
    double ArgdCp_dphi_1phi =-ArgCp_1phi +
                              2.0*Argsech(0.25*sigp_1phi, 0.25*PI) +
                              Arg(dsigp_dphi_1phi, 0.0);

    double dArgCp_dphi_1phi=(-(cosh(0.5*sigp_1phi))/
                              (sinh(0.5*sigp_1phi)*sinh(0.5*sigp_1phi)+ 1.0) )*
                              0.25*dsigp_dphi_1phi;

    double RedCp_dB_Bphi   = AbsdCp_dB_Bphi * cos(ArgdCp_dB_Bphi);
    double ImdCp_dB_Bphi   = AbsdCp_dB_Bphi * sin(ArgdCp_dB_Bphi);
    double RedCp_dphi_Bphi = AbsdCp_dphi_Bphi * cos(ArgdCp_dphi_Bphi);
    double ImdCp_dphi_Bphi = AbsdCp_dphi_Bphi * sin(ArgdCp_dphi_Bphi);
    /* double RedCp_dB_1phi   = AbsdCp_dB_1phi * cos(ArgdCp_dB_1phi);
       double ImdCp_dB_1phi   = AbsdCp_dB_1phi * sin(ArgdCp_dB_1phi); */
    double RedCp_dphi_1phi = AbsdCp_dphi_1phi * cos(ArgdCp_dphi_1phi);
    double ImdCp_dphi_1phi = AbsdCp_dphi_1phi * sin(ArgdCp_dphi_1phi);
    
    double dXdA = -(ReCp_Bphi - B*ReCp_1phi)*dApdA -
                  B*sin(PIh*Ap + (1.0-Ap)*ArgCp_1phi)*(PIh-ArgCp_1phi)*dApdA;
    double dRdA = -(ImCp_Bphi - B*ImCp_1phi)*dApdA +
                  B*cos(PIh*Ap + (1.0-Ap)*ArgCp_1phi)*(PIh-ArgCp_1phi)*dApdA +
                  (1.0-B)*dApdA;
    double dXdB = (1.0-Ap)*(RedCp_dB_Bphi-ReCp_1phi) +
                  cos(PIh*Ap + (1.0-Ap)*ArgCp_1phi);
    double dRdB = (1.0-Ap)*(ImdCp_dB_Bphi-ImCp_1phi) +
                  sin(PIh*Ap + (1.0-Ap)*ArgCp_1phi) - Ap;
    double dXdphi=(1.0-Ap)*(RedCp_dphi_Bphi-B*RedCp_dphi_1phi) -
                  B*sin(PIh*Ap + (1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    double dRdphi=(1.0-Ap)*(ImdCp_dphi_Bphi-B*ImdCp_dphi_1phi) +
                  B*cos(PIh*Ap + (1.0-Ap)*ArgCp_1phi)*(1.0-Ap)*
                  dArgCp_dphi_1phi;
    double nenner, dAdX, dAdR, dAdphi, dBdX, dBdR, dBdphi;
    /* M = {{dXdA, dXdB, dXdphi}, {dRdA, dRdB, dRdphi}, {0,0,1}} 
       nenner = dRdB*dXdA - dRdA*dXdB
      In[4]:= Inverse[M]*nenner
      Out[4]= {{ dRdB, -dXdB,   dRdphi dXdB  - dRdB dXdphi},
               {-dRdA,  dXdA, -(dRdphi dXdA) + dRdA dXdphi},
               {0, 0, nenner}}    */
    if(A!=1.0)
    {
      nenner = dRdB*dXdA - dRdA*dXdB;
      dAdX   = dRdB/nenner;
      dAdR   =-dXdB/nenner;
      dAdphi = (dRdphi*dXdB - dRdB*dXdphi)/nenner;
      dBdX   =-dRdA/nenner;
      dBdR   = dXdA/nenner;
      dBdphi = (-(dRdphi*dXdA) + dRdA*dXdphi)/nenner;
    }
    else /* factor dRdB out of nenner and Zaehler */
    {
      double dXdBodRdB = (2*ArgCp_1phi - PI + 2*ReCp_1phi - 2*RedCp_dB_Bphi)/
                         (-2 + 2*ImCp_1phi - 2*ImdCp_dB_Bphi);
      dAdX   = 1.0/(dXdA - dRdA*dXdBodRdB);
      dAdR   =-dXdBodRdB/(dXdA - dRdA*dXdBodRdB);
      dAdphi = (dRdphi*dXdBodRdB - dXdphi)/(dXdA - dRdA*dXdBodRdB);
      dBdX   = dBdR = dBdphi = 0.0; /* allowed since  du/dB=0 at A=1 */
    }
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

    /* use Eq. (23) */
    X = (1.0-Ap)*(ReCp_Bphi - B*ReCp_1phi) + 
        B*cos(PIh*Ap + (1.0-Ap)*ArgCp_1phi);
    R = (1.0-Ap)*(ImCp_Bphi - B*ImCp_1phi) + 
        B*sin(PIh*Ap + (1.0-Ap)*ArgCp_1phi) + Ap*(1.0-B);
  }


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
  if(Rsqr_p_Xsqr>0.0) /* if not at infinity */
  {
    double ooRsqr_p_Xsqr_sqr = 1.0/(Rsqr_p_Xsqr*Rsqr_p_Xsqr);
    double ooRsqr_p_Xsqr_cube = ooRsqr_p_Xsqr_sqr/Rsqr_p_Xsqr;
    double rho = b*(ooRsqr_p_Xsqr_sqr - 1.0)*R*X;
    /* compute derivs */
    double dxdX = (-2*b*X*(-Rsqr + Xsqr))*ooRsqr_p_Xsqr_cube +
                   b*X*(1 + ooRsqr_p_Xsqr_sqr);
    double dxdR = (-2*b*R*(-Rsqr + Xsqr))*ooRsqr_p_Xsqr_cube -
                   b*R*(1 + ooRsqr_p_Xsqr_sqr);
    double drhodX=(-4*b*R*Xsqr)*ooRsqr_p_Xsqr_cube +
                   b*R*(-1 + ooRsqr_p_Xsqr_sqr);
    double drhodR=(-4*b*Rsqr*X)*ooRsqr_p_Xsqr_cube +
                   b*X*(-1 + ooRsqr_p_Xsqr_sqr);
    double det, dXdx, dXdrho, dRdx, dRdrho;
    
    if( !(dequal(R,0.0)&&dequal(X,1.0)) && !(dequal(R,1.0)&&dequal(X,0.0)) )
    {
      det = -(drhodX*dxdR) + drhodR*dxdX;
      /* M = {{dxdX, dxdR}, {drhodX, drhodR}}
         Inverse[M]*det = {{drhodR, -dxdR}, {-drhodX, dxdX}}  */
      dXdx   = drhodR/det;
      dXdrho =-dxdR/det;
      dRdx   =-drhodX/det;
      dRdrho = dxdX/det;
    }
    else /* we are at (X,R)=(1,0) or  (X,R)=(0,1) */
    {
      /* ad hoc regularization. CHECK THIS!!! */
      dXdx = dXdrho = dRdx = dRdrho = 0.0; // this must be wrong!
      printf("dABphi_dxyz_AnsorgNS: incorrect regularization at "
             "(X,R)=(1,0) or (X,R)=(0,1):\n"
             " dXdx = dXdrho = dRdx = dRdrho = 0.0\n");
    }
    dXRphi_dxyz[1][1]=dXdx;
    dXRphi_dxyz[1][2]=dXdrho*cos(phi);
    dXRphi_dxyz[1][3]=dXdrho*sin(phi);
    dXRphi_dxyz[2][1]=dRdx;
    dXRphi_dxyz[2][2]=dRdrho*cos(phi);
    dXRphi_dxyz[2][3]=dRdrho*sin(phi);
    dXRphi_dxyz[3][1]=0.0;
    if( !dequal(rho, 0.0) )
    {
      dXRphi_dxyz[3][2]=-sin(phi)/rho;
      dXRphi_dxyz[3][3]=cos(phi)/rho;
    }
    else /* regularize */
    {
      /* du/dx^m = du/dX dX/dx^m + du/dR dR/dx^m + du/dphi dphi/dx^m */
      /* thus if du/dphi=0 at rho=0 we can safely say dphi/dx^m = 0 */
      dXRphi_dxyz[3][2]=0.0;
      dXRphi_dxyz[3][3]=0.0;
    }
    /* compute x,y,z */
    *x = b*(ooRsqr_p_Xsqr_sqr + 1.0)*(Xsqr - Rsqr)*0.5;
    *y = rho*cos(phi);
    *z = rho*sin(phi);
  }
  else  /* if we are at infinity */
  {
    //printf("we are at infinity!!! Probably all dXRphi_dxyz are zero???\n");
    //printf("    Rsqr_p_Xsqr=%g:%d\n",Rsqr_p_Xsqr, Rsqr_p_Xsqr>0);
    /* Assume that all dXRphi_dxyz[l][m] are zero at infinity */
    // check this later!
    dXRphi_dxyz[1][1] = dXRphi_dxyz[1][2] = dXRphi_dxyz[1][3] = 0.0;
    dXRphi_dxyz[2][1] = dXRphi_dxyz[2][2] = dXRphi_dxyz[2][3] = 0.0;
    dXRphi_dxyz[3][1] = dXRphi_dxyz[3][2] = dXRphi_dxyz[3][3] = 0.0;
    dABphi_dXRphi[1][1] = dABphi_dXRphi[1][2] = dABphi_dXRphi[1][3] = 0.0;
    dABphi_dXRphi[2][1] = dABphi_dXRphi[2][2] = dABphi_dXRphi[2][3] = 0.0;
    dABphi_dXRphi[3][1] = dABphi_dXRphi[3][2] = dABphi_dXRphi[3][3] = 0.0;

    *x = *y = *z = 1.0e300;
  }
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
  /* *dphidx=*dphidy=*dphidz=0.0; */
  for(l=1; l<=3; l++)
  {
    *dAdx+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][1];
    *dAdy+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][2];
    *dAdz+=dABphi_dXRphi[1][l]*dXRphi_dxyz[l][3];

    *dBdx+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][1];
    *dBdy+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][2];
    *dBdz+=dABphi_dXRphi[2][l]*dXRphi_dxyz[l][3];

    /*  *dphidx+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][1];
        *dphidy+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][2];
        *dphidz+=dABphi_dXRphi[3][l]*dXRphi_dxyz[l][3]; */
  }
  *dphidx=dXRphi_dxyz[3][1]; /* dABphi_dXRphi[3][l]=1 if l=3 otherwise 0*/
  *dphidy=dXRphi_dxyz[3][2];
  *dphidz=dXRphi_dxyz[3][3];

  /* and save */
  xsav=*x; ysav=*y; zsav=*z;
  dAdxsav=*dAdx;     dAdysav=*dAdy;     dAdzsav=*dAdz;
  dBdxsav=*dBdx;     dBdysav=*dBdy;     dBdzsav=*dBdz;
  dphidxsav=*dphidx; dphidysav=*dphidy; dphidzsav=*dphidz;

//if( !finite(*dAdx) || !finite(*dAdy) || !finite(*dAdz) ||  
//    !finite(*dBdx) || !finite(*dBdy) || !finite(*dBdz) ||  
//    !finite(*dphidx) || !finite(*dphidy) || !finite(*dphidz) )
//{
////int k,l;
////for(k=1; k<=3; k++)
////for(l=1; l<=3; l++)
//{
////printf("dABphi_dXRphi[%d][%d]=%f ",k,l, dABphi_dXRphi[k][l]);
////printf("dXRphi_dxyz[%d][%d]=%f\n",k,l, dXRphi_dxyz[k][l]);
////printf("dXRphi_dxyz[%d][1]=%f rho=%g\n", l, dXRphi_dxyz[l][1], rho);
//printf("Rsqr_p_Xsqr=%g  %g %g %g  %g %g %g  %g %g %g\n",Rsqr_p_Xsqr, 
//*dAdx,*dAdy,*dAdz, *dBdx,*dBdy,*dBdz, *dphidx,*dphidy,*dphidz);
//}
//}
}

/* compute d^2(x,y,z)/(d(A,B,phi)d(A,B,phi)) */
void ddxyz_ddABphi_AnsorgNS(tBox *box, int ind, int domain,
                          double A, double B, double phi, 
                          double ddxyzddABphi[4][4][4])
{
/*
#define Pi PI
#define Ap(A) Ap_A
#define dApdA(A) dApdA_A
#define ddApddA(A) ddApddA_A
#define sigp1(phi) sigp_1phi
#define dsigpdB1(phi) dsigpdB_1phi
#define dsigpdphi1(phi) dsigpdphi_1phi
#define ddsigpddB1(phi) dsigpddB_1phi
#define ddsigpdBdphi1(phi) dsigpdBdphi_1phi
#define ddsigpddphi1(phi) dsigpddphi_1phi
#define sigp(B,phi) sigp_Bphi
#define dsigpdB(B,phi) dsigpdB_Bphi
#define dsigpdphi(B,phi) dsigpdphi_Bphi
#define ddsigpddB(B,phi) dsigpddB_Bphi
#define ddsigpdBdphi(B,phi) dsigpdBdphi_Bphi
#define ddsigpddphi(B,phi) dsigpddphi_Bphi
double b;                // we still have to set all those!
double Ap_A;
double dApdA_A;
double ddApddA_A;
double sigp_1phi;
double dsigpdB_1phi;
double dsigpdphi_1phi;
double dsigpddB_1phi;
double dsigpdBdphi_1phi;
double dsigpddphi_1phi;
double sigp_Bphi;
double dsigpdB_Bphi;
double dsigpdphi_Bphi;
double dsigpddB_Bphi;
double dsigpdBdphi_Bphi;
double dsigpddphi_Bphi;
// include math code 
#include "ddxyzddABphi_AnsorgNS0.c"
*/
}

/* compute d^2(A,B,phi)/(d(x,y,z)d(x,y,z)) */
void ddABphi_ddxyz_AnsorgNS(tBox *box, int ind, int domain,
                          double A, double B, double phi, 
                          double ddABphi_ddxyz[4][4][4])
{
  double ddxyz_ddABphi[4][4][4];
  double x,y,z;
  double dABphi_dxyz[4][4];
                                                      
  ddxyz_ddABphi_AnsorgNS(box, ind, domain, A,B,phi, ddxyz_ddABphi);
  dABphi_dxyz_AnsorgNS(box, ind, domain, A,B,phi, &x,&y,&z,
                  &dABphi_dxyz[1][1], &dABphi_dxyz[1][2], &dABphi_dxyz[1][3],
                  &dABphi_dxyz[2][1], &dABphi_dxyz[2][2], &dABphi_dxyz[2][3],
                  &dABphi_dxyz[3][1], &dABphi_dxyz[3][2], &dABphi_dxyz[3][3]);

  ddXdxdx_from_dXdx_ddxdXdX(ddABphi_ddxyz, dABphi_dxyz, ddxyz_ddABphi);
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
/* functions to treat cartesian derivs at singular points */
#if 0
void copy_var_along_x_AnsorgNS(tBox *box, int vi, int K0)
{
  double *u = box->v[vi];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;
  double *uline, *cline;

  /* take deriv d/dy at phi=pi/2 (k=K0) <=> y=0 and use it everywhere */
  uline = (double *) calloc(n1, sizeof(double));
  cline = (double *) calloc(n1, sizeof(double));

  /* loop over all points with B=0 and B=1 */  
  for(j=0; j<n2; j=j+n2-1)
  {
    int m;
    double a1=box->bbox[0];
    double b1=box->bbox[1];
        
    /* get coeffs of u along B=0 or B=1 */
    get_memline(u, uline, 1, j,K0, n1,n2,n3);
    matrix_times_vector(box->Mcoeffs1, uline, cline, n1);

    for(k=0; k<n3; k++)
    {
      if(k==K0) continue;
      
      for(i=0; i<n1; i++)
      {
        double x = box->v[Ind("x")][Index(i,j,k)];
        double sum = 0.0;
        double X = box->v[Ind("X")][Index(i,j,k)];
        double Y = box->v[Ind("Y")][Index(i,j,k)];
        double Z = 0.0;
        double x1, x2, x3, x4, xswap;
        x1 = box->v[Ind("x")][Index(n1-1,n2-1,K0)];
        x2 = box->v[Ind("x")][Index(0,n2-1,K0)];
        x3 = box->v[Ind("x")][Index(0,0,K0)];
        x4 = box->v[Ind("x")][Index(n1-1,0,K0)];
        if(x2<0) { xswap = x1;  x1=x2;  x2=xswap; }
        if(x3<0) { xswap = x3;  x3=x4;  x4=xswap; }
        if(x3<0 && x4>1e+299) x4=-x4;
        if( (x>x2 || x<x1) && (x>x4 || x<x3) )
        {
          errorexit("x out of range");
        }
        else
        {
          if(X==0.0) X+=1e-7;
          if(X==1.0) X-=1e-7;
          X_of_x_forgiven_YZ(box, &X, x, Y,Z);
          /* interpolate u to X */
          for(m=0; m<n1; m++)  sum += cline[m] * box->basis1(a1,b1, m,n1, X);
          u[Index(i,j,k)] = sum;
        }
      }
    }
  }
  free(uline);
  free(cline);
}
#endif
void set_d_dy_at_rhoEQzero_AnsorgNS(void *bo, void *va, 
                                    void *v1,void *v2,void *v3)
{
  tBox *box = (tBox *) bo;
  double *vy = (double *) v2;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;
  double rho_0 = box->v[Ind("y")][Index(0,0,0)];
  
  if(fabs(rho_0)>0.0) return;

  /* take deriv d/dy at phi=0 (k=0) <=> z=0 and use it everywhere */
  for(k=0; k<n3; k++)
    for(j=0; j<n2; j=j+n2-1)
      for(i=0; i<n1; i++)
        vy[Index(i,j,k)] = vy[Index(i,j,0)];
}
void set_d_dz_at_rhoEQzero_AnsorgNS(void *bo, void *va, 
                                    void *v1,void *v2,void *v3)
{
  tBox *box = (tBox *) bo;
  double *vz = (double *) v3;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;
  double rho_0 = box->v[Ind("y")][Index(0,0,0)];

  if(fabs(rho_0)>0.0) return;

  if(n3 % 4)
    errorexiti("set_d_dz_at_rhoEQzero_AnsorgNS: "
               "box[%d]->n3 has to be divisible by 4.", box->b);
  /* take deriv d/dy at phi=pi/2 (k=n3/4) <=> y=0 and use it everywhere */
  for(k=0; k<n3; k++)
    for(j=0; j<n2; j=j+n2-1)
      for(i=0; i<n1; i++)
        vz[Index(i,j,k)] = vz[Index(i,j,n3/4)];
}
void set_d_dz_at_rhoEQzero_AnsorgNS_new(void *bo, void *va, 
                                    void *v1,void *v2,void *v3)
{
  tBox *box = (tBox *) bo;
  double *u = (double *) v3;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;
  double rho_0 = box->v[Ind("y")][Index(0,0,0)];
  double *M, *uline, *cline;
  void (*get_coeffs)(double *,double *, int)=NULL;
  
  if(fabs(rho_0)>0.0) return;

  if(n3 % 4)
    errorexiti("set_d_dz_at_rhoEQzero_AnsorgNS: "
               "box[%d]->n3 has to be divisible by 4.", box->b);

  /* take deriv d/dy at phi=pi/2 (k=n3/4) <=> y=0 and use it everywhere */
  uline = (double *) calloc(n1, sizeof(double));
  cline = (double *) calloc(n1, sizeof(double));
  M = (double*) calloc(n1*n1, sizeof(double));

  get_spec_functionpointerTO_get_coeffs(box, 1, &get_coeffs);
  initMatrix_ForCoeffs(M, n1, get_coeffs);

  /* loop over all points with B=0 and B=1 */  
  for(j=0; j<n2; j=j+n2-1)
  {
    int m;
    double a1=box->bbox[0];
    double b1=box->bbox[1];
        
    /* get coeffs of u along B=0 or B=1 */
    get_memline(u, uline, 1, j,n3/4, n1,n2,n3);
    matrix_times_vector(M, uline, cline, n1);

    for(k=0; k<n3; k++)
    {
      if(k==n3/4) continue;
      
      for(i=0; i<n1; i++)
      {
        double x = box->v[Ind("x")][Index(i,j,k)];
//        double xmin0 = box->v[Ind("x")][Index(i,j,n3/4)];
        double sum = 0.0;
        double X = box->v[Ind("X")][Index(i,j,k)];
        double Y = box->v[Ind("Y")][Index(i,j,k)];
        double Z = 0.0;
        if(X==0.0) X+=1e-7;
        if(X==1.0) X-=1e-7;
        X_of_x_forgiven_YZ(box, &X, x, Y,Z);
        /* interpolate u to X */
        for(m=0; m<n1; m++)  sum += cline[m] * box->basis1(a1,b1, m,n1, X);
        u[Index(i,j,k)] = sum;
      }
    }
  }
  free(uline);
  free(cline);
  free(M);
}
/* second coord. derivs are currently not implemented */

/* end: _AnsorgNS coordinates: */

