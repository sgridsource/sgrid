/* radiative.c */
/* Bernd Bruegmann 10/02, Wolfgang Tichy 12/2003 */

/* Radiative boundary condition without corotation */

#include "sgrid.h"
#include "boundary.h"


/* set radiative boundary for one variable */
void set_boundary_radiative(tGrid *grid, 
	                    int unew, int upre, double c, int ucur,
			    double var0, double v) 
{
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *xp = box->v[Ind("x")];
    double *yp = box->v[Ind("y")];
    double *zp = box->v[Ind("z")];
    double *var = box->v[ucur];
    double *varnew = box->v[unew];
    double *varpre = box->v[upre];
    double *dvar_dx = box->v[Ind("temp1")];
    double *dvar_dy = box->v[Ind("temp2")];
    double *dvar_dz = box->v[Ind("temp3")];
    double r, vx, vy, vz, x, y, z, rhs;

    cart_partials(box, var, dvar_dx, dvar_dy, dvar_dz); 

    forPointList_inbox(radiativeBoundaryPointList, box, pi , ijk)
    {
      x = xp[ijk];
      y = yp[ijk];
      z = zp[ijk];
      r = sqrt(x*x + y*y + z*z);
      
      vx = v * x/r;
      vy = v * y/r;
      vz = v * z/r;

      rhs = -v*(var[ijk] - var0)/r; 

      rhs -= vx*dvar_dx[ijk] + vy*dvar_dy[ijk] + vz*dvar_dz[ijk];

      if (c != 0.0) 
        varnew[ijk] = varpre[ijk] + c*rhs;
      else
        varnew[ijk] = rhs;
    }
  }
}

