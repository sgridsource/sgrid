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

//if(v==0) printf("v=0 ");

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


/* set u according to u-ua = w(r-vt)/r,
   where ua is the analytic soln at t=0  */
void set_boundary_radiative_analytic(tVarList *unew, tVarList *upre,
                                     double c, tVarList *ucur)
{
  tGrid *grid = unew->grid;
  int bi, pi, ijk;
  int i,j;
  double v;
  static tVarList *uinitial=NULL;

  if(ucur->time == 0.0)
  {
    printf("set_boundary_radiative_analytic: "
           "saving upre at t=%g\n", grid->time);
    /* duplicate upre into uinitial and turn on memory */
    vlfree(uinitial);
    uinitial = AddDuplicateEnable(upre, "_initial");
    for (j = 0; j < upre->n; j++)
    {
      forallboxes(grid,bi)
      {
        tBox *box  = grid->box[bi];
        double *up = box->v[upre->index[j]];
        double *ui = box->v[uinitial->index[j]];

        /* loop over all points and write into uinitial */
        forallpoints(box, ijk)  ui[ijk] = up[ijk];
      }
    }
  }

  /* for all variables */
  for (j = 0; j < unew->n; j++)
  {
    char str[1000];
    int slen;

    i   = unew->index[j];
    v   = VarPropSpeed(i);

    /* keep this variable constant, just copy it */
    snprintf(str, 999, "%s", VarName(i));  slen=strlen(str);
    if(Getv("boundary_radconstant", str)) v = 0;

    /* remove endings like "_r" from rk4 from var name in str */
    if(slen>=3)
      if(str[slen-2]=='_') str[slen-2]=0;
    if(Getv("boundary_radconstant", str)) v = 0;

    /* do it */
    /* u-ua = w(r-vt)/r  =>  w = r (u-ua)  =>  w' = (u-ua) + r (u'-ua') 
       udot-uadot = -v w'/r 
                  = -v (u-ua)/r -v (u'-ua')                              */
    forallboxes(grid,bi)
    {
      tBox *box=grid->box[bi];
      double *xp = box->v[Ind("x")];
      double *yp = box->v[Ind("y")];
      double *zp = box->v[Ind("z")];
      double *var = box->v[ucur->index[j]];
      double *varnew = box->v[unew->index[j]];
      double *varpre = box->v[upre->index[j]];
      double *dvar_dx = box->v[Ind("temp1")];
      double *dvar_dy = box->v[Ind("temp2")];
      double *dvar_dz = box->v[Ind("temp3")];
      double *varinit = box->v[uinitial->index[j]];
      double *dvarinit_dx = box->v[Ind("ADMvars_dgxxx")]; /* we again abuse */
      double *dvarinit_dy = box->v[Ind("ADMvars_dgxxy")]; /* ADMvars_dg     */
      double *dvarinit_dz = box->v[Ind("ADMvars_dgxxz")];
      double r, vx, vy, vz, x, y, z, rhs;

      cart_partials(box, var, dvar_dx, dvar_dy, dvar_dz); 
      cart_partials(box, varinit, dvarinit_dx, dvarinit_dy, dvarinit_dz); 

      forPointList_inbox(radiativeBoundaryPointList, box, pi , ijk)
      {
        x = xp[ijk];
        y = yp[ijk];
        z = zp[ijk];
        r = sqrt(x*x + y*y + z*z);
        
        vx = v * x/r;
        vy = v * y/r;
        vz = v * z/r;

        /* u-ua = w(r-vt)/r  =>  w = r (u-ua)  =>  w' = (u-ua) + r (u'-ua') 
           udot-uadot = -v w'/r 
                      = -v (u-ua)/r -v (u'-ua')                            */
        rhs = -v*(var[ijk] - varinit[ijk])/r; 

        rhs -= vx*(dvar_dx[ijk]-dvarinit_dx[ijk]) + 
               vy*(dvar_dy[ijk]-dvarinit_dy[ijk]) + 
               vz*(dvar_dz[ijk]-dvarinit_dz[ijk]);

        if (c != 0.0) 
          varnew[ijk] = varpre[ijk] + c*rhs;
        else
          varnew[ijk] = rhs;
      }
    }
  }
}
