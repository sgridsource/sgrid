/* utility.c */
/* Wolfgang Tichy, 2005 */

#include "sgrid.h"
#include "ADMvars.h"



double ADMvars_detg(double g11, double g12, double g13, 
	    double g22, double g23, double g33)
{
  double gginv11, gginv12, gginv13;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  return g11*gginv11 + g12*gginv12 + g13*gginv13;
}



double ADMvars_invg(double g11, double g12, double g13, 
	    double g22, double g23, double g33,
	    double *i11, double *i12, double *i13, 
	    double *i22, double *i23, double *i33)
{
  double detg, gginv11, gginv12, gginv13, gginv22, gginv23, gginv33;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  gginv22 = g11*g33 - g13*g13;
  gginv23 = g12*g13 - g11*g23;
  gginv33 = g11*g22 - g12*g12;
  detg = g11*gginv11 + g12*gginv12 + g13*gginv13;
  *i11 = gginv11/detg;
  *i12 = gginv12/detg;
  *i13 = gginv13/detg;
  *i22 = gginv22/detg;
  *i23 = gginv23/detg;
  *i33 = gginv33/detg;
  return detg;
}




/* compute and save initial value of the trace of K */
int set_K_initial(tGrid *grid)
{
  double ixx, ixy, ixz, iyy, iyz, izz;
  int i,b;

  enablevar(grid, Ind("K_initial"));
  
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *gxx = box->v[Ind("gxx")];
    double *gxy = box->v[Ind("gxy")];
    double *gxz = box->v[Ind("gxz")];
    double *gyy = box->v[Ind("gyy")];
    double *gyz = box->v[Ind("gyz")];
    double *gzz = box->v[Ind("gzz")];
  
    double *Kxx = box->v[Ind("Kxx")];
    double *Kxy = box->v[Ind("Kxy")];
    double *Kxz = box->v[Ind("Kxz")];
    double *Kyy = box->v[Ind("Kyy")];
    double *Kyz = box->v[Ind("Kyz")];
    double *Kzz = box->v[Ind("Kzz")];
  
    double *psi = box->v[Ind("psi")];
    double *K_initial = box->v[Ind("K_initial")];
      
    forallpoints(box, i)
    {
      ADMvars_invg(gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
        	   &ixx, &ixy, &ixz, &iyy, &iyz, &izz);
    
      K_initial[i] = pow(psi[i],-4.0) * 
        (ixx*Kxx[i] + iyy*Kyy[i] + izz*Kzz[i] 
         + 2*(ixy*Kxy[i] + ixz*Kxz[i] + iyz*Kyz[i]));
    }
  }
  return 0;
}


/* undo conformal transformation:
   g = psi^4 g
   psi = 1, dpsi = 0, ddpsi = 0 */
int ADMvars_undo_conformal_split(tGrid *grid)
{
  tBox *box;
  double psi4;
  int ig = Ind("gxx");
  int ip = Ind("psi");
  int idpop = Ind("dpsiopsix");
  int iddpop = Ind("ddpsiopsixx");
  int b, i, n;

  printf("ADMvars_undo_conformal_split: g -> psi^4 g,  psi -> 1 \n");

  forallgridpoints(grid,box, b, i)
  {
    double *psi = box->v[ip];

    /* set psi4 and then set psi=1 */
    psi4 = pow(psi[i], 4.0);
    psi[i] = 1;

    /* rescale metric */
    for(n = 0; n < 6; n++) box->v[ig+n][i] *= psi4;

    if (box->v[idpop])
    {
      /* set all derivs of psi to zero */
      for(n = 0; n < 3; n++) box->v[idpop+n][i] = 0;
      for(n = 0; n < 6; n++) box->v[iddpop+n][i] = 0;
    }
  }

  return 0;
}
