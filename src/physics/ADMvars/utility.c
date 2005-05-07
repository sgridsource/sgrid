/* utility.c */
/* Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "ADMvars.h"



double detg(double g11, double g12, double g13, 
	    double g22, double g23, double g33)
{
  double gginv11, gginv12, gginv13;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  return g11*gginv11 + g12*gginv12 + g13*gginv13;
}



double invg(double g11, double g12, double g13, 
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
  double *gxx = Ptr(level, "gxx");
  double *gxy = Ptr(level, "gxy");
  double *gxz = Ptr(level, "gxz");
  double *gyy = Ptr(level, "gyy");
  double *gyz = Ptr(level, "gyz");
  double *gzz = Ptr(level, "gzz");

  double *Kxx = Ptr(level, "Kxx");
  double *Kxy = Ptr(level, "Kxy");
  double *Kxz = Ptr(level, "Kxz");
  double *Kyy = Ptr(level, "Kyy");
  double *Kyz = Ptr(level, "Kyz");
  double *Kzz = Ptr(level, "Kzz");

  double *psi = Ptr(level, "psi");
  double *K_initial = PtrEnable(level, "K_initial");

  double ixx, ixy, ixz, iyy, iyz, izz;
  int i,b;

  /* K_initial is zero upon enabling, return if that is what we want */
  if (!Getv("physics", "KerrSchild")) 
    return;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
      
    forallpoints(box, i)
    {
      invg(gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
	   &ixx, &ixy, &ixz, &iyy, &iyz, &izz);
    
      K_initial[i] = pow(psi[i],-4.0) * 
        (ixx*Kxx[i] + iyy*Kyy[i] + izz*Kzz[i] 
         + 2*(ixy*Kxy[i] + ixz*Kxz[i] + iyz*Kyz[i]));
    }
  }
}
