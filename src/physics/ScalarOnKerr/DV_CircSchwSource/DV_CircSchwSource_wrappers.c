/* DV_CircSchwSource_wrappers.c */
/* Wolfgang Tichy 1/2009 */

#include "sgrid.h"
#include "DV_CircSchwSource.h"
#include "Constants.h"
#include "Source.h"


/***************************************************/
/* wrappers for functions we use in ScalarOnKerr.c */
/***************************************************/
/* initialize internal pars in Ian's code */
void DV_set_parameters(double M, double q, double r0)
{
 set_parameters(M, q, r0, Getd("DV_CircSchwSource_DVWindow_n"),
                          Getd("DV_CircSchwSource_DVWindow_width"));
}

/* get Source at (tKS, xKS,yKS,zKS) */
double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS)
{
  return SourceInKerrSchild(tKS, xKS,yKS,zKS);
}

/* get current KS coords of particle position */
void DV_current_particle_xyzKS(double t_KS, double *x0, double *y0, double *z0)
{
  *x0 = current_x(t_KS);
  *y0 = current_y(t_KS);
  *z0 = current_z(t_KS);
}
