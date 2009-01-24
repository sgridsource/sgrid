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

/* get current xKS od particle */
double DV_current_xKS(double t_KS)
{
  return current_x(t_KS);
}
/* get current yKS od particle */
double DV_current_yKS(double t_KS)
{
  return current_y(t_KS);
}
/* get current zKS od particle */
double DV_current_zKS(double t_KS)
{
  return current_z(t_KS);
}
