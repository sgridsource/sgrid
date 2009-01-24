/* sgrid_DV_CircSchwSource.h */
/* Wolfgang Tichy  5/2008 */

/***************************************************/
/* wrappers for functions we use in ScalarOnKerr.c */
/***************************************************/
void DV_set_parameters(double M, double q, double r0);
double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS);
double DV_current_xKS(double t_KS);
double DV_current_yKS(double t_KS);
double DV_current_zKS(double t_KS);
