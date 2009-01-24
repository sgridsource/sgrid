/* sgrid_DV_CircSchwSource2.h */
/* Wolfgang Tichy  5/2008 */

/***************************************************/
/* wrappers for functions we use in ScalarOnKerr.c */
/***************************************************/
void DV_set_parameters(double M, double q, double r0);
void DV_show_parameters(void);
double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS);
void DV_current_particle_xyzKS(double t_KS, double *x0, double *y0, double *z0);
