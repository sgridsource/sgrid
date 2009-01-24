/* DV_CircSchwSource2_wrappers.cc */
/* Wolfgang Tichy 1/2009 */

#include <math.h>
#include "WolfSource.hh"
//#include "sgrid.h"
//#include "DV_CircSchwSource2.h"
/* do not include "sgrid.h", just declare the funcs we are using here */
extern "C" double Getd(char *name);

/* tell the C++ compiler to create functions that are callabe from C */
extern "C" void DV_set_parameters(double M, double q, double r0);
extern "C" void DV_show_parameters(void);
extern "C" double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS);
extern "C" void DV_current_particle_xyzKS(double t_KS, double *x0, double *y0, double *z0);

#define PI  3.14159265358979323846264338327950

/* global structure var DV_params that contains constants */
constants DV_params[1];


/***************************************************/
/* wrappers for functions we use in ScalarOnKerr.c */
/***************************************************/
/* initialize internal pars in Ian's code */
void DV_set_parameters(double M, double q, double r0)
{
  double q1=1.2, s1=1.9, r1=0, r2=0;
  double q3=0,   s3=0,   r3=0, r4=0;
  int smooth=0;

  set_orbit(DV_params, r0, M, q);
  // set_NoWindow(DV_params); // sets the window function to 1 everyhere
  set_OrigWindow(DV_params, Getd("DV_CircSchwSource_DVWindow_n"),
                 Getd("DV_CircSchwSource_DVWindow_width"), smooth);
  //set_WolfWindow(DV_params, q1,s1,r1,r2, q3,s3,r3,r4, smooth);
}

void DV_show_parameters(void)
{
  show_parameters(DV_params);
}
  
/* get Source at (tKS, xKS,yKS,zKS) */
double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS)
{
  return SourceInKerrSchild(tKS, xKS,yKS,zKS, DV_params);
}

/* get current KS coords of particle position */
void DV_current_particle_xyzKS(double tKS, double *x0, double *y0, double *z0)
{
  double M     = DV_params->M;
  double r0    = DV_params->R;
  double Omega = DV_params->Omega;
  double tSchw = tKS - 2.0*M*log(r0/(2.0*M)-1.0);
  double phi0  = Omega*tSchw;
  //double theta0;
  int n;

  /* particle position??? */
  //theta0 = PIh;
  n = phi0/(2.0*PI);
  phi0 -= n*2*PI; /* stay in [0,2PI] */
  *x0 = r0*cos(phi0);
  *y0 = r0*sin(phi0);
  *z0 = 0.0;
}
