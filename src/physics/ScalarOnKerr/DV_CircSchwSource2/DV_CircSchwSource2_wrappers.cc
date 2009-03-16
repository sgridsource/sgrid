/* DV_CircSchwSource2_wrappers.cc */
/* Wolfgang Tichy 1/2009 */

#include <math.h>
#include "WolfSource.hh"
//#include "sgrid.h"
//#include "DV_CircSchwSource2.h"
/* do not include "sgrid.h", just declare the funcs we are using here */
extern "C" double Getd(char *name);
extern "C" int    Geti(char *name);
extern "C" int    Getv(char *name, char *value);

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
  int I_n = Geti("DV_CircSchwSource_DVWindow_n");
  double I_width = Getd("DV_CircSchwSource_DVWindow_width");
  int smooth = Geti("DV_CircSchwSource_Window_smooth");
  double q1 = Getd("DV_CircSchwSource_Window_q1");
  double s1 = Getd("DV_CircSchwSource_Window_s1");
  double r1 = Getd("DV_CircSchwSource_Window_r1");
  double r2 = Getd("DV_CircSchwSource_Window_r2");
  double q3 = Getd("DV_CircSchwSource_Window_q3");
  double s3 = Getd("DV_CircSchwSource_Window_s3");
  double r3 = Getd("DV_CircSchwSource_Window_r3");
  double r4 = Getd("DV_CircSchwSource_Window_r4");

  set_orbit(DV_params, r0, M, q);
  
  if(Getv("DV_CircSchwSource_useWindow", "no") ||
     Getv("DV_CircSchwSource_Window_type", "no"))
    set_NoWindow(DV_params); // sets the window function to 1 everyhere
  else if(Getv("DV_CircSchwSource_Window_type", "orig"))
    set_OrigWindow(DV_params, I_n, I_width, smooth);
  else // if(Getv("DV_CircSchwSource_DVWindow_type", "wolf"))
    set_WolfWindow(DV_params, q1,s1,r1,r2, q3,s3,r3,r4, smooth);
}

void DV_show_parameters(void)
{
  show_parameters(DV_params);
}
  
/* get Source at (tKS, xKS,yKS,zKS) */
double DV_SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS)
{
  return -SourceInKerrSchild(tKS, xKS,yKS,zKS, DV_params);
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
