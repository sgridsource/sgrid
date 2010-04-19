/* sgrid_PN_CircularOrbit_GWs.c */
/* Wolfgang Tichy, 4/2010 */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"


int sgrid_PN_CircularOrbit_GWs() 
{
  if (!Getv("physics", "PN_CircularOrbit_GWs")) return 0;
  printf("Adding PN_CircularOrbit_GWs\n");

  /* functions: */
  AddFun(PRE_GRID, PN_CircularOrbit_GWs_set_boxsize, "setup box size");
  AddFun(INITIALDATA, PN_CircularOrbit_GWs_startup, "startup");
  AddFun(ANALYZE, PN_CircularOrbit_GWs, "PN_CircularOrbit_GWs: read data and compute modes");

  /* variables */
  AddVar("PN_CircularOrbit_GWs_hplus" ,  "", "h+ on a sphere, note H = (h+ - ihx)");
  AddVar("PN_CircularOrbit_GWs_hcross",  "", "hx on a sphere, note H = (h+ - ihx)");
  AddVar("PN_CircularOrbit_GWs_Re_Hmode","", "real part of H = (h+ - ihx) mode");
  AddVar("PN_CircularOrbit_GWs_Im_Hmode","", "imag part of H = (h+ - ihx) mode");
  AddVar("PN_CircularOrbit_GWs_Re_Psi4",     "", "real part of Psi4");
  AddVar("PN_CircularOrbit_GWs_Im_Psi4",     "", "imag part of Psi4");
  AddVar("PN_CircularOrbit_GWs_Re_Psi4mode", "", "real part of Psi4 mode");
  AddVar("PN_CircularOrbit_GWs_Im_Psi4mode", "", "imag part of Psi4 mode");
  AddVar("PN_CircularOrbit_GWs_Re_sYlm", "", "real part of spin-weighted Ylm");
  AddVar("PN_CircularOrbit_GWs_Im_sYlm", "", "imaginary part spin-weighted Ylm");

  /* Parameters: */
  AddPar("PN_CircularOrbit_GWs_outfile_prefix", "h_",
         "imaginary part of data on sphere");
  AddPar("PN_CircularOrbit_GWs_lmax", "6", "max l mode we compute");
  AddPar("PN_CircularOrbit_GWs_sphere_Lmax", "13",
         "max l mode we can represent on grid. We need Lmax>=lmax");
  AddPar("PN_CircularOrbit_GWs_spinweight", "-2", "spin weight we use in sYlm");

  return 0;
}
