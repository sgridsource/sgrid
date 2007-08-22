/* sgrid_SingleBH.c */
/* Wolfgang Tichy, 12/02 */

#include "sgrid.h"
#include "SingleBH.h"


int sgrid_SingleBH() 
{
  if (!Getv("physics", "SingleBH")) return 0;
  printf("Adding SingleBH\n");

  /* functions: */
  AddFun(INITIALDATA, SingleBH, "SingleBH: set initial data for a single BH");

  /* Parameters: */
  AddPar("BHmass1", "0.0", "mass of black hole 1");
  AddPar("BHsx1", "0.0", "spin_x of black hole 1");
  AddPar("BHsy1", "0.0", "spin_y of black hole 1");
  AddPar("BHsz1", "0.0", "spin_z of black hole 1");

  AddPar("SingleBH_type", "KerrSchild",
         "on what slice and in which coords we compute BH data"
         " [KerrSchild, isotropic]");
  AddPar("SingleBH_initial_lapse", "donothing",
         "initial lapse [donothing,one]");
  AddPar("SingleBH_initial_shift", "donothing",
         "initial shift [donothing,zero]");
  AddPar("SingleBH_ConformalFactor","no","set cconformal factor psi=1+2M/r");

  return 0;
}
