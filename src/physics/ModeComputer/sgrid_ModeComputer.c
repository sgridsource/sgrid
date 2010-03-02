/* sgrid_ModeComputer.c */
/* Wolfgang Tichy, 2/2010 */

#include "sgrid.h"
#include "ModeComputer.h"


int sgrid_ModeComputer() 
{
  if (!Getv("physics", "ModeComputer")) return 0;
  printf("Adding ModeComputer\n");

  /* functions: */
  AddFun(PRE_GRID, ModeComputer_set_boxsize, "setup box size");
  AddFun(INITIALDATA, ModeComputer_startup, "startup");
  AddFun(ANALYZE, ModeComputer, "ModeComputer: read data and compute modes");

  /* variables */
  AddVar("ModeComputer_Re_var",  "", "real part of variable");
  AddVar("ModeComputer_Im_var",  "", "imaginary part of variable");
  AddVar("ModeComputer_Re_mode", "", "real part of mode");
  AddVar("ModeComputer_Im_mode", "", "imaginary part of mode");

  /* Parameters: */
  AddPar("ModeComputer_Re_sphere_data", "rpsi4.r3.l3",
         "real part of data on sphere");
  AddPar("ModeComputer_Im_sphere_data", "ipsi4.r3.l3",
         "imaginary part of data on sphere");
  AddPar("ModeComputer_outfile_prefix", "psi4_r3_",
         "imaginary part of data on sphere");
  AddPar("ModeComputer_lmax", "6", "max l mode we compute");
  AddPar("ModeComputer_spinweight", "-2", "spin weight we use in sYlm");
  AddPar("ModeComputer_Re_bitantsym","+1", "symmetry factor for bitant of real part");
  AddPar("ModeComputer_Re_rotantsym","+1", "symmetry factor for rotant of real part");
  AddPar("ModeComputer_Re_octantsym","+1", "symmetry factor for octant of real part");
  AddPar("ModeComputer_Im_bitantsym","-1", "symmetry factor for bitant of imaginary part");
  AddPar("ModeComputer_Im_rotantsym","+1", "symmetry factor for rotant of imaginary part");
  AddPar("ModeComputer_Im_octantsym","-1", "symmetry factor for octant of imaginary part");

  return 0;
}
