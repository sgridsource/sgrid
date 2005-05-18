/* sgrid_boundary.c */
/* Bernd Bruegmann 01/00 */

#include "sgrid.h"
#include "boundary.h"


/* global point lists to mark boundries */
tPointList *radiativeBoundaryPointList;
tPointList *constantBoundaryPointList;
tPointList *selectedconstantBoundaryPointList;
tPointList *ExcisionBoundaryPointList;
tPointList *boxBoundaryPointList;


int sgrid_boundary() 
{
  printf("Adding boundary\n");

  /* functions */
  AddFun(PRE_INITIALDATA, initialize_BoundaryPointLists,
         "initialize all Boundary PointLists");
         
  /* parameters */
  AddPar("boundary", "", "boundary condition"
         "[radiative,constant,constantExcision,simpleExcision,"
          "selectedconstantExcision]");

  /* radiative boundary */
  if (Getv("boundary", "radiative"))
  {
    AddPar("boundary_radpower", "0", "exponent of non-wave term");
    AddPar("boundary_radconstant", "", 
	   "keep these variables constant for radiative boundary");
  }

  /* selectedconstantExcision boundary */
  if (Getv("boundary", "selectedconstantExcision"))
  {
    AddPar("boundary_selectedconstantExcision", "", 
	   "keep these variables constant on selectedconstantExcision boundary");
  }

  return 0;
}
