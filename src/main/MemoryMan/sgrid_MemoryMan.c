/* sgrid_MemoryMan.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "MemoryMan.h"



int sgrid_MemoryMan() 
{
  int nboxes;
  int b;
  
  printf("Adding MemoryMan\n");

  /* variables */
  AddConstantVar("X", "",
  "coordinate 1 used for output and in spectral expansion e.g. rho");
  AddConstantVar("Y", "",
  "coordinate 2 used for output and in spectral expansion e.g. theta");
  AddConstantVar("Z", "",
  "coordinate 3 used for output and in spectral expansion e.g. z");

  /* parameters */
  AddPar("storage_verbose", "no", 
	 "verbose mode for memory allocation [no,yes]");
  AddPar("nboxes", "1", "number of boxes in grid");
  nboxes = Geti("nboxes");
  for(b=0; b<nboxes; b++)
  {
    char str[1000];

    /* number of popint in all 3 directions */
    snprintf(str, 999, "box%d_n1", b);
    AddPar(str, "1", "number of points in x-direction");
    snprintf(str, 999, "box%d_n2", b);
    AddPar(str, "1", "number of points in y-direction");
    snprintf(str, 999, "box%d_n3", b);
    AddPar(str, "1", "number of points in z-direction");

    /* min and max values in all 3 coords */
    snprintf(str, 999, "box%d_min1", b);
    AddPar(str, "-1", "lower boundary in x-direction");
    snprintf(str, 999, "box%d_max1", b);
    AddPar(str, "+1", "upper boundary in x-direction");
    
    snprintf(str, 999, "box%d_min2", b);
    AddPar(str, "-1", "lower boundary in y-direction");
    snprintf(str, 999, "box%d_max2", b);
    AddPar(str, "+1", "upper boundary in y-direction");
    
    snprintf(str, 999, "box%d_min3", b);
    AddPar(str, "-1", "lower boundary in z-direction");
    snprintf(str, 999, "box%d_max3", b);
    AddPar(str, "+1", "upper boundary in z-direction");
    
    /* basis functions used */
    snprintf(str, 999, "box%d_basis1", b);
    AddPar(str, "ChebExtrema", 
                "basis functions in x-direction [ChebExtrema, ChebZeros]");
    snprintf(str, 999, "box%d_basis2", b);
    AddPar(str, "ChebExtrema", 
                "basis functions in y-direction [ChebExtrema, ChebZeros]");
    snprintf(str, 999, "box%d_basis3", b);
    AddPar(str, "ChebExtrema", 
                "basis functions in z-direction [ChebExtrema, ChebZeros]");
  }
  return 0;
}
