/* sgrid_MemoryMan.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "MemoryMan.h"



int sgrid_MemoryMan(void) 
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

    /* number of points in all 3 directions */
    snprintf(str, 999, "box%d_n1", b);
    AddPar(str, "n1", "number of points in X-direction");
    snprintf(str, 999, "box%d_n2", b);
    AddPar(str, "n2", "number of points in Y-direction");
    snprintf(str, 999, "box%d_n3", b);
    AddPar(str, "n3", "number of points in Z-direction");

    /* min and max values in all 3 coords */
    snprintf(str, 999, "box%d_min1", b);
    AddPar(str, "min1", "lower boundary in X-direction");
    snprintf(str, 999, "box%d_max1", b);
    AddPar(str, "max1", "upper boundary in X-direction");
    
    snprintf(str, 999, "box%d_min2", b);
    AddPar(str, "min2", "lower boundary in Y-direction");
    snprintf(str, 999, "box%d_max2", b);
    AddPar(str, "max2", "upper boundary in Y-direction");
    
    snprintf(str, 999, "box%d_min3", b);
    AddPar(str, "min3", "lower boundary in Z-direction");
    snprintf(str, 999, "box%d_max3", b);
    AddPar(str, "max3", "upper boundary in Z-direction");
    
    /* basis functions used */
    snprintf(str, 999, "box%d_basis1", b);
    AddPar(str, "basis1", 
                "basis functions in X-direction [ChebExtrema, ChebZeros]");
    snprintf(str, 999, "box%d_basis2", b);
    AddPar(str, "basis2", 
                "basis functions in Y-direction [ChebExtrema, ChebZeros]");
    snprintf(str, 999, "box%d_basis3", b);
    AddPar(str, "basis3", 
                "basis functions in Z-direction [ChebExtrema, ChebZeros]");

    /* simple filter amounts used */
    snprintf(str, 999, "box%d_filter1", b);
    AddPar(str, "0", 
                "filter amount in X-direction [some number up to n1]");
    snprintf(str, 999, "box%d_filter2", b);
    AddPar(str, "0", 
                "filter amount in Y-direction [some number up to n2]");
    snprintf(str, 999, "box%d_filter3", b);
    AddPar(str, "0", 
                "filter amount in Z-direction [some number up to n3]");
  }
  return 0;
}
