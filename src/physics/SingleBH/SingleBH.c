/* SingleBH.c */
/* Wolfgang Tichy 12/02 */

#include "sgrid.h"
#include "SingleBH.h"



/* compute SingleBH data */
int SingleBH(tGrid *grid)
{
  /* allocate memory for psi and its derivs, and for Var List psiandderivs */
  psiandderivs->grid = grid;
  enablevarlist(psiandderivs);

  /* allocate memory for ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));

  printf("Computing Scwarzschild initial data:\n");
  printf("  BHmass = %e\n", Getd("BHmass1"));

  /* compute */
  SingleBHKS(grid, Ind("x"), Ind("gxx"), Ind("Kxx"), 
                        Ind("psi"), Ind("dpsiopsix"), Ind("ddpsiopsixx"), 
                        Ind("alpha"), Ind("betax") );

  /* set initial TrK for (needed) for some gauges in BSSN */
  set_K_initial(grid);

  /* temporary hack */
  // test_bampi_getdata(level);

  return 0;
}
