/* TestID.c */
/* Wolfgang Tichy 09/2009 */

#include "sgrid.h"
#include "TestID.h"



/* compute TestID data */
int TestID(tGrid *grid)
{
  /* allocate memory for psi and its derivs, and for Var List psiandderivs */
  psiandderivs->grid = grid;
  enablevarlist(psiandderivs);

  /* allocate memory for ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));

  printf("TestID: Computing initial data for testing:\n");
  printf("  TestID_type = %s\n", Gets("TestID_type"));

  if(Getv("TestID_type", "1dGaugeWave"))
  TestID_1dGaugeWave(grid, Ind("x"), Ind("gxx"), Ind("Kxx"), 
                     Ind("psi"), Ind("dpsiopsix"), Ind("ddpsiopsixx"), 
                     Ind("alpha"), Ind("betax") );
  else
    errorexit("unknown TestID_type");
    
  /* set initial TrK for (needed) for some gauges in BSSN */
  set_K_initial(grid);

  return 0;
}
