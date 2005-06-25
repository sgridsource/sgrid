/* Z4secondO.c */
/* Wolfgang Tichy 6/2005 */

#include "sgrid.h"
#include "Z4secondO.h"

tVarList *Z4secondOvars;

/* evolve in the interior and 
   for those boundary points set by special evolution
*/ 
void Z4secondO_evolve(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)
{
  Z4secondO_rhs(unew, upre, dt, ucur);
}



/* initialize Z4secondO after initial data has been 
   computed in POST_INITIALDATA */
int Z4secondO_startup(tGrid *grid)
{
  double vgauge;

  printf("Initializing Z4secondO:\n");

  /* set gauge speed for lapse and related quantities */
  if (Getv("Z4secondO_lapse", "1+log"))
    vgauge = sqrt(Getd("Z4secondO_lapseharmonicf"));
  else
    vgauge = 1;

  /* set boundary information for Z4secondO evolution: 
     farlimit, falloff, propagation speed 
  */
  VarNameSetBoundaryInfo("gxx", 1, 1, 1.0);
  VarNameSetBoundaryInfo("gxy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("gxz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("gyy", 1, 1, 1.0);
  VarNameSetBoundaryInfo("gyz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("gzz", 1, 1, 1.0);

  VarNameSetBoundaryInfo("Kxx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Kxy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Kxz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Kyy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Kyz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Kzz", 0, 1, 1.0);

  VarNameSetBoundaryInfo("Z4secondO_Zx",  0, 1, 1.0);
  VarNameSetBoundaryInfo("Z4secondO_Zy",  0, 1, 1.0);
  VarNameSetBoundaryInfo("Z4secondO_Zz",  0, 1, 1.0);
  VarNameSetBoundaryInfo("Z4secondO_Theta",   0, 1, 1.0);

  VarNameSetBoundaryInfo("alpha",    1, 1, vgauge);
  VarNameSetBoundaryInfo("betax",    0, 1, 1.0); 
  VarNameSetBoundaryInfo("betay",    0, 1, 1.0);
  VarNameSetBoundaryInfo("betaz",    0, 1, 1.0);
  VarNameSetBoundaryInfo("betadotx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("betadoty", 0, 1, 1.0);
  VarNameSetBoundaryInfo("betadotz", 0, 1, 1.0);

  /* create a variable list for Z4secondO evolutions 
     note that we include lapse and shift directly
  */
  Z4secondOvars = vlalloc(grid);
  vlpush(Z4secondOvars, Ind("gxx"));
  vlpush(Z4secondOvars, Ind("Kxx"));
  vlpush(Z4secondOvars, Ind("Z4secondO_Theta"));
  vlpush(Z4secondOvars, Ind("Z4secondO_Zx"));
  vlpush(Z4secondOvars, Ind("alpha"));
  vlpush(Z4secondOvars, Ind("betax"));
  vlpush(Z4secondOvars, Ind("betadotx"));
    
  if (0) prvarlist(Z4secondOvars);
  enablevarlist(Z4secondOvars);

  /* register evolved variables */
  evolve_vlregister(Z4secondOvars);
  
  /* register evolution routine */
  evolve_rhsregister(Z4secondO_evolve);

  /* enable all derivative vars */
  enablevar(grid, Ind("Z4secondO_dThetax"));
  enablevar(grid, Ind("Z4secondO_dZxx"));
  enablevar(grid, Ind("Z4secondO_dalpx"));
  enablevar(grid, Ind("Z4secondO_ddalpxx"));
  enablevar(grid, Ind("Z4secondO_dbetaxx"));
  enablevar(grid, Ind("Z4secondO_ddbetaxxx"));

  return 0;
}
