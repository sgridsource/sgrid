/* testwave.h */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 6/02 */


void testwave_evo(tVarList *vlunew,
                  tVarList *vlupre, double dt, tVarList *vlucur);
int testwave_startup(tGrid *grid);
int testwave_analyze(tGrid *grid);
int testwave_filter(tGrid *grid);
