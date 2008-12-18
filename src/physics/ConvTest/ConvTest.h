/* ConvTest.h */
/* Wolfgang Tichy, April 2005 */


void ConvTest_evo(tVarList *vlunew,
                  tVarList *vlupre, double dt, tVarList *vlucur);
int ConvTest_startup(tGrid *grid);
int ConvTest_analyze(tGrid *grid);
int ConvTest_filter(tGrid *grid);
