/* ScalarOnKerr.h */
/* Wolfgang Tichy  8/2005 */


int ScalarOnKerr_startup(tGrid *grid);
void ScalarOnKerr_evolve(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur);
int ScalarOnKerr_analyze(tGrid *grid);
void Kerr(tGrid *grid, int i_x, int i_g, int i_gup, int i_Gam);
