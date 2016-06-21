/* Poisson2.h */
/* Wolfgang Tichy, June 2016 */


int Poisson2_startup(tGrid *grid);
int Poisson2_analyze(tGrid *grid);
int Poisson2_solve(tGrid *grid);

void F_Poisson2(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_Poisson2(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);
