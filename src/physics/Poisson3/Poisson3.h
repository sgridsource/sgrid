/* Poisson3.h */
/* Wolfgang Tichy, Oct. 2017 */


int Poisson3_startup(tGrid *grid);
int Poisson3_analyze(tGrid *grid);
int Poisson3_solve(tGrid *grid);

void F_Poisson3(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_Poisson3(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);
