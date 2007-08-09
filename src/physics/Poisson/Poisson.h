/* Poisson.h */
/* Wolfgang Tichy, April 2005 */


int Poisson_startup(tGrid *grid);
int Poisson_analyze(tGrid *grid);
int Poisson_solve(tGrid *grid);

void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_Poisson(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);
void Precon_I(tVarList *vlJdu, tVarList *vldu,
              tVarList *vlduDerivs, tVarList *vlu);
