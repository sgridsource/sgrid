/* ADMvars.h */
/* Bernd Bruegmann, 6/02 */

void ADMconstraints(tVarList *u);
int computeADMconstraints(tGrid *grid);
int set_K_initial(tGrid *grid);
int CheckIfFinite(tGrid* grid, char *varname);
int ExitIfNAN(tGrid* grid);
