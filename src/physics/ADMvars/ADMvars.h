/* ADMvars.h */
/* Bernd Bruegmann, 6/02 */

void ADMconstraints(tVarList *u);
int computeADMconstraints(tGrid *grid);
int allocateADMvars(tGrid *grid);
int ExitIfNAN(tGrid* grid);
void ADMenergy_spheric_intergrand(tVarList *u);
