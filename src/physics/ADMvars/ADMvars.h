/* ADMvars.h */
/* Wolfgang Tichy, 2005 */

void ADMconstraints(tVarList *u);
int allocateADMvars(tGrid *grid);
int ExitIfNAN(tGrid* grid);
void ADMenergy_spheric_intergrand(tVarList *u);
