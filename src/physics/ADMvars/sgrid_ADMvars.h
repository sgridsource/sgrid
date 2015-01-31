/* sgrid_ADMvars.h */
/* Wolfgang Tichy, 2005 */

extern tVarList *psiandderivs, *K_initial;


/* functions to deal with ADM vars */
int set_K_initial(tGrid *grid);
int ADMvars_undo_conformal_split(tGrid *grid);
int CheckIfFinite(tGrid* grid, char *varname);
int computeADMconstraints(tGrid *grid);
