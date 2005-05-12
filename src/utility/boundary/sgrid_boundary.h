/* sgrid_boundary.h */
/* Wolfgang Tichy 5/2005 */

/* global point lists to mark boundries */
extern tPointList *radiativeBoundaryPointList;
extern tPointList *simpleExcisionBoundaryPointList;

/* functions to set BCs */
void set_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void set_boundary_radiative(tGrid *grid, 
			    int unew, int upre, double c, int ucur,
			    double var0, double v);
void set_boundary_excision(tGrid *grid, int unew, int upre);
