/* sgrid_boundary.h */
/* Wolfgang Tichy 5/2005 */

/* global point lists to mark boundries */
extern tPointList *radiativeBoundaryPointList;
extern tPointList *constantBoundaryPointList;
extern tPointList *selectedBoundaryPointList;
extern tPointList *ExcisionBoundaryPointList;
extern tPointList *boxBoundaryPointList;

/* functions to set BCs */
void set_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void set_boundary_radiative(tGrid *grid, 
			    int unew, int upre, double c, int ucur,
			    double var0, double v);
void set_boundary_simpleExcision(tGrid *grid, int unew, int upre);
void set_boundary_constant(tPointList *PL, int unew, int upre);
void set_boundary_value(tPointList *PL, int ui, double value);
void set_boundary_normalderiv_leftBound(tPointList *PL, int direc, int ui, 
                                        double deriv);
void set_boundary_normalderiv_rightBound(tPointList *PL, int direc, int ui, 
                                         double deriv);
void set_boundary_VonNeumannExcision(tPointList *PL, int unew);
