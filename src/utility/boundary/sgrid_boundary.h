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
void set_boundary_radiative_analytic(tVarList *unew, tVarList *upre,
                                     double c, tVarList *ucur);

/* from BCs_from_bfaces.c */
void boxface_normal_at_ijk(tBox *box, int f, int ijk, double n[4]);
void FPsi_1Dinterp_for_bface(int iFPsi, tBface *bface, int idir,
                             int iPsi, int idPsi[4]);
void FPsi_2Dinterp_for_bface(int iFPsi, tBface *bface, int plN,
                             int iPsi, int idPsi[4]);
void FPsi_3Dinterp_for_bface(int iFPsi, tBface *bface, int iPsi, int idPsi[4]);
int ijk_in_other_box_if_same_fpts(tBface *bface, int pi);
void FPsi_copy_for_bface(int iFPsi, tBface *bface, int iPsi, int idPsi[4]);
