/* sgrid_Coordinates.h */
/* Wolfgang Tichy, 4/2005 */

#include "Coordinates.h"


/* Coordinates.c */
int init_CoordTransform_And_Derivs(tGrid *grid);
int compute_xyz_dXYZdxyz_ddXYZddxyz(tGrid *grid);

/* cartesianDerivs.c */
void cart_partials(tBox *box, double *u, double *u1, double *u2, double *u3);
void cart_partial_all(tBox *box, double *u, double *u1, double *u2, double *u3,
                      double *u11, double *u12, double *u13,
                      double *u22, double *u23, double *u33 );

/* Helper functions to set partial derivatives of a symmetric tensor */
void FirstDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc);
void FirstAndSecondDerivsOf_Sab(tBox *box, int i_Sab,
                                int i_dSabc, int i_ddSabcd);
void allDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc, int i_ddSabcd);
void D_and_DD_of_Sab(tBox *box, int i_Sab, int i_dSabc, int i_ddSabcd);

/* Helper functions to set partial derivatives of a vector */
void FirstDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab);
void FirstAndSecondDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc);
void allDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc);
void D_and_DD_of_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc);

/* Helper functions to set partial derivatives of a scalar */
void FirstDerivsOf_S(tBox *box, int i_S, int i_dSa);
void FirstAndSecondDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab);
void allDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab);
void D_and_DD_of_S(tBox *box, int i_S, int i_dSa, int i_ddSab);

/* doubleCovering.c */
void reset_doubleCoveredPoints(tVarList *unew);
void copy_to_doubleCoveredPoints_SphericalDF(tBox *box, int vind);

/* coordFilters.c */
void coordinateDependentFilter(tVarList *unew);

/* from ComplexFunctions.c */
double BaseAngle(double p, double peri, double p0);
double Arg(double x, double y);
double Arg_plus(double x, double y);
double dArgdx(double x, double y);
double dArgdy(double x, double y);
double ddArgdxdx(double x, double y);
double ddArgdxdy(double x, double y);
double ddArgdydy(double x, double y);
double Retanh(double x, double y);
double Imtanh(double x, double y);
double Argtanh(double x, double y);
double Abstanh(double x, double y);
double Argsech(double x, double y);
double Abssech(double x, double y);
double dArgtanhdx(double x, double y);
double dArgtanhdy(double x, double y);
double ddArgtanhdxdx(double x, double y);
double ddArgtanhdxdy(double x, double y);
double ddArgtanhdydy(double x, double y);
double dAbstanhdx(double x, double y);
double dAbstanhdy(double x, double y);
double ddAbstanhdxdx(double x, double y);
double ddAbstanhdxdy(double x, double y);
double ddAbstanhdydy(double x, double y);

/* from findXYZ_of_xyz.c */
int XYZ_of_xyz(tBox *box, double *X, double *Y, double *Z,
               double x, double y, double z);
int b_XYZ_of_xyz(tGrid *grid, double *X, double *Y, double *Z,
                 double x, double y, double z);
int b_XYZ_of_xyz_inboxlist(tGrid *grid, int *blist, int nb,
                           double *X, double *Y, double *Z,
                           double x, double y, double z);
double nearestXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                         double x, double y, double z);
double nearest_b_XYZ_of_xyz(tGrid *grid,  int *b, int *ind,
                            double *X, double *Y, double *Z,
                            double x, double y, double z);
double nearest_b_XYZ_of_xyz_inboxlist(tGrid *grid, int *blist, int nb, 
                            int *b, int *ind,
                            double *X, double *Y, double *Z,
                            double x, double y, double z);
double nearestinnerXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                              double x, double y, double z);
int XYZ_on_face(tBox *box, int *face, double X, double Y, double Z);
double guessXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                       double x, double y, double z);
double nearestXYZ_of_xyz_inplane(tBox *box, int *ind, 
                                 double *X, double *Y, double *Z,
                                 double x, double y, double z,
                                 int plane, int pind);
int X_of_x_forgiven_YZ(tBox *box, double *X, double x, double Y, double Z);
int b_X_of_x_forgiven_YZ(tGrid *grid, double *X, double x, double Y, double Z);
int b_X_of_x_forgiven_YZ_inboxlist(tGrid *grid, int *blist, int nb, 
                                   double *X, double x, double Y, double Z);
int Y_of_y_forgiven_XZ(tBox *box, double *Y, double y, double X, double Z);

/* from Coordinates_set_bfaces.c */
int Coordinates_set_bfaces(tGrid *grid);
int set_ofi_in_all_bfaces(tGrid *grid);
int set_bits_in_all_bfaces(tGrid *grid);
int set_touching_bfaces_of_boxes_with_same_facepoints(tGrid *grid, int b0, int nb);
