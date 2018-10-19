/* sgrid_Coordinates.h */
/* Wolfgang Tichy, 4/2005 */

#include "Coordinates.h"


/* For coordtrans_CubedSphere.c :
   Type of cubed sphere or rather sphered cube coord transform */
enum
{
  CoordInfoNotSet,  /* if box->CI->type is not set, box->CI->type=0 */
  PyramidFrustum,   /* both inner & outer surfaces are flat */  
  innerCubedSphere, /* inner surface is curved, but outer surface is flat */
  outerCubedSphere, /* outer surface is curved, but inner surface is flat */
  CubedShell        /* both inner & outer surfaces are curved */
};


/* Coordinates.c */
int init_CoordTransform_And_Derivs(tGrid *grid);
int compute_xyz_dXYZdxyz_ddXYZddxyz(tGrid *grid);
int set_box_CI_struct_from_pars(tGrid *grid);
int set_box_CI_pars_from_box_CI_struct(tGrid *grid);

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
int moveXYZ_off_face(tBox *box, double *X, double *Y, double *Z);
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
int set_bits_in_all_bfaces_VER0(tGrid *grid);
int set_bits_in_all_bfaces(tGrid *grid);
int set_touching_bfaces_of_boxes_with_same_facepoints(tGrid *grid, int b0, int nb);
void set_all_bfaces_with_ob_minus1_to_outerbound(tGrid *grid, int b0, int nb);
int set_oXi_oYi_oZi_in_all_bfaces(tGrid *grid);
int set_oX_oY_oZ_vars_for_bfaces(tGrid *grid);
void find_external_faces_of_box(tBox *box, int *extface, int inlcOuterBound);
int remove_bfacepoints_that_cause_inconsistent_touch_bits(tGrid *grid);
int set_consistent_flags_in_all_bfaces(tGrid *grid);

/* from populate_bfaces_AR.c */
int populate_bfaces(tGrid *grid);

/* from coordtrans_CubedSphere.c */
int r_dr_dlam_of_lamAB_CubSph(tBox *box, int ind, double lam,
                              double A, double B, double *r, double *drdlam);
int ThetaPhi_of_AB_CubSph(tBox *box, double A, double B,
                          double *Theta, double *Phi);
int ThetaPhi_dThetaPhidAB_of_AB_CubSph(tBox *box, double A, double B,
                                       double *Theta,    double *Phi,
                                       double *dThetadA, double *dThetadB,
                                       double *dPhidA,   double *dPhidB);

/* from setup_CubedSpheres.c */
void init_6CubedSphereBoxes_from_CI_iFS(tGrid *grid, int bi_dom0);
void copy_CubedSphere_sigma01_inplane(tBox *box, int si,
                                      int isigma_src, int isigma_dest);
void disable_and_reset_CI_iSurf_vars(tBox *box);
void disable_Coordinates_CubedSphere_sigma01(tBox *box);
void compute_CubedSphere_dsigma01(tBox *box, int isigma,
                                  int isigma_dA, int isigma_dB);
int arrange_1box12CubSph_into_full_cube(tGrid *grid, int b0, double *xc,
                                        double din, double dmid, double dout);
int two_full_cubes_touching_at_x0(tGrid *grid, int b0, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2);
int sphere_around_two_full_cubes_touching_at_x0(tGrid *grid, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0);
int two_spheres_around_two_full_cubes(tGrid *grid, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1);

/* from find_extrema.c */
int box_extremum_of_F(tBox *box, int Fi,
                      double *X, double *Y, double *Z, double *Fextr);
int box_extremum_of_F_in_dir(tBox *box, int Fi, int dir, double C1, double C2,
                             double *C, double *Fextr);

/* from find_ijk.c */
int find_i_Of_X(tBox *box, double X0);
int find_j_Of_Y(tBox *box, double Y0);
int find_k_Of_Z(tBox *box, double Z0);

/* from FSurf_CubedSpheres.c */
int FSurf_CubSph_init6Boxes_from_CI_iFS(tGrid *grid, int bi_dom0);
