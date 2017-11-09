/* Coordinates.h */
/* Wolfgang Tichy 4/2005 */

/* use we can use long double in some places */
#define LDOUBLE long double


/* from cartesianDerivs.c */
double dX_dx(int i, int m, int dXd,
             tBox *box, int ind, double X, double Y, double Z);
double ddX_dxdx(int i, int m, int n, int ddXdd,
                tBox *box, int ind, double X, double Y, double Z);

/* Singularities.c */
void prSingInfo(tSingInfo *si);
int isSing_AnsorgNS12(void *aux, double X, double Y, double Z,
                      int update, tSingInfo *si);
int isSing_AnsorgNS03(void *aux, double X, double Y, double Z,
                      int update, tSingInfo *si);
/* from findXYZ_of_xyz.c */
int check_xyz_error(tBox *box, double *X, double *Y, double *Z,
            double x, double y, double z,
            void *p, double tol, tSingInfo *si, double *err, int pr);
int recover_if_start_on_singularity(tBox *box,
        double *X, double *Y, double *Z, double x, double y, double z,
        void *p, double tol, tSingInfo *si, double *err);

/* from doubleCovering.c */
void reset_doubleCoveredPoints_SphericalDF(tBox *box, tVarList *unew);

/* from coordFilters.c */
void coordinateDependentFilter_SphericalDF(tBox *box, tVarList *unew);
void coordinateDependentFilter_Spherical(tBox *box, tVarList *unew);

/* from transfDerivs.c */
void dXdx_from_dxdX(double dXdx[4][4], double dxdX[4][4]);
void dxdX_from_dXdx(double dxdX[4][4], double dXdx[4][4]);
void ddXdxdx_from_dXdx_ddxdXdX(double ddXdxdx[4][4][4],
                               double dXdx[4][4], double ddxdXdX[4][4][4]);
void dxdX_from_dxdU_dUdX(double dxdX[4][4],
                         double dxdU[4][4], double dUdX[4][4]);
void ddxdXdX_from_dxdU_dUdX_ddxdUdU_ddUdXdX(double ddxdXdX[4][4][4],
                         double dxdU[4][4], double dUdX[4][4],
                         double ddxdUdU[4][4][4], double ddUdXdX[4][4][4]);
double check_box_dx_dX(tBox *box, double X, double Y, double Z);

/* Some trivial functions */
double zero_of_xyz(void *aux, int ind, double X, double Y, double Z);
double one_of_xyz(void *aux, int ind, double X, double Y, double Z);
double x_equals_X(void *aux, int ind, double X, double Y, double Z);
double y_equals_Y(void *aux, int ind, double X, double Y, double Z);
double z_equals_Z(void *aux, int ind, double X, double Y, double Z);

/* Polar coordinates: */
double x_ofPolar(void *aux, int ind, double rho, double phi, double Z);
double y_ofPolar(void *aux, int ind, double rho, double phi, double Z);
double z_ofPolar(void *aux, int ind, double rho, double phi, double Z);
double drho_dx(void *aux, int ind, double rho, double phi, double Z);
double drho_dy(void *aux, int ind, double rho, double phi, double Z);
double dphi_dx(void *aux, int ind, double rho, double phi, double Z);
double dphi_dy(void *aux, int ind, double rho, double phi, double Z);
void set_d_dx_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
void set_d_dy_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
/*
double drho_dxdx(void *aux, int ind, double rho, double phi, double Z);
double drho_dxdy(void *aux, int ind, double rho, double phi, double Z);
double drho_dydy(void *aux, int ind, double rho, double phi, double Z);
double dphi_dxdx(void *aux, int ind, double rho, double phi, double Z);
double dphi_dxdy(void *aux, int ind, double rho, double phi, double Z);
double dphi_dydy(void *aux, int ind, double rho, double phi, double Z);
*/

/* PolarCE coordinates: */
double x_ofPolarCE(void *aux, int ind, double rho, double Y, double Z);
double y_ofPolarCE(void *aux, int ind, double rho, double Y, double Z);
double dYPolarCE_dx(void *aux, int ind, double rho, double Y, double Z);
double dYPolarCE_dy(void *aux, int ind, double rho, double Y, double Z);

/* SphericalDF coordinates: */
double x_ofSphericalDF(void *aux, int ind, double r, double thm, double phi);
double y_ofSphericalDF(void *aux, int ind, double r, double thm, double phi);
double z_ofSphericalDF(void *aux, int ind, double r, double thm, double phi);
double dr_dx(void *aux, int ind, double r, double thm, double phi);
double dr_dy(void *aux, int ind, double r, double thm, double phi);
double dr_dz(void *aux, int ind, double r, double thm, double phi);
double dphiSphericalDF_dx(void *aux, int ind, double r, double thm, double phi);
double dphiSphericalDF_dy(void *aux, int ind, double r, double thm, double phi);
double dthm_dx(void *aux, int ind, double r, double thm, double phi);
double dthm_dy(void *aux, int ind, double r, double thm, double phi);
double dthm_dz(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi);
double ddr_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi);
double ddthm_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dxdx(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dxdy(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dxdz(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dydy(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dydz(void *aux, int ind, double r, double thm, double phi);
double ddphi_SphericalDF_dzdz(void *aux, int ind, double r, double thm, double phi);

/* compactSphericalDF coordinates: */
double r_of_xi(double xi);
double dxi_dr(double xi);
double x_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi);
double y_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi);
double z_ofcompactSphericalDF(void *aux, int ind, double xi, double thm, double phi);
double dxi_dx(void *aux, int ind, double xi, double thm, double phi);
double dxi_dy(void *aux, int ind, double xi, double thm, double phi);
double dxi_dz(void *aux, int ind, double xi, double thm, double phi);
double dthmcompactSphericalDF_dx(void *aux, int ind, double xi, double thm, double phi);
double dthmcompactSphericalDF_dy(void *aux, int ind, double xi, double thm, double phi);
double dthmcompactSphericalDF_dz(void *aux, int ind, double xi, double thm, double phi);
double dphicompactSphericalDF_dx(void *aux, int ind, double xi, double thm, double phi);
double dphicompactSphericalDF_dy(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdx(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdy(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdz(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dydy(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dydz(void *aux, int ind, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dzdz(void *aux, int ind, double xi, double thm, double phi);

/* Spherical coordinates */
double x_ofSpherical(void *aux, int ind, double r, double theta, double phi);
double y_ofSpherical(void *aux, int ind, double r, double theta, double phi);
double z_ofSpherical(void *aux, int ind, double r, double theta, double phi);
double drSpherical_dx(void *aux, int ind, double r, double theta, double phi);
double drSpherical_dy(void *aux, int ind, double r, double theta, double phi);
double drSpherical_dz(void *aux, int ind, double r, double theta, double phi);
double dthetaSpherical_dx(void *aux, int ind, double r, double theta, double phi);
double dthetaSpherical_dy(void *aux, int ind, double r, double theta, double phi);
double dthetaSpherical_dz(void *aux, int ind, double r, double theta, double phi);
double dphiSpherical_dx(void *aux, int ind, double r, double theta, double phi);
double dphiSpherical_dy(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dxdx(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dxdy(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dxdz(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dydy(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dydz(void *aux, int ind, double r, double theta, double phi);
double ddr_Spherical_dzdz(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dxdx(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dxdy(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dxdz(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dydy(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dydz(void *aux, int ind, double r, double theta, double phi);
double ddtheta_Spherical_dzdz(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dxdx(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dxdy(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dxdz(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dydy(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dydz(void *aux, int ind, double r, double theta, double phi);
double ddphi_Spherical_dzdz(void *aux, int ind, double r, double theta, double phi);

/* Spherical2 coordinates */
double x_ofSpherical2(void *aux, int ind, double r, double U, double phi);
double y_ofSpherical2(void *aux, int ind, double r, double U, double phi);
double z_ofSpherical2(void *aux, int ind, double r, double U, double phi);
double dr_Spherical2_dx(void *aux, int ind, double r, double U, double phi);
double dr_Spherical2_dy(void *aux, int ind, double r, double U, double phi);
double dr_Spherical2_dz(void *aux, int ind, double r, double U, double phi);
double dU_Spherical2_dx(void *aux, int ind, double r, double U, double phi);
double dU_Spherical2_dy(void *aux, int ind, double r, double U, double phi);
double dU_Spherical2_dz(void *aux, int ind, double r, double U, double phi);
double dphi_Spherical2_dx(void *aux, int ind, double r, double U, double phi);
double dphi_Spherical2_dy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dxdx(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dxdy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dxdz(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dydy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dydz(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical2_dzdz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dxdx(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dxdy(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dxdz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dydy(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dydz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical2_dzdz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dxdx(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dxdy(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dxdz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dydy(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dydz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical2_dzdz(void *aux, int ind, double r, double U, double phi);

/* Spherical3 coordinates */
double x_ofSpherical3(void *aux, int ind, double r, double U, double phi);
double y_ofSpherical3(void *aux, int ind, double r, double U, double phi);
double z_ofSpherical3(void *aux, int ind, double r, double U, double phi);
double dr_Spherical3_dx(void *aux, int ind, double r, double U, double phi);
double dr_Spherical3_dy(void *aux, int ind, double r, double U, double phi);
double dr_Spherical3_dz(void *aux, int ind, double r, double U, double phi);
double dU_Spherical3_dx(void *aux, int ind, double r, double U, double phi);
double dU_Spherical3_dy(void *aux, int ind, double r, double U, double phi);
double dU_Spherical3_dz(void *aux, int ind, double r, double U, double phi);
double dphi_Spherical3_dx(void *aux, int ind, double r, double U, double phi);
double dphi_Spherical3_dy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dxdx(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dxdy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dxdz(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dydy(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dydz(void *aux, int ind, double r, double U, double phi);
double ddr_Spherical3_dzdz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dxdx(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dxdy(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dxdz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dydy(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dydz(void *aux, int ind, double r, double U, double phi);
double ddU_Spherical3_dzdz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dxdx(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dxdy(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dxdz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dydy(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dydz(void *aux, int ind, double r, double U, double phi);
double ddphi_Spherical3_dzdz(void *aux, int ind, double r, double U, double phi);

/* tan_stretch coordinates: */
double x_of_xs(double xs);
double dxs_dx(double xs);
double ddxs_dxdx(double xs);
double x_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double y_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double z_of_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double dxs_dx_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double dys_dy_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double dzs_dz_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double ddxs_dxdx_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double ddys_dydy_tan_stretch(void *aux, int ind, double xs, double ys, double zs);
double ddzs_dzdz_tan_stretch(void *aux, int ind, double xs, double ys, double zs);

/* AnsorgNS coordinates:                                               */
/* 4 domains: 0=inside NS+, 1=outside NS+, 2=outside NS-, 3=inside NS- */
void xyz_of_AnsorgNS(tBox *box, int ind, int domain,
                     double A, double BB, double phi, 
                     double *x, double *y, double *z, double *Xp, double *Rp);
double ddphi_dydy_AnsorgNS(void *aux, int ind, double A, double B, double phi);
double ddphi_dydz_AnsorgNS(void *aux, int ind, double A, double B, double phi);
double ddphi_dzdz_AnsorgNS(void *aux, int ind, double A, double B, double phi);
double AnsorgNS_sigma_pm(tBox *box, int ind, double B, double phi);
double AnsorgNS_dsigma_pm_dB(tBox *box, int ind, double B, double phi);
double AnsorgNS_dsigma_pm_dphi(tBox *box, int ind, double B, double phi);
double AnsorgNS_ddsigma_pm_dBdB(tBox *box, int ind, double B, double phi);
double AnsorgNS_ddsigma_pm_dBdphi(tBox *box, int ind, double B, double phi);
double AnsorgNS_ddsigma_pm_dphidphi(tBox *box, int ind, double B, double phi);
double AnsorgNS_sigma_p_one(tBox *box, int ind, double B, double phi);
double AnsorgNS_dsigma_zero(tBox *box, int ind, double B, double phi);
double AnsorgNS_sigma_m_one(tBox *box, int ind, double B, double phi);
double x_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double y_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double z_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dA_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dA_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dA_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dB_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dB_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dB_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dphi_dx_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dphi_dy_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double dphi_dz_AnsorgNS0(void *aux, int ind, double A, double B, double phi);

double x_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double y_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double z_of_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dA_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dA_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dA_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dB_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dB_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dB_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dphi_dx_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dphi_dy_AnsorgNS1(void *aux, int ind, double A, double B, double phi);
double dphi_dz_AnsorgNS1(void *aux, int ind, double A, double B, double phi);

double x_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double y_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double z_of_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dA_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dA_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dA_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dB_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dB_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dB_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dphi_dx_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dphi_dy_AnsorgNS2(void *aux, int ind, double A, double B, double phi);
double dphi_dz_AnsorgNS2(void *aux, int ind, double A, double B, double phi);

double x_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double y_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double z_of_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dA_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dA_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dA_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dB_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dB_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dB_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dphi_dx_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dphi_dy_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
double dphi_dz_AnsorgNS3(void *aux, int ind, double A, double B, double phi);
void set_d_dy_at_rhoEQzero_AnsorgNS(void *bo, void *va, 
                                    void *v1,void *v2,void *v3);
void set_d_dz_at_rhoEQzero_AnsorgNS(void *bo, void *va, 
                                    void *v1,void *v2,void *v3);

/* generic coord trafos */
void init_dXdx_generic(tBox *box);
void init_ddXdxdx_generic(tBox *box);
double dX_dx_generic(void *aux, int ind, double X, double Y, double Z);
double dX_dy_generic(void *aux, int ind, double X, double Y, double Z);
double dX_dz_generic(void *aux, int ind, double X, double Y, double Z);
double dY_dx_generic(void *aux, int ind, double X, double Y, double Z);
double dY_dy_generic(void *aux, int ind, double X, double Y, double Z);
double dY_dz_generic(void *aux, int ind, double X, double Y, double Z);
double dZ_dx_generic(void *aux, int ind, double X, double Y, double Z);
double dZ_dy_generic(void *aux, int ind, double X, double Y, double Z);
double dZ_dz_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dxdx_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dxdy_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dxdz_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dydy_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dydz_generic(void *aux, int ind, double X, double Y, double Z);
double ddX_dzdz_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dxdx_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dxdy_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dxdz_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dydy_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dydz_generic(void *aux, int ind, double X, double Y, double Z);
double ddY_dzdz_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dxdx_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dxdy_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dxdz_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dydy_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dydz_generic(void *aux, int ind, double X, double Y, double Z);
double ddZ_dzdz_generic(void *aux, int ind, double X, double Y, double Z);

/* from math */
void xyz_dxyz_ddxyz_of_AnsorgXRphi_XRphi(tBox *box, int ind, double X, double R, double phi, double xyz[4],double dxyz[4][4],double ddxyz[4][4][4]);
void XRphi_dXRphi_ddXRphi_of_AnsorgNS0_ABphi(tBox *box, int ind, double A, double B, double phi, double XRphi[4],double dXRphi[4][4],double ddXRphi[4][4][4]);
void XRphi_dXRphi_ddXRphi_of_AnsorgNS1_ABphi(tBox *box, int ind, double A, double B, double phi, double XRphi[4],double dXRphi[4][4],double ddXRphi[4][4][4]);
void XRphi_dXRphi_ddXRphi_of_AnsorgNS2_ABphi(tBox *box, int ind, double A, double B, double phi, double XRphi[4],double dXRphi[4][4],double ddXRphi[4][4][4]);
void XRphi_dXRphi_ddXRphi_of_AnsorgNS3_ABphi(tBox *box, int ind, double A, double B, double phi, double XRphi[4],double dXRphi[4][4],double ddXRphi[4][4][4]);
void ABphi_dABphi_ddABphi_of_Bfunc_AfBfPhi(tBox *box, int ind, double Af, double Bf, double Phi, double ABphi[4],double dABphi[4][4],double ddABphi[4][4][4]);

/* from Coordinates_AnsorgNS.c */
double   x_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double   y_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double   z_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dA_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dA_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dA_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dB_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dB_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dB_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dphi_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dphi_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dphi_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddA_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddB_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddphi_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi);

double   x_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double   y_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double   z_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dA_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dA_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dA_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dB_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dB_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dB_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dphi_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dphi_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dphi_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddA_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddB_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddphi_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi);

double   x_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double   y_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double   z_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dA_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dA_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dA_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dB_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dB_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dB_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dphi_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dphi_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dphi_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddA_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddB_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddphi_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi);

double   x_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double   y_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double   z_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dA_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dA_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dA_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dB_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dB_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dB_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dphi_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dphi_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dphi_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddA_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddB_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddphi_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi);

// IAnsorgNS? functions
double   x_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double   y_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double   z_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dx_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dx_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dx_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dy_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dy_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dy_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dz_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dz_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  dz_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddx_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddy_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);
double  ddz_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi);

double   x_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double   y_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double   z_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dx_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dx_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dx_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dy_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dy_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dy_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dz_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dz_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  dz_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddx_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddy_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);
double  ddz_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi);

double   x_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double   y_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double   z_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dx_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dx_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dx_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dy_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dy_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dy_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dz_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dz_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  dz_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddx_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddy_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);
double  ddz_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi);

double   x_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double   y_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double   z_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dx_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dx_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dx_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dy_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dy_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dy_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dz_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dz_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  dz_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddx_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddy_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);
double  ddz_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi);

/* from coordtrans_CubedSphere.c */
double x_of_CubedSphere(void *aux, int ind, double lam, double A, double B);
double y_of_CubedSphere(void *aux, int ind, double lam, double A, double B);
double z_of_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dlam_dx_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dlam_dy_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dlam_dz_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dA_dx_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dA_dy_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dA_dz_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dB_dx_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dB_dy_CubedSphere(void *aux, int ind, double lam, double A, double B);
double dB_dz_CubedSphere(void *aux, int ind, double lam, double A, double B);
