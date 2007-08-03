/* Coordinates.h */
/* Wolfgang Tichy 4/2005 */

/* from doubleCovering.c */
void reset_doubleCoveredPoints_SphericalDF(tBox *box, tVarList *unew);

/* from coordFilters.c */
void coordinateDependentFilter_SphericalDF(tBox *box, tVarList *unew);

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
double x_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double y_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
double z_of_AnsorgNS0(void *aux, int ind, double A, double B, double phi);
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
