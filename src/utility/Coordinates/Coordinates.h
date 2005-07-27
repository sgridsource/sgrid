/* Coordinates.h */
/* Wolfgang Tichy 4/2005 */

/* from doubleCovering.c */
void reset_doubleCoveredPoints_SphericalDF(tBox *box, tVarList *unew);

/* from coordFilters.c */
void coordinateDependentFilter_SphericalDF(tBox *box, tVarList *unew);

/* Some trivial functions */
double zero_of_xyz(void *aux, double X, double Y, double Z);
double one_of_xyz(void *aux, double X, double Y, double Z);
double x_equals_X(void *aux, double X, double Y, double Z);
double y_equals_Y(void *aux, double X, double Y, double Z);
double z_equals_Z(void *aux, double X, double Y, double Z);

/* Polar coordinates: */
double x_ofPolar(void *aux, double rho, double phi, double Z);
double y_ofPolar(void *aux, double rho, double phi, double Z);
double z_ofPolar(void *aux, double rho, double phi, double Z);
double drho_dx(void *aux, double rho, double phi, double Z);
double drho_dy(void *aux, double rho, double phi, double Z);
double dphi_dx(void *aux, double rho, double phi, double Z);
double dphi_dy(void *aux, double rho, double phi, double Z);
void set_d_dx_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
void set_d_dy_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
/*
double drho_dxdx(void *aux, double rho, double phi, double Z);
double drho_dxdy(void *aux, double rho, double phi, double Z);
double drho_dydy(void *aux, double rho, double phi, double Z);
double dphi_dxdx(void *aux, double rho, double phi, double Z);
double dphi_dxdy(void *aux, double rho, double phi, double Z);
double dphi_dydy(void *aux, double rho, double phi, double Z);
*/

/* PolarCE coordinates: */
double x_ofPolarCE(void *aux, double rho, double Y, double Z);
double y_ofPolarCE(void *aux, double rho, double Y, double Z);
double dYPolarCE_dx(void *aux, double rho, double Y, double Z);
double dYPolarCE_dy(void *aux, double rho, double Y, double Z);

/* SphericalDF coordinates: */
double x_ofSphericalDF(void *aux, double r, double thm, double phi);
double y_ofSphericalDF(void *aux, double r, double thm, double phi);
double z_ofSphericalDF(void *aux, double r, double thm, double phi);
double dr_dx(void *aux, double r, double thm, double phi);
double dr_dy(void *aux, double r, double thm, double phi);
double dr_dz(void *aux, double r, double thm, double phi);
double dphiSphericalDF_dx(void *aux, double r, double thm, double phi);
double dphiSphericalDF_dy(void *aux, double r, double thm, double phi);
double dthm_dx(void *aux, double r, double thm, double phi);
double dthm_dy(void *aux, double r, double thm, double phi);
double dthm_dz(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dxdx(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dxdy(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dxdz(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dydy(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dydz(void *aux, double r, double thm, double phi);
double ddr_SphericalDF_dzdz(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dxdx(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dxdy(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dxdz(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dydy(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dydz(void *aux, double r, double thm, double phi);
double ddthm_SphericalDF_dzdz(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dxdx(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dxdy(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dxdz(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dydy(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dydz(void *aux, double r, double thm, double phi);
double ddphi_SphericalDF_dzdz(void *aux, double r, double thm, double phi);

/* compactSphericalDF coordinates: */
double r_of_xi(double xi);
double dxi_dr(double xi);
double x_ofcompactSphericalDF(void *aux, double xi, double thm, double phi);
double y_ofcompactSphericalDF(void *aux, double xi, double thm, double phi);
double z_ofcompactSphericalDF(void *aux, double xi, double thm, double phi);
double dxi_dx(void *aux, double xi, double thm, double phi);
double dxi_dy(void *aux, double xi, double thm, double phi);
double dxi_dz(void *aux, double xi, double thm, double phi);
double dthmcompactSphericalDF_dx(void *aux, double xi, double thm, double phi);
double dthmcompactSphericalDF_dy(void *aux, double xi, double thm, double phi);
double dthmcompactSphericalDF_dz(void *aux, double xi, double thm, double phi);
double dphicompactSphericalDF_dx(void *aux, double xi, double thm, double phi);
double dphicompactSphericalDF_dy(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdx(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdy(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dxdz(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dydy(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dydz(void *aux, double xi, double thm, double phi);
double ddxi_compactSphericalDF_dzdz(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdx(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdy(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dxdz(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dydy(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dydz(void *aux, double xi, double thm, double phi);
double ddthm_compactSphericalDF_dzdz(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdx(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdy(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dxdz(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dydy(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dydz(void *aux, double xi, double thm, double phi);
double ddphi_compactSphericalDF_dzdz(void *aux, double xi, double thm, double phi);

/* tan_stretch coordinates: */
double x_of_xs(double xs);
double dxs_dx(double xs);
double x_of_tan_stretch(void *aux, double xs, double ys, double zs);
double y_of_tan_stretch(void *aux, double xs, double ys, double zs);
double z_of_tan_stretch(void *aux, double xs, double ys, double zs);
double dxs_dx_tan_stretch(void *aux, double xs, double ys, double zs);
double dys_dy_tan_stretch(void *aux, double xs, double ys, double zs);
double dzs_dz_tan_stretch(void *aux, double xs, double ys, double zs);
