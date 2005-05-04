/* Coordinates.h */
/* Wolfgang Tichy 4/2005 */


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
