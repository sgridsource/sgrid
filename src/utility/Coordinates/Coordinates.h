/* Coordinates.h */
/* Wolfgang Tichy 4/2005 */


/* Some trivial functions */
double zero_of_xyz(double X, double Y, double Z);
double one_of_xyz(double X, double Y, double Z);
double x_equals_X(double X, double Y, double Z);
double y_equals_Y(double X, double Y, double Z);
double z_equals_Z(double X, double Y, double Z);

/* Polar coordinates: */
double x_ofPolar(double rho, double phi, double Z);
double y_ofPolar(double rho, double phi, double Z);
double z_ofPolar(double rho, double phi, double Z);
double drho_dx(double x, double y, double z);
double drho_dy(double x, double y, double z);
double dphi_dx(double x, double y, double z);
double dphi_dy(double x, double y, double z);
void set_d_dx_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
void set_d_dy_at_rhoEQzero(void *bo, void *va, void *v1,void *v2,void *v3);
/*
double drho_dxdx(double x, double y, double z);
double drho_dxdy(double x, double y, double z);
double drho_dydy(double x, double y, double z);
double dphi_dxdx(double x, double y, double z);
double dphi_dxdy(double x, double y, double z);
double dphi_dydy(double x, double y, double z);
*/
