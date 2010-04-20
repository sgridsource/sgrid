/* PN_CircularOrbit_GWs.h */
/* Wolfgang Tichy 2/2010 */

/* funcs in PN_CircularOrbit_GWs.c */
int PN_CircularOrbit_GWs_set_boxsize(tGrid *grid);
int PN_CircularOrbit_GWs_startup(tGrid *grid);
int PN_CircularOrbit_GWs(tGrid *grid);

/* Petr's funcs in compute_graw_wave.c and xodeint.c */
void xodeint(double m1_in, double m2_in, double t1, double t2, double ystart_in[]);
void compute_hcross_hplus(double y[], double *hcross, double *hplus, double D, double theta, double phi, double m1, double m2);
