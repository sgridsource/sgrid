/* PN_CircularOrbit_GWs.h */
/* Wolfgang Tichy 2/2010 */

/* funcs in PN_CircularOrbit_GWs.c */
int PN_CircularOrbit_GWs_set_boxsize(tGrid *grid);
int PN_CircularOrbit_GWs_startup(tGrid *grid);
int PN_CircularOrbit_GWs(tGrid *grid);
void compute_hplus_hcross_on_sphere(tBox *box, int hpind, int hxind,
                                    double *yvec, 
                                    double D, double m1, double m2,
                                    int imin, int imax,
                                    int set_doublecovered_points);
void compute_psi4_and_hplus_hcross(tBox *box, int Rpsi4ind, int Ipsi4ind, 
                                   int hpind, int hxind, double *yin, 
                                   double D, double m1, double m2,
                                   double t, double dt,
                                   int imin, int imax,
                                   int set_doublecovered_points);
void compute_sYlmModes_of_PN_H(tBox *box, int ReHind, int ImHind,
                               int Re_sYlmind, int Im_sYlmind, int lmax,
                               int Re_Hmind, int Im_Hmind, 
                               int ImHmodeSign);


/* Petr's funcs in compute_graw_wave.c and xodeint.c */
void xodeint(double m1_in, double m2_in, double t1, double t2, double ystart_in[]);
void compute_hcross_hplus(double y[], double *hcross, double *hplus, double D, double theta, double phi, double m1, double m2);
