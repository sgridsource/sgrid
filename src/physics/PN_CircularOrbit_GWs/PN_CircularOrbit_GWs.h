/* PN_CircularOrbit_GWs.h */
/* Wolfgang Tichy 2/2010 */

/* types of EOM we can use */
enum
{
  Kidder1995,       /* as in Kidder, PRD 52, 821 (1995) */
  BuonannoEtAl2003, /* Petr's implementation */
  TaylorT4_bug,     /* M. Boyle: this old T4 from SpEC has a bug  */
  TaylorT4          /* Taylor T4 as used in SpEC */
};


/* funcs in PN_CircularOrbit_GWs.c */
int PN_CircularOrbit_GWs_set_boxsize(tGrid *grid);
int PN_CircularOrbit_GWs_startup(tGrid *grid);
int PN_CircularOrbit_GWs(tGrid *grid);
void compute_hplus_hcross_on_sphere(tBox *box, int hpind, int hxind,
                                    double *yvec, 
                                    double D, double m1, double m2,
                                    int imin, int imax,
                                    int set_doublecovered_points);
void compute_FDpsi4_and_hplus_hcross_on_sphere(tBox *box,
                                               int Rpsi4ind, int Ipsi4ind, 
                                               int hpind, int hxind, double *yin, 
                                               double D, double m1, double m2,
                                               double t, double dt,
                                               int imin, int imax,
                                               int set_doublecovered_points);
void compute_sYlmModes_of_H(tBox *box, int ReHind, int ImHind,
                            int Re_sYlmind, int Im_sYlmind, int lmax,
                            int Re_Hmind, int Im_Hmind, 
                            int ImHmodeSign);
void output_sYlmModes_of_PN_H(char *outname, double time,
                              double *Re_Hmodep, double *Im_Hmodep,
                              int lmax, int s, int ImHmodeSign);


/* Petr's funcs in compute_graw_wave.c and PN_CircOrbit_xodeint.c */
void PN_CircOrbit_compute_constants(double m1_in, double m2_in);
double PN_CircOrbit_compute_r(double y[]);
void PN_CircOrbit_xodeint(double m1_in, double m2_in, double t1, double t2, double ystart_in[]);
void compute_hcross_hplus(double y[], double *hcross, double *hplus, double D, double theta, double phi, double m1, double m2);

/* funcs in minimize_PN-NR_diff.c */
int minimize_PN_NR_diff(tGrid *grid);
