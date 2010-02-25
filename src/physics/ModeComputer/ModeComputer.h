/* ModeComputer.h */
/* Wolfgang Tichy 2/2010 */


double Re_sYlm(int l, int m, int s, double costheta, double phi);
double Im_sYlm(int l, int m, int s, double costheta, double phi);
int ModeComputer_read_spheregrid(FILE *fp, int *n2, int *n3, char *gridtype,
                                 int *ntheta, int *nphi);
int ModeComputer_set_boxsize(tGrid *grid);
int ModeComputer_startup(tGrid *grid);
int ModeComputer(tGrid *grid);
