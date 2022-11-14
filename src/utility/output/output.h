/* output.h */
/* Wolfgang Tichy 4/2005 */

#include <sys/stat.h>
#include <sys/types.h>

void gnuplot_out2d_boxvar(tBox *box, char *name);
void xgraph_out1_boxvar(tBox *box, char *name);
void dump3d_boxvar(tBox *box, char *name);

void output0d_value(char *filename, double time, double val);
void output0d_boxvar(tBox *box, char *name);

/* from output3d.c and VTK_out.c */
void write3d_boxvar(tBox *box, char *name);
void write_raw_vtk_data(FILE *fp, double *buffer, int n, int stride, int offset,
		        int dbl, int flt, int text, int binary);
void write_raw_vtk_points(FILE *fp, double *px, double *py, double *pz,
                          int n, int stride, int offset,
                        int dbl, int flt, int text, int binary);
FILE *fopen_vtk(char *varname, char *suffix, int b, int n);
