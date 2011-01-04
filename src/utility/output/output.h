/* output.h */
/* Wolfgang Tichy 4/2005 */

void gnuplot_out2d_boxvar(tBox *box, char *name);
void xgraph_out1_boxvar(tBox *box, char *name);
void dump3d_boxvar(tBox *box, char *name);

void output0d_value(char *filename, double time, double val);
void output0d_boxvar(tBox *box, char *name);
