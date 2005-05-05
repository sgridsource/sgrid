/* output.h */
/* Wolfgang Tichy 4/2005 */

void gnuplot_out2d_boxvar(tBox *box, char *name);
void xgraph_out1_boxvar(tBox *box, char *name);

int find_ind_closest_to_X0(tBox *box, double X0);
int find_ind_closest_to_Y0(tBox *box, double Y0);
int find_ind_closest_to_Z0(tBox *box, double Z0);
