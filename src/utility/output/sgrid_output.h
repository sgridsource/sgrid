/* sgrid_output.h */
/* Wolfgang Tichy, April 2005 */

    
/* output.c */
int timeforoutput_di_dt(tGrid *grid, int di, double dt);
int timeforoutput_index(tGrid *grid, int index);
int timeforoutput(tGrid *grid, tVarList *vl);
int timeforoutput_any(tGrid *grid);
int find_ind_closest_to_X0(tBox *box, double X0);
int find_ind_closest_to_Y0(tBox *box, double Y0);
int find_ind_closest_to_Z0(tBox *box, double Z0);

int write_grid(tGrid *grid);

/* quick_out.c for debugging */
void quick_Var_output(tBox *box, char *name,
                      double fake_t, int fake_i, int stop);
void quick_Array_output(tBox *box, double *Ar, char *name,
                        double fake_t, int fake_i, int stop);
