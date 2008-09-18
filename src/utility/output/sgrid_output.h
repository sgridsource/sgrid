/* sgrid_output.h */
/* Wolfgang Tichy, April 2005 */

    
/* this value marks points in the output for which no data is available */


/* output.c */
int timeforoutput_di_dt(tGrid *grid, int di, double dt);
int timeforoutput_index(tGrid *grid, int index);
int timeforoutput(tGrid *grid, tVarList *vl);
int timeforoutput_any(tGrid *grid);

int write_grid(tGrid *grid);
