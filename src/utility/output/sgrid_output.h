/* sgrid_output.h */
/* Wolfgang Tichy, April 2005 */

    
/* this value marks points in the output for which no data is available */


/* output.c */
int timeforoutput(tGrid *grid, tVarList *vl);
int timeforoutput_any(tGrid *grid);

int write_grid(tGrid *grid);
