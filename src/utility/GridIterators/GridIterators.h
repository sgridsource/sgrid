/* GridIterators.h */
/* Wolfgang Tichy 8/2008 */


/* from utility.c */
double dot(tVarList *u, tVarList *v);
double norm2(tVarList *u);

/* from WTsolver.c */
void copy_varlist_into_array(tVarList *vlx, double *x);
void copy_array_into_varlist(double *x, tVarList *vlx);
int WTiterator(tSparseVector **Aline, int nlines,
               double *x, double *b, int itmax, double tol, double *normres);
