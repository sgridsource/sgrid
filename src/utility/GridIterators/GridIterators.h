/* GridIterators.h */
/* Wolfgang Tichy 8/2008 */


/* abbreviations used in GridIterators for funcs from utility.c */
#define dot(vlu,vlv)  GridDotProduct((vlu), (vlv))
#define norm2(vlu)    GridL2Norm((vlu))


/* from WTsolver.c */
void copy_varlist_into_array(tVarList *vlx, double *x);
void copy_array_into_varlist(double *x, tVarList *vlx);
int WTiterator(tSparseVector **Aline, int nlines,
               double *x, double *b, int itmax, double tol, double *normres);
