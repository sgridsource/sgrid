/* GridIterators.h */
/* Wolfgang Tichy 8/2008 */


/* abbreviations used in GridIterators for funcs from utility.c */
#define dot(vlu,vlv)  GridDotProduct((vlu), (vlv))
#define norm2(vlu)    GridL2Norm((vlu))

/* init Newton */
int Init_Newton(tGrid *grid);

/* from WTsolver.c */
int WTiterator(tSparseVector **Aline, int nlines,
               double *x, double *b, int itmax, double tol, double *normres);

/* from utility.c */
void copy_varlist_into_array(tVarList *vlx, double *x);
void copy_array_into_varlist(double *x, tVarList *vlx);
void copy_varlist_into_array_LEVEL3(tVarList *vlx, double *x);
void copy_array_into_varlist_LEVEL3(double *x, tVarList *vlx);
void copy_varlist_into_array_forSortedVars(tVarList *vlx, double *x);
void copy_array_into_varlist_forSortedVars(double *x, tVarList *vlx);
void copy_vl_into_array_outervlloop(tVarList *vlx, double *xa);
void copy_array_into_vl_outervlloop(double *xa, tVarList *vlx);
void copy_varlistCompInSubbox_into_array(tVarList *vlx, int vi, int bi,
                                         int sbi, int sbj, int sbk,
                                         int nsb1, int nsb2, int nsb3,
                                         double *x);
void copy_array_into_varlistCompInSubbox(double *x, 
                                         tVarList *vlx, int vi, int bi,
                                         int sbi, int sbj, int sbk,
                                         int nsb1, int nsb2, int nsb3);

/* from wrappers_for_JacobiPrecon.c */
void Jacobi_Preconditioner_from_DiagM(tVarList *vlx, tVarList *vlb,
                                      tVarList *vlc1, tVarList *vlc2);
