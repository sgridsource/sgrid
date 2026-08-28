/* numerics.h */
/* Wolfgang Tichy, June 2022 */

/* make sure symbol rename is done here as well: */
#include "../../main/main/rename_symbols.h"

/* from main/main/utilities.c */
int finit(double x);


/* rtbrent_brak.c */
int rtbrent_brak_1dVF(double *x0,
                      void (*vecfuncP)(int n,double x[], double f[],void *par),
                      double x1, double x2, void *par, int vecfuncP_ilow,
                      int itmax, double xacc, int pr);
