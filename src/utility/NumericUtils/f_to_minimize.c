/* (c) Wolfgang Tichy 6.6.2003 */
/* adapted from numrec fmin.c */

#define NRANSI
#include "nrutil.h"

extern int newton_lnsrch_nn;
extern double *newton_lnsrch_fvec;
extern void (*newton_lnsrch_nrfuncv)(int n, double v[], double f[]);

double f_to_minimize(double x[])
{
	int i;
	double sum;

	(*newton_lnsrch_nrfuncv)(newton_lnsrch_nn,x,newton_lnsrch_fvec);
	for (sum=0.0,i=1;i<=newton_lnsrch_nn;i++) 
	  sum += SQR(newton_lnsrch_fvec[i]);
	return 0.5*sum;
}
#undef NRANSI
