/* (c) Wolfgang Tichy 6.6.2003 */
/* adapted from numrec lnsrch.c, tolerances changed */

#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ALF 1.0e-4
#define TOLX 1.0e-14

void lnsrch_double(int n, double xold[], double fold, 
	double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(double []))
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) nrerror("Roundoff problem in lnsrch_double.");
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

/* same as above but with 3 more args, so we can use f_to_minimize_in_linesrch */
void linesrch(int n, double xold[], double fold, 
	double g[], double p[], double x[],
	double *f, double stpmax, int *check, 
	double (*func)(double [], int, double *newton_linesrch_fvec,
	               void (*vecfunc)(int, double [], double [])),
	double *newton_linesrch_fvec, 
	void (*vecfunc)(int, double [], double []))
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x, n, newton_linesrch_fvec, vecfunc);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) nrerror("Roundoff problem in lnsrch_double.");
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

/* same as linesrch but pass on void *par */
void linesrchP(int n, double xold[], double fold, 
	double g[], double p[], double x[],
	double *f, double stpmax, int *check, 
	double (*func)(double [], int, double *newton_linesrch_fvec,
	               void (*vecfunc)(int, double [], double [], void *par),
	               void *par),
	double *newton_linesrch_fvec, 
	void (*vecfunc)(int, double [], double [], void *par),
	void *par)
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x, n, newton_linesrch_fvec, vecfunc, par);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) nrerror("Roundoff problem in linesrchP.");
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

#undef ALF
#undef TOLX
#undef NRANSI
