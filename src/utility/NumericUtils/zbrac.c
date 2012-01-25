#include <stdio.h>
#include <math.h>
#define FACTOR 1.6
#define NTRY 50

/* same as zbrac from numrec, but with pars for func
  returns number of tries j    if ok
  returns <0                   if failure */
int zbrac_P(double (*func)(double,void *par), double *x1, double *x2, void *par)
{
	void nrerror(char error_text[]);
	int j;
	double f1,f2;

	if (*x1 == *x2)
	{
	  printf("Bad initial range in zbracP\n");
	  return -1;
	}
	f1=(*func)(*x1, par);
	f2=(*func)(*x2, par);
	for (j=1;j<=NTRY;j++) {
		if (f1*f2 < 0.0) return j;
		if (fabs(f1) < fabs(f2))
			f1=(*func)(*x1 += FACTOR*(*x1-*x2), par);
		else
			f2=(*func)(*x2 += FACTOR*(*x2-*x1), par);
	}
	printf("Did not find bracket in zbracP\n");
	return -j-1;
}
#undef FACTOR
#undef NTRY
