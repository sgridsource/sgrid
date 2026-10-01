/* rt_brak.c */
/* (c) Wolfgang Tichy 6.6.2022 */

#include <stdio.h>
#include <math.h>

/* find a bracket [x1,x2] that contains a root of a 1d function func(x, par)
   returns number of tries j    if ok
   returns <0                   if failure */
int rt_brak_f(double (*func)(double x, void *par), void *par,
              double *x1, double *x2, double Fac, int maxits)
{
  int j;
  double f1,f2;

  if (*x1 == *x2)
  {
    printf("rt_brak_f: Warning: *x1 and *x2 are equal!\n");
    return -1;
  }

  f1 = (*func)(*x1, par);
  f2 = (*func)(*x2, par);

  for(j=1; j<=maxits; j++)
  {
    /* return if we have a bracket */
    if(f1*f2 < 0.) return j;

    if(fabs(f1) < fabs(f2))
    {
      /* move bound x1 */
      *x1 += Fac * (*x1 - *x2);
      f1 = (*func)(*x1, par);
    }
    else
    {
      /* move bound x2 */
      *x2 += Fac * (*x2 - *x1);
      f2 = (*func)(*x2, par);
    }
  }
  printf("rt_brak_f: Could not find bracket!\n");
  return -j-1;
}

/* widen bracket with standard factor of 1.6 */
int rt_brak(double (*func)(double x, void *par), void *par,
            double *x1, double *x2, int maxits)
{
  return rt_brak_f(func, par, x1,x2, 1.6, maxits);
}
