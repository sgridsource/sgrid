/* newton1d_brak.c,  Wolfgang Tichy 6/2019 */

#include <stdio.h>
#include <math.h>

//#include "nmesh.h"
#include "numerics.h"


/* Newton-Raphson with bracketing in [x1,x2]:

   It is similar to sgrid's rtsafe_itsP
   It allows for par and maxits args and returns error code.
   fdf has same signature as fdf in GSL
   The actual root is returned in x0.
   returns j>=0   if ok,        j = number of iterations done
   returns <0     if error!!!
   returns root in *x0
   *x0 is used as starting guess if x1<=x0<=x2 || x2<=x0<=x1 */
int newton1d_brak(double *x0,
                  void (*fdf)(double x, void *par, double *f, double *df),
                  double x1, double x2, void *par, int maxits, double xacc,
                  int pr)
{
  int j, bisect;
  double f,fh,fl, df;
  double xh,xl, xrt, dx,dxold;
  double tmp;

  /* check bracket */
  if( (!finit(x1)) || (!finit(x2)) )
  {
    if(pr) printf("newton1d_brak: Bracket is not finite!  "
                  "x1=%g x2=%g\n", x1,x2);
    return -2*maxits-2;
  }

  /* check if *x0 is within x1,x2 */
  if( (x1<=*x0 && *x0<=x2) || (x2<=*x0 && *x0<=x1) )  xrt=*x0;
  else  xrt=0.5*(x1+x2); /* we also get this if *x0 is NAN */

  /* get func values */
  (*fdf)(x1, par, &fl, &df);
  (*fdf)(x2, par, &fh, &df);

  /* check if root is bracketed and if func values are not NAN */
  j=0;
  if(fl <= 0.0 && fh >= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(fl >= 0.0 && fh <= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(!j)
  {
    if(pr) printf("newton1d_brak: Root is not bracketed!  "
                  "fl=%g fh=%g\n", fl,fh);
    return -2*maxits-3;
  }

  /* early returns */
  if(fl == 0.0) { *x0=x1; return 0; }
  if(fh == 0.0) { *x0=x2; return 0; }

  /* set xl and xh */
  if(fl < 0.0) { xl=x1; xh=x2; }
  else         { xh=x1; xl=x2; }

  /* set dx ... */
  dxold=fabs(x2-x1);
  dx=dxold;

  /* iterate */
  for(j=1; j<=maxits; j++)
  {
    /* eval func at new x=xrt */
    (*fdf)(xrt, par, &f, &df);

    /* catch special cases */
    if(!finit(f)) { return -j; }
    if(f == 0.0) { *x0=xrt; return j; }

    if(f < 0.0) xl=xrt;
    else        xh=xrt;

    /* test if Newton step could be possible */
    if( !finit(df) ) bisect = 1;
    else                bisect = (fabs(2.0*f) > fabs(dxold*df));
    /* note: if df=0 and f!=0 the "else" case sets bisect=1 */

    if(bisect==0) /* try Newton step */
    {
      dxold=dx;
      dx=f/df;
      tmp=xrt;
      xrt -= dx;
      //if(tmp == xrt) { *x0=xrt; return j; } // if df is large this can fail
      /* check if Newton step brings us outside interval,
         or if step was so small that it didn't change xrt */
      if( ((xrt-xh)*(xrt-xl) > 0.0) || (tmp == xrt) )
      {
        bisect=1; xrt=tmp; dx=dxold;
      }
    }
    if(bisect) /* Bisection step */
    {
      dxold=dx;
      dx=0.5*(xh-xl);
      xrt=xl+dx;
      if(xl == xrt) { *x0=xrt; return j; }
    }
    /* check accuracy goal */
    if(fabs(dx) < xacc) { *x0=xrt; return j; }
  }

  if(pr) printf("newton1d_brak: Maximum number of iterations exceeded!  "
                "j=%d > maxits=%d\n", j, maxits);
  *x0=xrt;
  return j;
}


/* Bisection method with bracket [x1,x2]:

   The actual root is returned in x0.
   returns j>=0   if ok,        j = number of iterations done
   returns <0     if error!!!
   returns root in *x0
   *x0 is used as starting guess if x1<=x0<=x2 || x2<=x0<=x1 */
int rtbisect(double *x0, double (*func)(double,void *par),
             double x1, double x2, void *par, int maxits, double xacc, int pr)
{
  int j;
  double f,fh,fl, xh,xl, xrt, dx,dxold;

  /* check bracket */
  if( (!finit(x1)) || (!finit(x2)) )
  {
    if(pr) printf("rtbisect: Bracket is not finite!  "
                  "x1=%g x2=%g\n", x1,x2);
    return -2*maxits-2;
  }

  /* check if *x0 is within x1,x2 */
  if( (x1<=*x0 && *x0<=x2) || (x2<=*x0 && *x0<=x1) )  xrt=*x0;
  else  xrt=0.5*(x1+x2); /* we also get this if *x0 is NAN */

  /* get func values */
  fl = (*func)(x1, par);
  fh = (*func)(x2, par);

  /* check if root is bracketed and if func values are not NAN */
  j=0;
  if(fl <= 0.0 && fh >= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(fl >= 0.0 && fh <= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(!j)
  {
    if(pr) printf("rtbisect: Root is not bracketed!  "
                  "fl=%g fh=%g\n", fl,fh);
    return -2*maxits-3;
  }

  /* early returns */
  if(fl == 0.0) { *x0=x1; return 0; }
  if(fh == 0.0) { *x0=x2; return 0; }

  /* set xl and xh */
  if(fl < 0.0) { xl=x1; xh=x2; }
  else         { xh=x1; xl=x2; }

  /* set dx ... */
  dxold=fabs(x2-x1);
  dx=dxold;

  /* iterate */
  for(j=1; j<=maxits; j++)
  {
    /* eval func at new x=xrt */
    f= (*func)(xrt, par);

    /* catch special cases */
    if(!finit(f)) { return -j; }
    if(f == 0.0) { *x0=xrt; return j; }

    if(f < 0.0) xl=xrt;
    else        xh=xrt;

    /* Bisection step */
    dxold = dx;
    dx = 0.5*(xh - xl);
    xrt = xl + dx;
    if(xl == xrt) { *x0=xrt; return j; }
    /* check accuracy goal */
    if(fabs(dx) < xacc) { *x0=xrt; return j; }
  }

  if(pr) printf("rtbisect: Maximum number of iterations exceeded!  "
                "j=%d > maxits=%d\n", j, maxits);
  *x0=xrt;
  return j;
}
