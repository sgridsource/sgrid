/*rtsafe_itsP.c,  Wolfgang Tichy 9/2015 */

#include <stdio.h>
#include <math.h>

/* Newton-Raphson with bracketing:

   similar to numrec's rtsafe, but it allows for par and MAXIT 
   args and returns error code. The actual root is returned in x0.
   returns j>=0   if ok,        j = number of iterations done
   returns <0     if error!!! 
   returns root in *x0
   *x0 is used as starting guess if x1<=x0<=x2 || x2<=x0<=x1  */
int rtsafe_itsP(double *x0, 
                void (*funcd)(double x, double *f,double *df, void *par),
                double x1, double x2, void *par, int MAXIT, double xacc)
{
  int j, bisect;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  /* check bracket */
  if(!isfinite(x1))  return -2*MAXIT-1; /* bracket must be finite */
  if(!isfinite(x2))  return -2*MAXIT-2;

  /* check if *x0 is within x1,x2 */
  if( (x1<=*x0 && *x0<=x2) || (x2<=*x0 && *x0<=x1) )  rts=*x0;
  else  rts=0.5*(x1+x2); /* we also get this if *x0 is NAN */

  /* get func values */  
  (*funcd)(x1,&fl,&df, par);
  (*funcd)(x2,&fh,&df, par);

  /* check if root is bracketed and if func values are not NAN */
  j=0;
  if(fl <= 0.0 && fh >= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(fl >= 0.0 && fh <= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(!j)
  {
    /* printf("Root must be bracketed in rtsafe_itsP\n"); */
    return -2*MAXIT-3;
  }
  
  /* early returns */
  if(fl == 0.0) { *x0=x1; return 0;}
  if(fh == 0.0) { *x0=x2; return 0;}

  if(fl < 0.0)
  {
    xl=x1;
    xh=x2;
  }
  else
  {
    xh=x1;
    xl=x2;
  }

  /* set dx ... */
  dxold=fabs(x2-x1);
  dx=dxold;

  /* iterate */
  for (j=1;j<=MAXIT;j++)
  {
    /* eval func at new x=rts */
    (*funcd)(rts,&f,&df, par);
    if(!isfinite(f)) { return -j; }
    if(f == 0.0) { *x0=rts; return j; }
    if(f < 0.0)
      xl=rts;
    else
      xh=rts;

    /* test if Newton step could be possible */
    if( !isfinite(df) )
      bisect = 1;
    else 
      bisect = (fabs(2.0*f) > fabs(dxold*df));

    if(bisect==0) /* try Newton step */
    {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if(temp == rts) { *x0=rts; return j; }
      /* check if Newton step brings us outside interval */
      if( (rts-xh)*(rts-xl) > 0.0 ) { bisect=1; rts=temp; dx=dxold; }
    }
    if(bisect) /* Bisection step */
    {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if(xl == rts) { *x0=rts; return j; }
    }
    /* check accuracy goal */
    if(fabs(dx) < xacc) { *x0=rts; return j; }
  }
  /* printf("Maximum number of iterations exceeded in rtsafe_itsP\n"); */
  *x0=rts;
  return j;
}
