/* minbrent_brak.c */
/* WT 6/2022 */

#include <stdio.h>
#include <math.h>

#define A_WITH_SIGN_B(A,B) ((B)>=0. ? fabs(A) : -fabs(A))
#define MOVE4(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/* Golden ration and related things */
const double golden = 0.3819660;  /* golden = (3 - sqrt(5))/2 */
const double omgolden = 0.618034; /* 1 - golden */
const double tmgolden = 1.618034; /* 1 + (1 - golden) */


/* find x=xmin where f(x, p) is minimal using Brent's method.
   returns iter>0    if ok,     iter = number of iterations done
   you need to check if iter>maxits !!!
   returns <0        if error!!! */
#define ZERO_EPS 1.0e-10
int minbrent_brak(double (*f)(double x, void *p), void *par,
                  double ax, double bx, double cx,
                  int maxits, double tol,
                  double *xmin, double *fmin, int pr)
{
  double a,b;
  double fu,fv,fw,fx, u,v,w,x;
  double e;
  int iter;

  /* order a,b */
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);

  /* init */
  x  =  w =  v = bx;
  fx = fw = fv = (*f)(x, par);

  /* do iterations */
  e = 0.;
  for(iter=1; iter<=maxits; iter++)
  {
    double d, etemp,p,q,r;
    double xm = 0.5*(a+b);
    double tol1 = tol*fabs(x) + ZERO_EPS;
    double tol2 = 2.*tol1;

    if(fabs(x-xm) <= (tol2 - 0.5*(b-a)))
    {
      *xmin = x;
      *fmin = fx;
      return iter;
    }

    if(fabs(e) > tol1)
    {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.*(q-r);
      if(q > 0.) p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
      {
        e = (x >= xm ? a-x : b-x);
        d = golden * e;
      }
      else
      {
      	d = p/q;
      	u = x+d;
      	if(u-a < tol2 || b-u < tol2) d = A_WITH_SIGN_B(tol1, xm-x);
      }
    }
    else
    {
      e = (x >= xm ? a-x : b-x);
      d = golden * e;
    }

    u = (fabs(d) >= tol1 ? x+d : x+A_WITH_SIGN_B(tol1, d));
    fu = (*f)(u, par);
    if (fu <= fx)
    {
      if(u >= x) a = x;
      else       b = x;
      MOVE4(v,w,x,u)
      MOVE4(fv,fw,fx,fu)
    }
    else
    {
      if(u < x) a = u;
      else      b = u;

      if (fu <= fw || w == x)
      {
      	v = w;
      	w = u;
      	fv = fw;
      	fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
      	v = u;
      	fv = fu;
      }
    }
  }
  if(pr) printf("minbrent_brak: Too many iterations!");
  *xmin = x;
  *fmin = fx;
  return iter;
}
#undef ZERO_EPS


/* find a bracket for a min */
#define GLIMIT 100.0
#define TINY 1.0e-20
int min_brak(double (*func)(double x, void *p), void *par,
             double *ax, double *bx, double *cx,
             double *fa, double *fb, double *fc, int maxits)
{
  double ulim,u,r,q,fu,dum;
  int iter;

  *fa = (*func)(*ax, par);
  *fb = (*func)(*bx, par);
  if(*fb > *fa)
  {
    MOVE4(dum,*ax,*bx,dum)
    MOVE4(dum,*fb,*fa,dum)
  }

  *cx = (*bx) + tmgolden*(*bx-*ax);
  *fc = (*func)(*cx, par);

  for(iter=1; iter<=maxits; iter++)
  {
    if(*fb <= *fc) break;

    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);

    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r) /
              (2.*A_WITH_SIGN_B(fmax(fabs(q-r),TINY), q-r));

    ulim = (*bx) + GLIMIT*(*cx-*bx);

    if((*bx-u)*(u-*cx) > 0.)
    {
      fu = (*func)(u, par);
      if(fu < *fc)
      {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
        return iter;
      }
      else if (fu > *fb)
      {
        *cx = u;
        *fc = fu;
        return iter;
      }
      u = (*cx) + tmgolden*(*cx-*bx);
      fu = (*func)(u, par);
    }
    else if((*cx-u)*(u-ulim) > 0.)
    {
      fu = (*func)(u, par);
      if(fu < *fc)
      {
        MOVE4(*bx, *cx,  u, *cx + tmgolden*(*cx-*bx))
        MOVE4(*fb, *fc, fu, (*func)(u, par))
      }
    }
    else if((u-ulim)*(ulim-*cx) >= 0.)
    {
      u = ulim;
      fu = (*func)(u, par);
    }
    else
    {
      u = (*cx) + tmgolden*(*cx-*bx);
      fu = (*func)(u, par);
    }
    MOVE4(*ax, *bx, *cx, u)
    MOVE4(*fa, *fb, *fc, fu)
  }
  return iter;
}
#undef GLIMIT
#undef TINY
#undef MOVE4
#undef A_WITH_SIGN_B
