/* newton_linesearch.c */
/* (c) Wolfgang Tichy 6.6.2022 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numerics.h"

#define SQR(a) (a)*(a)


/* external funcs newton_linesearch needs: */
int lu_decomp(int n, double A[][n], int *p_idx, int *parity);
void lu_solve(int n, double A[][n], int *p_idx, double b[]);


/* helper funcs for newton_linesearch in this file: */
void fd_jacobianP(int n, double x[],
        double newton_linesrch_fvec[], double(*df)[n],
        void (*vecfuncP)(int, double [], double [], void *par),
        void *par, int vecfuncP_ilow);
double f_to_minimize_in_linesrchP(int n, double x[],
        double *newton_linesrch_fvec,
        void (*vecfuncP)(int, double [], double [], void *par),
        void *par, int vecfuncP_ilow);
void linesrchP(int n, double xold[], double fold, double g[], double p[],
        double x[], double *f, double stpmax, int *check, 
        double (*func)(int, double [], double *newton_linesrch_fvec,
                       void (*vecfuncP)(int, double [], double [], void *par), 
                       void *par, int vecfuncP_ilow),
        double *newton_linesrch_fvec, 
        void (*vecfuncP)(int, double [], double [], void *par),
        void *par, int vecfuncP_ilow);


/* Newton linesearch with vecfuncP that has an argument par for
   parameters (a pointer to void):
   returns its>=0   if ok, 
   returns -its-1<0 if error!!!    its = number of iterations done
   returns -its-itmax*10<0  if f=NAN.
   here x[i0...n-1+i0] with i0 = vecfuncP_ilow (normally vecfuncP_ilow=0)
   The function vecfuncP has to be of the form:
    void vecfuncP(int n, double *x, double *fvec, void *par)
                  where x[i0...n-1+i0], fvec[i0...n-1+i0],
                  with i0 = vecfuncP_ilow, 
                  par can point to something that contains parameters.
   Then newton_linesearch will find x s.t. vecfuncP = 0.               */
#define TOLMIN 1.0e-13
#define TOLX 1.0e-14
#define STPMX 100.0

#define FREE_RETURN {free(newton_linesrch_fvec); free(xold);\
                     free(p); free(g); free(fjac); free(p_idx);\
                     if(!finit(f)) return -its-itmax*10;\
                     return its;}
#define FREE_RETURN_ERR {free(newton_linesrch_fvec); free(xold);\
                         free(p); free(g); free(fjac); free(p_idx);\
                         if(!finit(f)) return -its-itmax*10;\
                         return -its-1;}
int newton_linesearch(int n, double x[], int *check,
                      void (*vecfuncP)(int n, double x[], double f[],
                                       void *par),
                      void *par, int vecfuncP_ilow, int itmax, double tolf)
{
  int i,j;
  int its=0; /* no iterations yet */
  double f,fold, stpmax, sum,temp,test;
  int *p_idx;
  double (*fjac)[n];
  double *g, *p, *xold, *newton_linesrch_fvec;

  /* shift x pointer in case the first enrty in the x that was passed in
     is at index vecfuncP_ilow */
  x = x + vecfuncP_ilow;

  p_idx = malloc(n * sizeof(int));
  fjac  = malloc(n*n * sizeof(double));
  g     = malloc(n * sizeof(double));
  p     = malloc(n * sizeof(double));
  xold  = malloc(n * sizeof(double));
  newton_linesrch_fvec = malloc(n * sizeof(double));
  if(!p_idx || !fjac || !g || !p || !xold || !newton_linesrch_fvec)
  {
    printf("allocation failure!!!");
    abort();
    exit(1);
  }

  f=f_to_minimize_in_linesrchP(n, x, newton_linesrch_fvec,
                               vecfuncP, par, vecfuncP_ilow);

  test = 0.0;
  for(i=0; i<n; i++)
    if(fabs(newton_linesrch_fvec[i]) > test)
      test = fabs(newton_linesrch_fvec[i]);
  if(test < 0.01*tolf)
  {
    *check=0;
    FREE_RETURN
  }

  for(sum=0.0,i=0; i<n; i++) sum += SQR(x[i]);
  stpmax = STPMX*fmax(sqrt(sum),(double)n);
  for(its=1; its<=itmax; its++)
  {
    int parity;

    fd_jacobianP(n,x,newton_linesrch_fvec,fjac,vecfuncP, par, vecfuncP_ilow);
    for(i=0; i<n; i++)
    {
      for(sum=0.0,j=0; j<n; j++) sum += fjac[j][i]*newton_linesrch_fvec[j];
      g[i]=sum;
    }
    for(i=0; i<n; i++) xold[i]=x[i];
    fold = f;
    for(i=0; i<n; i++) p[i] = -newton_linesrch_fvec[i];
    /* now p contains RHS of fjac[i][j] dx[i] = -fvec[i] = p[i] */
    j = lu_decomp(n, fjac, p_idx, &parity);
    if(j) /* matrix fjac has all zeros in col j (or row -j) */
    { /* throw away all off-diag elements if fjac is singlar */
      /*
      int k;
      for(k=0; k<n; k++)
      {
        double fjac_kk=fjac[k][k];
        if(fjac_kk!=0.0)  p[k] = p[k]/fjac_kk;
        else              p[k] = 0.0;
      }
      */
      int k,l;
      double ma;
      int *idx = calloc(n, sizeof(int));
      if(!idx) { printf("allocation failure!!!"); abort(); exit(1); }
      for(l=0; l<n; l++)
      {
        for(ma=0.0, k=0; k<n; k++)
        {
          int li, fl;
          double ffkl=fabs(fjac[k][l]);
          if(ffkl>ma)
          {
            for(fl=1, li=0; li<l; li++)
              if(idx[li]==k) { fl=0; break; }
            if(fl) // if k is not in idx yet, mark max & put k into idx
            {
              ma = ffkl;
              idx[l] = k;
            }
          }
        }
        // now we have max ma in fjac[k][l], where k = idx[l]
        k = idx[l];
        if(k>0) p[l] = p[l]/fjac[k][l];
        else    p[l] = 0.0;
      } // end l-loop
      free(idx);
    } /* end Sing. treatment */
    else
    {
      lu_solve(n, fjac, p_idx, p);
    }
    
    linesrchP(n,xold,fold,g,p,x,&f,stpmax,check,
              f_to_minimize_in_linesrchP, 
              newton_linesrch_fvec, vecfuncP, par, vecfuncP_ilow);
    test = 0.0;
    for(i=0; i<n; i++)
      if(fabs(newton_linesrch_fvec[i]) > test)
        test=fabs(newton_linesrch_fvec[i]);
    if(test < tolf)
    {
      *check=0;
      FREE_RETURN
    }
    if(*check)
    {
      double denom = fmax(f,0.5*n);
      test = 0.0;
      for(i=0; i<n; i++)
      {
        temp = fabs(g[i])*fmax(fabs(x[i]),1.0)/denom;
        if(temp > test) test = temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      FREE_RETURN
    }
    test = 0.0;
    for(i=0; i<n; i++)
    {
      temp=(fabs(x[i]-xold[i]))/fmax(fabs(x[i]),1.0);
      if(temp > test) test = temp;
    }
    if(test < TOLX) FREE_RETURN
  }
  /* printf("itmax exceeded in newton_linesearch\n"); */
  FREE_RETURN_ERR
}
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREE_RETURN
#undef FREE_RETURN_ERR



/*********************************************************************/
/* funcs needed by newton_linesearch */
/*********************************************************************/

/* compute Jacobian with finite differences */
#define EPS 1.0e-7    /* should be approx square root of machine precision */
#define HMIN 1.0e-10  /* min h we use in fin. diff. computation of derivs */
void fd_jacobianP(int n, double x[], double fvec[], double(*df)[n],
                  void (*vecfunc)(int, double [], double [], void *par),
                  void *par, int vecfunc_ilow)
{
  int i0=vecfunc_ilow;
  int i,j, k, df_nonzero;
  double h,tmp,*f;

  f = malloc(n * sizeof(double));
  if(!f)
  {
    printf("allocation failure!!!");
    abort();
    exit(1);
  }

  for(j=0; j<n; j++)
  {
    tmp = x[j];
    h = EPS*fabs(tmp);
    if(h < HMIN) h = HMIN; /* make sure h is never below HMIN */
    for(k=0; k<8; k++, h=h*10) /* increase h at most by 8 orders of mag. */
    {
      x[j] = tmp+h;
      h    = x[j]-tmp;
      (*vecfunc)(n, x-i0, f-i0, par);
      x[j] = tmp;
      for(df_nonzero=0, i=0; i<n; i++)
      {
        df[i][j]=(f[i]-fvec[i])/h;
        if(df[i][j] != 0.0) df_nonzero=1;
      }
      if(df_nonzero) break;
    }
  }
  free(f);
}
#undef HMIN
#undef EPS


/* func we minimize in linesrchP */
double f_to_minimize_in_linesrchP(int n, double x[],
        double *newton_linesrch_fvec,
        void (*vecfunc)(int, double [], double [], void *par),
        void *par, int vecfunc_ilow)
{
  int i0=vecfunc_ilow;
  int i;
  double sum;

  (*vecfunc)(n, x-i0, newton_linesrch_fvec-i0, par);
  for(sum=0.0,i=0;i<n;i++)
    sum += SQR(newton_linesrch_fvec[i]);
  return 0.5*sum;
}


#define ALPHA 1.0e-4
#define TOLX 1.0e-14
/* line-search with arg for pars */
void linesrchP(int n, double xold[], double fold, 
        double g[], double p[], double x[],
        double *f, double stpmax, int *check, 
        double (*func)(int, double [], double *newton_linesrch_fvec,
                       void (*vecfunc)(int, double [], double [], void *par),
                       void *par, int vecfunc_ilow),
        double *newton_linesrch_fvec, 
        void (*vecfunc)(int, double [], double [], void *par),
        void *par, int vecfunc_ilow)
{
  int i;
  double a, alam,alam2,alamin, b, discrimi, term1,term2;
  double f2, fold2, g_slope, sum, tmp, test, tmplam;

  *check=0;
  for(sum=0.,i=0;i<n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++)
      p[i] *= stpmax/sum;

  for(g_slope=0.,i=0;i<n;i++)
    g_slope += g[i]*p[i];

  test=0.;
  for(i=0;i<n;i++)
  {
    tmp = fabs(p[i])/fmax(fabs(xold[i]),1.);
    if(tmp > test) test = tmp;
  }

  alamin = TOLX/test;
  alam   = 1.;
  for(;;)
  {
    for(i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(n, x, newton_linesrch_fvec, vecfunc, par, vecfunc_ilow);
    if(alam < alamin)
    {
      for(i=0;i<n;i++) x[i] = xold[i];
      *check = 1;
      return;
    } 
    else if(*f <= fold+ALPHA*alam*g_slope)
    {
      return;
    }
    else
    {
      if(alam == 1.)
      {
        tmplam = -g_slope/(2.*(*f-fold-g_slope));
      }
      else
      {
        term1 = *f-fold-alam*g_slope;
        term2 = f2-fold2-alam2*g_slope;
        a = (term1/(alam*alam)-term2/(alam2*alam2))/(alam-alam2);
        b = (-alam2*term1/(alam*alam)+alam*term2/(alam2*alam2))/(alam-alam2);
        if(a == 0.)
        {
          tmplam = -g_slope/(2.*b);
        }
        else
        {
          discrimi = b*b-3.*a*g_slope;
          if(discrimi<0.)
          { printf("Roundoff problem in linesrchP\n"); abort(); exit(1); }
          else
          { tmplam = (-b+sqrt(discrimi))/(3.*a); }
        }
        if(tmplam>0.5*alam) tmplam = 0.5*alam;
      }
    }
    alam2 = alam;
    f2    = *f;
    fold2 = fold;
    alam  = fmax(tmplam,0.1*alam);
  }
}
#undef ALPHA
#undef TOLX



/***********************************************************/
/* functions to determine if we are within a certain range */
/***********************************************************/

/* check if vec is in a certain range: lo[i] <= vec[i] <= hi[i] */
int newton_lnsrch_vec_within_range(int n, double vec[],
                                   double hi[], double lo[])
{
  int i;
  int ret=1;
  
  for(i=0; i<n; i++)
  {
    if(vec[i]>hi[i]) ret=0;
    if(vec[i]<lo[i]) ret=0;
  }
  return ret;
}

/* use vector vec and range lo[i] <= vec[i] <= hi[i] to obtain vector 
   vb on boundary and vi slightly inside this range */
void newton_lnsrch_set_vecs_for_lininterp(int n, double vec[], 
                double hi[], double lo[], double *vb, double *vi) 
{
  int i;
  double *cent;
  double *dir;
  double s, sc;

  /* alloc */
  cent= malloc(n * sizeof(double));
  dir = malloc(n * sizeof(double));
  if(!cent || !dir)
  {
    printf("allocation failure!!!");
    abort();
    exit(1);
  }

  s=1;
  for(i=0; i<n; i++) 
  {
    /* pos. vec. of center */
    cent[i] = 0.5*(hi[i]+lo[i]);

    /* set directional vec dir */
    dir[i] = vec[i] - cent[i];
    
    /* find min s such that vb[i] = cent[i] + dir[i]*s is on boundary */
    if(dir[i]!=0)
    {
      if(vec[i]>hi[i])  sc = (hi[i] - cent[i])/dir[i];
      if(vec[i]<lo[i])  sc = (lo[i] - cent[i])/dir[i];
      if(sc<s) s=sc;
    }
  }

  /* set vb on boundary and vi a little bit inside boundary */
  for(i=0; i<n; i++)
  {
    vb[i] = cent[i] + dir[i]*s;
    vi[i] = cent[i] + dir[i]*s*0.999999;
    if(vb[i]>hi[i]) vb[i]=hi[i]; /* make sure we stay in range */
    if(vb[i]<lo[i]) vb[i]=lo[i];
  }
  /*
  printf("newton_lnsrch_set_vecs_for_lininterp:\n");
  printf("vec=( ");
  for(i=0; i<n; i++) printf("%g ", vec[i]);
  printf(")\n");
  printf(" vb=( ");
  for(i=0; i<n; i++) printf("%g ", vb[i]);
  printf(")\n");
  printf(" vi=( ");
  for(i=0; i<n; i++) printf("%g ", vi[i]);
  printf(")\n");
  */
  /* free */
  free(cent);
  free(dir);
}

/* use function values fvb and fvi at vb and vi to find function value
   fvec at vec by linear interolation */
void newton_lnsrch_get_fvec_by_lininterp(int n, double vec[], 
          double vb[], double vi[], double *fvec, double fvb[], double fvi[]) 
{
  int i;
  double h, s;

  /* get length of vb[i] - vi[i] and vec[i] - vb[i]*/
  h=s=0;
  for(i=0; i<n; i++)
  {
    double dum;
    
    dum = vb[i] - vi[i];
    h += dum*dum;
    
    dum = vec[i] - vb[i];
    s += dum*dum; 
  }
  h=sqrt(h);
  s=sqrt(s);

  /* use lengths to find fvec from fvb and fvi by lin. interp. */
  for(i=0; i<n; i++) fvec[i] = fvb[i] + ((fvb[i]-fvi[i])/h)*s;
  /*
  printf("newton_lnsrch_get_fvec_by_lininterp:\n");
  printf("fvec=( ");
  for(i=0; i<n; i++) printf("%g ", fvec[i]);
  printf(")\n");
  printf(" fvb=( ");
  for(i=0; i<n; i++) printf("%g ", fvb[i]);
  printf(")\n");
  printf(" fvi=( ");
  for(i=0; i<n; i++) printf("%g ", fvi[i]);
  printf(")\n");
  */
}
