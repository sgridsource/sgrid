/* (c) Wolfgang Tichy 11.5.2013 */
/* adapted from numrec newt.c + tolerances changed */

#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"




/* the function f = F.F/2 */
double f_to_minimize_in_WT_linesrch(double x[], int n, 
	                            double *newton_linesrch_fvec,
	                            void (*vecfunc)(int, double [], double [], void *par),
	                            void *par)
{
	int i;
	double sum;

	(*vecfunc)(n,x,newton_linesrch_fvec, par);
	for (sum=0.0,i=1;i<=n;i++) 
	  sum += SQR(newton_linesrch_fvec[i]);
	return 0.5*sum;
}

/* line search but with g eliminated */
#define ALF 1.0e-4
#define TOLX 1.0e-14
void WT_linesrch(int n, double xold[], double fold, 
  double p[], double x[],
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

  /* Note: g[j] = F_k J_kj 
     => slope = g[j]*p[j] = d_j f * p[j] = F_i J_ij p[j] = -F_i F_i = -2fold
     before we scale p[j] */
  slope = -2.0*fold;
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  /* scale p[i] and thus slope? */
  if (sum > stpmax)
  {
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
    slope *= stpmax/sum;
  }
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
          if (disc<0.0) nrerror("Roundoff problem in WT_linesrch.");
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



#define TOLMIN 1.0e-13
#define TOLX 1.0e-14
#define STPMX 100.0

/* similar to newton_linesrch_itsP but now we also hand in a linear solver
   and function that gives the result of Jacobian times vector
   returns its>=0   if ok, 
   returns -its-1<0 if error!!!    its = number of iterations done
   returns -its-MAXITS*10<0  if f=NAN.
   The funtion F_x has to be of the form:
    void F_x(int n, double *x, double *fvec, void *par)
                  here x[1...n], fvec[1...n], par can point to something
                  that contains parameters.
   Then WT_newton will find x s.t. F_x = 0.               */
#undef FREERETURN
#define FREERETURN {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return its;}
#define FREERETURNERROR {free_vector(newton_linesrch_fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);\
	if(!finite(f)) return -its-MAXITS*10;\
	return -its-1;}
int WT_newton(double *x, int n, int *check,
        void (*F_x)(int, double *x, double *Fx, void *par),
        void (*J_x_dx)(int, double *dx, double *Jdx, double *x, void *par),
        void *par, int MAXITS, double TOLF,
        int (*linSol)(int n, double *b, double *dx,
            void (*Jdx)(int, double *, double *, double *, void *),
            int (*precon)(int n, double *b, double *dx, double *x, void *par),
            double *x, void *par, int itmax, double tol),
        int (*precon)(int n, double *b, double *dx, double *x, void *par),
        int linitmax, double lintolfac)
{
  int i,j;
  int its=0;
  double lintol,f,fold,stpmax,sum,temp,test,*p,*xold; // ,den
  double *newton_linesrch_fvec;

  p=vector(1,n);
  xold=vector(1,n);
  newton_linesrch_fvec=vector(1,n);
  f=f_to_minimize_in_WT_linesrch(x, n, newton_linesrch_fvec, F_x, par);
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
  if (test < 0.01*TOLF) {
    *check=0;
    FREERETURN
  }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);

  /* main iteration loop */
  for (its=1;its<=MAXITS;its++)
  {
    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;


    /* get tol for linear solver using the precon */
    /* precon(n, newton_linesrch_fvec, p, x, par);
       test=0.0;
       for(i=1;i<=n;i++) if(p[i]>test) test=fabs(p[i]);
       lintol = test*lintolfac; */
    /* get tol for linear solver using newton_linesrch_fvec */
    test=0.0;
    for(i=1;i<=n;i++) 
      if(newton_linesrch_fvec[i]>test) test=fabs(newton_linesrch_fvec[i]);
    lintol = test*lintolfac;

    /* find p from linear solver */
    /* use newton_linesrch_fvec[i] as RHS instead of -newton_linesrch_fvec[i]
       this gives -p[i], so afterwards we flip the sign. */
    j=linSol(n, newton_linesrch_fvec, p, 
             J_x_dx, precon, x, par, linitmax, lintol);
    for(i=1;i<=n;i++) p[i] = -p[i];
    if(j<0) FREERETURNERROR

    /* do line search */
    WT_linesrch(n,xold,fold, p,x,&f,stpmax,check,
                f_to_minimize_in_WT_linesrch, 
                newton_linesrch_fvec, F_x, par);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(newton_linesrch_fvec[i]) > test) test=fabs(newton_linesrch_fvec[i]);
    if (test < TOLF) {
      *check=0;
      FREERETURN
    }
    if (*check) {
      test=0.0;
      //den=FMAX(f,0.5*n);
      for (i=1;i<=n;i++)
      {
        /* we do not have g[i], so set temp=0 for now... */
        //temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
        temp=0.0; /* <- this always results in check=1 */
        if (temp > test) test=temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      FREERETURN
    }
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) FREERETURN
  }
  /* nrerror("MAXITS exceeded in newt"); */
  /* FREERETURNERROR */
  FREERETURN
}
#undef FREERETURNERROR


#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI

/*********************************************************/

/* some test code: */
/*
void F_x(int n, double *x, double *Fx, void *par)
{
  int i;
  for(i=1;i<=n;i++)
    Fx[i] = pow(x[i],i) - i;
  Fx[1] = pow(x[1],3) - 4*x[1]*x[1] + 5;
  for(i=1;i<=n;i++) printf(" i=%d x=%g Fx=%g\n", i, x[i], Fx[i]);
}
void J_x_dx(int n, double *dx, double *Jdx, double *x, void *par)
{
  int i;
  for(i=1;i<=n;i++)
    Jdx[i] = i*pow(x[i],i-1)*dx[i];
  Jdx[1] = 3*x[1]*x[1]*dx[1] - 8*x[1]*dx[1];
  for(i=1;i<=n;i++) printf("  i=%d dx=%g Jdx=%g x=%g\n", i, dx[i], Jdx[i], x[i]);
}

int linSol(int n, double *b, double *dx,
           void (*J_dx)(int, double *, double *, double *, void *),
           int (*precon)(int n, double *b, double *dx, double *x, void *par),
           double *x, void *par, int itmax, double tol)
{
  int i;
  double Jdx[9999];
  //preconI(n, b, dx, x, par);
  for(i=1;i<=n;i++) dx[i] = 1.0;
  J_dx(n, dx, Jdx, x, par);
  for(i=1;i<=n;i++) dx[i] = b[i]/Jdx[i];
  for(i=1;i<=n;i++) printf("linSol -> dx[%d]=%g\n", i, dx[i]);
  return i;
}
int preconI(int n, double *b, double *dx, double *x, void *par)
{
  int i;
  for(i=1;i<=n;i++) dx[i] = b[i];
  for(i=1;i<=n;i++) printf("preconI (*par=%g) -> dx[%d]=%g\n",
                           *((double *) par), i, dx[i]);
  return i;
}

#define NRANSI
#include "nrutil.c"


int main()
{
  int i, ret, check, n=3;
  double x[9];
  void *par = (void *) x;
  
  x[0]=-1;
  for(i=1;i<=n;i++) x[i]=1+0.01*i;
  x[1]=2.67;
  for(i=1;i<=n;i++) printf("x[%d]=%g\n", i, x[i]);

  ret = WT_newton(x, n, &check, F_x, J_x_dx, par, 11, 1e-11,
                  linSol, preconI, 22, 0.1);
  printf("\nret=%d\n", ret);
  for(i=1;i<=n;i++) printf("x[%d]=%.15g\n", i, x[i]);
  return 0;
}
*/
