/* ODE as an Initial Value Problem  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include "nr.h"       // Funktionen aus numrec ohne nrutil */
#include "nrutil.h"  

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));

/* meine Funktionen */
void derivs(double x, double *y, double *dy);

/* meine Variablen */
double q=0.0;          /* ein parameter in meiner Dgl. */
const double pi=3.1415926535897932384626;

int main()
{
 /* Variablen fuer odeint.c */
 int kmax=100;            /* max # of points outputed by odeint */
 int kount;               /* # of points outputed by odeint */
 double *xp,**yp;          /* points outputed by odeint  */
 double dxsav;             /* approx. x-distance between points */
 /* meine Variablen */
 double x1,x2;   /* intial and final x-point */
 double *y;      /* The functions y1, y2, ...  */
 double *dy;     /* The functions' derivs dy1, dy2, ... */
 int   nvar=3;  /* The number of functions  */
 double eps, h1, hmin;   /* error, first step, minimum step  */
 int   nok,nbad;        /* # of ok steps, # of bad steps    */

 y =vector(1,nvar);   /* The functions y1, y2, ... */
 dy=vector(1,nvar);   /* The functions' derivs dy1, dy2, ... */

 xp=vector(1,kmax);         /* output of odeint */
 yp=matrix(1,nvar,1,kmax);


 /* Programmanfang*/
 printf("#  ODE as an Initial Value Problem: \n\n");
 printf("integrate using rkqs:\n");
 printf("parameter q=%f :\n",q);
 x1=0.0;   
 x2=pi;
 y[1]=1.0;   /* initial values */
 y[2]=0.0;
 y[3]=1.0;
 printf("x=x1=%f:  y[1]=%f  y[2]=%f  y[3]=%f \n",x1,y[1],y[2],y[3]);
 
 h1=0.0000001;      hmin=0.0000001;
 eps=0.000001;
 odeintegrate(y,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs,
             kmax,&kount,xp,yp,dxsav);  
 printf("x=x2=%f:  y[1]=%f  y[2]=%f  y[3]=%f \n",x2,y[1],y[2],y[3]);

 free_vector(dy,1,2);
 free_vector(y,1,2);
 return 0;
}

void derivs(double x, double *y, double *dy)
{
 dy[1]= y[2];
 dy[2]=-y[3]*y[1] +2*q*cos(2*x)*y[1];
 dy[3]= 0.0;
}
