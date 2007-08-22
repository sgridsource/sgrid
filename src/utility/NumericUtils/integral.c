/* Integrale mit Trapezregel */
/* Wolfgang Tichy 3.2003 */
/* integral benutzt nur die Trapezregel
   rombintegral ist Rombergintegration, d.h. Polynomextrapolation
   wird benutzt um die Konvergenz zu beschleunigen.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ************************************************************************ */

/* Funktionen in integral.c */
/* integrators */
double integral(double (*func)(double), double a, double b, double s, int max);
double rombintegral(double (*func)(double), double a, double b,
                    double eps, int max);

/* helper functions */
static double trapezit(double (*func)(double), double a, double b,
                       int n, double s_old);
double *double_vector(long nl, long nh);
void free_double_vector(double *v, long nl, long nh);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

/* ************************************************************************ */


/* Integral mit trapezit */
/* Es gibt integral und rombintegral. Beide werden gleich aufgerufen. */
/* mache maximal max Iterationen, aber breche ab wenn
             |itnew-itold|/|itold|<s                     */
double integral(double (*func)(double), double a, double b, double s, int max)
{
 double trapezit(double (*func)(double), double a, double b, int n, double s_old);
 double itold,itnew,d;
 int j;
 
 itnew=trapezit(func,a,b,1,0.0);  /* 1. Iterationsstartwert berechnen    */
 itold=itnew;                     /* 1. Iterationsstartwert in itold tun */
 for(j=2;j<=max;j++)
 {
   itnew=trapezit(func,a,b,j,itold);
   d=fabs(itnew-itold);
   if( d < s*fabs(itold) ) break;
   itold=itnew;
 }
 return itnew;
}



/* Trapezregel um Integrale zu berechnen */
/* trapezit brechnet nur den n. Anteil der Iteration und addiert diesen
   zum vorherigen. Der Vorherige muss jeweils schon berechnet sein.
   Deshalb muss trapezit so aufgerufen werden:
   for(j=1;j<=M+1;j++) integral=trapezit(&f,a,b,j);*/
static double trapezit(double (*func)(double), double a, double b, 
                       int n, double s_old)
{
 double s,x,h,tnm,sum;
 int it,j;
 
 if(n==1) { s= 0.5*(b-a)*( ((*func)(a)) + ((*func)(b)) );  return s; }

 it=1;
 for(j=1;j<n-1;j++)  it<<=1;  /* it=2^(n-2) */
 tnm=it;
 h=(b-a)/tnm;
 x=a+0.5*h;
 sum=0.0;
 for(j=1;j<=it;j++,x+=h)  sum+=((*func)(x));
 s=0.5*(s_old+(b-a)*sum/tnm);
 return s;
}





/* ***************** Stuff form NumRec ********************** */

#define NR_END 1
/* allocate a double_vector with subscript range v[nl..nh] */
double *double_vector(long nl, long nh)
{
 double *v;

 v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
 if (!v) printf("allocation failure in double_vector()\n");
 return v-nl+NR_END;
}

/* free a double_vector allocated with double_vector() */
void free_double_vector(double *v, long nl, long nh)
{
 free( (v+nl-NR_END) );
}
#undef NR_END



/* Romberg integration from NumRec, modified to take eps and max as additional
   arguments, also modified so that it call trapezit instead of the numrec 
   thing. trapezit works in 3D as well since it doesn't use a static var to 
   save s in trapezit.
*/
double rombintegral(double (*func)(double), double a, double b, 
                   double eps, int max)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapezit(double (*func)(double), double a, double b, int n, double s_old);
	
	int K=5;
	int j;
	double ss,dss;
	double *s;
	double *h;
	
	s=double_vector(0,max+2);
	h=double_vector(0,max+3);

	h[1]=1.0;
	for (j=1;j<=max;j++) {
		s[j]=trapezit(func,a,b,j,s[j-1]);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= eps*fabs(ss)) 
			{
			  free_double_vector(s,0,max+2);
			  free_double_vector(h,0,max+3);
			  return ss;
			}
		}
		h[j+1]=0.25*h[j];
	}
	printf("Accuracy goal not achieved in rombintegral!\n");
	free_double_vector(s,0,max+2);
	free_double_vector(h,0,max+3);
	return ss;
}


/* polint from NumRec modified for double prec. */
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=double_vector(1,n);
	d=double_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) printf("Error in routine polint\n");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_double_vector(d,1,n);
	free_double_vector(c,1,n);
}


/* ********* End of NumRec stuff ******************************* */

