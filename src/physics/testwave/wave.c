/* Wellengleichung : Welle auf einer Saite */

#define STRINGLEN 300
#define PI  3.1415926535897932

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "slow_Cheb_trafos.c"

/* meine Funktionen */
int vorwaerts(void);


/*Variablen*/
static double *y_old, *y, *y_new, *y2, *c, *c1, *c2; 
static double L,x1,x2,v,dx,dt;
static int N,i;


int main()
{
 int j; 

 /* Programmanfang*/
// printf("#  adv.c: \n");
 
 N=48;    /* Zahl der Gitterpunkte*/
 
 y_old = (double*) calloc(N+1,sizeof(double));
 y     = (double*) calloc(N+1,sizeof(double));
 y_new = (double*) calloc(N+1,sizeof(double));
 y2    = (double*) calloc(N+1,sizeof(double));
 c     = (double*) calloc(N+1,sizeof(double));
 c1    = (double*) calloc(N+1,sizeof(double));
 c2    = (double*) calloc(N+1,sizeof(double));
 
 x1=-1.0;
 x2=2.0;
 L=fabs(x2-x1);   
 /* kl. Abstand zwischen zwei Gitterpunkten */
 dx=0.5*L*(1-cos(PI/N));

 v=1.0;        /* Geschw. der Welle */
 dt=0.25*dx/v;  /* Zeitschrittweite */

// printf("# Ich setze jetzt die Anfangsbedingunen!\n");
 for(i=0;i<=N;i++)
 {
   double x=cos(i*PI/N); // X = 0.5*(a-b)*x +0.5*(a+b)
   double X= 0.5*(x2-x1)*x +0.5*(x1+x2);

   y[i]=exp( -10*X*X );
   y_old[i]=1.0*y[i];
 }
 
 for(i=0;i<=2*N*N;i++)  
 { 
   if(i%(N*N/16)==0)
   {
     printf("\"Time = %e\"\n", (double) i*dt);
     for(j=0;j<=N;j++)
     {
       double x=cos(j*PI/N); // X = 0.5*(a-b)*x +0.5*(a+b)
       double X= 0.5*(x2-x1)*x +0.5*(x1+x2);

       printf("%f  %f\n", X, y[j]);
     }
     printf("\n"); 
     /* printf("y[n-1]=%f   y[N]=%f \n",y[N-1],y[N]); */
   }
   
   vorwaerts();
 }

 printf("%s","Ende!\n");
 return 0; 
}



/* zu loesende Gleichung */
void varnew(double *var_new, double *var, double *var_old, 
              double *dvar, double *ddvar, 
              double dt)
{
  int i;
  for(i=0;i<=N;i++)
  {       		
    /* y(x,t+dt)= 2*y(x,t)-y(x,t-dt)+ v^2 yd*dt^2 */
    //var_new[i] = 2*var[i] - var_old[i] + (v*v * ddvar[i]) * dt*dt;
    var_new[i] = 2*var[i] - var_old[i] + ( v*v * ddvar[i]  + 
                 -1.0 * var[i]*var[i]*var[i] / (1 + var[i]*var[i]) ) *
                 dt*dt;
  }      
}


/* Einen Zeitschritt vorwaerts gehen */
int vorwaerts(void)
{
  int i;
  double *dummy;

  /* d^2 y/dt^2 = v^2 yd = v^2 d^2 y/dx^2 */
  // yd=(y[i+1]-2*y[i]+y[i-1])/(dx*dx);
  
  /* brechne y_new */
  cheb_coeffs_fromExtrema(c, y, N);
  cheb_deriv(x1,x2, c , c1, N);
  cheb_deriv(x1,x2, c1, c2, N);
  cheb_eval_onExtrema(c2, y2, N);
  varnew(y_new, y, y_old, dummy, y2, dt);

 
  /* festes Ende bei x=-1 : y[N]=0 */
  y_new[N]=0;  

  /* loses Ende bei x=1 : */
  /* d^2 y/dt^2 = v^2 /dx (- dy/dx) */
  /* d^2 y/dt^2 = v^2 /dx (-y2) =v^2 /dx (- dy/dx) */
//  y2[0]=(y[1]-y[0])/(dx*dx);
//  y_new[0]=2*y[0]-y_old[0] + v*v*y2[0]*dt*dt;

  /* Sommerfeld bei x=-1 */
  /* dy/dt = -v dy/dx */
  /* y2=dy/dx */
  cheb_eval_onExtrema(c1, y2, N);
  y_new[0] = y_old[0] - v* y2[0] * (2.0*dt);

  /* jetzt doch festes Ende bei x=-1 (dann braucht man keine Filter): */
  y_new[0] = 0.0;

  /* filter y_new */ 
  cheb_coeffs_fromExtrema(c, y_new, N);
//  cheb_filter(c, 2*N/3, N);
  cheb_eval_onExtrema(c, y_new, N);


  dummy=y_old;
  y_old=y;
  y=y_new;
  y_new=dummy; 

  return 0;
}

