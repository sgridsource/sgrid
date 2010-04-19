/* changed by WT */
#include "sgrid.h"
/*
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
*/
#define NRANSI

void compute_hcross_hplus(double y[], double *hcross, double *hplus, double D, double theta, double phi, double m1, double m2)
{
   int i,j;   
   double
     *es1, 
     *es2, 
     *n, 
     *lambda, 
     **h,  
     *etheta, 
     *ephi, 
     **Tcross, 
     **Tplus;   

//    printf("%s %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e\n","y in compute_GW:", 
//            y[0], y[1], y[2], y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10], y[11]);
   es1 = dvector(1,3); 
   es2 = dvector(1,3);
   n   = dvector(1,3);
   lambda = dvector(1,3); 
   h      = dmatrix(1,3,1,3);  
   etheta = dvector(1,3); 
   ephi   = dvector(1,3); 
   Tcross = dmatrix(1,3,1,3);  
   Tplus  = dmatrix(1,3,1,3);

//printf("%s \n","after allocation");
   es1[1] = -y[9]/sqrt(y[8]*y[8] + y[9]*y[9]);
   es1[2] =  y[8]/sqrt(y[8]*y[8] + y[9]*y[9]);
   es1[3] =  0.0;  

//printf("%s  %10.6e %10.6e %10.6e \n","es1",es1[1],es1[2],es1[3]);
  
   es2[1] = -y[8]*y[10]/sqrt(y[8]*y[8] + y[9]*y[9]);
   es2[2] = -y[9]*y[10]/sqrt(y[8]*y[8] + y[9]*y[9]);
   es2[3] = (1.0 - y[10]*y[10])/sqrt(y[8]*y[8] + y[9]*y[9]);

//printf("%s  %10.6e %10.6e %10.6e \n","es2",es2[1],es2[2],es2[3]);
   etheta[1] = cos(theta)*cos(phi); 
   etheta[2] = cos(theta)*sin(phi);
   etheta[3] = -sin(theta);

//printf("%s  %10.6e %10.6e %10.6e \n","etheta",etheta[1],etheta[2],etheta[3]);

   ephi[1] = -sin(phi); 
   ephi[2] =  cos(phi);
   ephi[3] = 0.0;

//printf("%s  %10.6e %10.6e %10.6e \n","ephi",ephi[1],ephi[2],ephi[3]);

   for(i=1;i<=3;i++)
   {  
     for(j=1;j<=3;j++)
     {
//       printf("%s %3d \n","i:", i);
//       printf("%s%3d \n","j:",j);
       Tcross[i][j] = etheta[i]*ephi[j] + ephi[i]*etheta[j];
       Tplus[i][j] = etheta[i]*etheta[j] - ephi[i]*ephi[j];
     }  
   }
    
//printf("%s \n","after Tcross, Tplus");
   for(i=1;i<=3;i++)
   {
//       printf("%s%3d\n","i:", i); 
       n[i] = es1[i]*cos(y[11]) + es2[i]*sin(y[11]);
       lambda[i] = -es1[i]*sin(y[11]) + es2[i]*cos(y[11]);
   }
//printf("%s\n","after n lambda");
   for(i=1;i<=3;i++)
   {
     for(j=1;j<=3;j++)
     {
//       printf("%s %3d\n","i:", i);
//       printf("%s%3d\n","j:",j); 
        h[i][j] = 2.0*(lambda[i]*lambda[j]-n[i]*n[j]);
     }
   }
//printf("%s\n", "after h[i][j]"); 
//   m   = m1 + m2;
//   mu  = m1*m2/m;
   for(i=1;i<=3;i++)
   {
     for(j=1;j<=3;j++)
     {
//         printf("%s %3d \n","i:", i);
//         printf("%s%3d \n","j:",j);          
         h[i][j] = 2.0*(m1*m2)/D/y[0]*h[i][j];
     }
   }
//printf("%s\n","after hij"); 

   *hcross=0.5*(h[1][1]*Tcross[1][1] + h[1][2]*Tcross[1][2] + h[1][3]*Tcross[1][3] 
                 +h[2][1]*Tcross[2][1] + h[2][2]*Tcross[2][2] + h[2][3]*Tcross[2][3] 
                 +h[3][1]*Tcross[3][1] + h[3][2]*Tcross[3][2] + h[3][3]*Tcross[3][3]);

   *hplus=0.5*(h[1][1]*Tplus[1][1] + h[1][2]*Tplus[1][2] + h[1][3]*Tplus[1][3] 
                 +h[2][1]*Tplus[2][1] + h[2][2]*Tplus[2][2] + h[2][3]*Tplus[2][3] 
                 +h[3][1]*Tplus[3][1] + h[3][2]*Tplus[3][2] + h[3][3]*Tplus[3][3]);
   
   free_dvector(es1,1,3); 
   free_dvector(es2,1,3);
   free_dvector(n,1,3);
   free_dvector(lambda,1,3); 
   free_dmatrix(h,1,3,1,3);  
   free_dvector(etheta,1,3); 
   free_dvector(ephi,1,3); 
   free_dmatrix(Tcross,1,3,1,3);  
   free_dmatrix(Tplus,1,3,1,3);
   return;
}
#undef NRANSI
