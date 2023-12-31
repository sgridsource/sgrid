//Window_DONUT.c

//Implements the window function of the effective source.

#include<stdio.h>
#include<math.h>
#include "Constants.h"
#include "Schwarzschild.h"
#include "THZ.h"
#include "Singular.h"
#include "Window_DONUT.h"


    extern double ian_M;
    extern double ian_R;
    extern double ian_n;
    extern double ian_width; 

    extern double Kdelta[4][4];

    double ian_W;
    double DW[4];
    double DDW[4][4];

double powr(double r, int n){
    double temp=1.0;
    int i;
    for(i=0;i<n;i++){
        temp=temp*r;
    }
    return temp;
}

void eval_window(double ts, double r, double theta, double phi)
{
  static int useWindow=-1;
  double M=ian_M;
  double R=ian_R;
  double width=ian_width;
  double n=ian_n;

  double W;
  //double width=2.0*M; // width used in paper
  double rl = (r-R)/width;
  //int n = 8;                  
  int i,j;

  if(useWindow<0) useWindow=Getv("DV_CircSchwSource_useWindow","yes");
  
  if(useWindow)
  {
    eval_Kdelta();
    eval_singular(ts,r,theta,phi);
    eval_THZ(ts,r,theta,phi);
    
    W=exp(-powr(rl,n));
    ian_W=W;
    
    for(i=0;i<4;i++){
        DW[i]=-W*(n*powr(rl,n-1)*(Kdelta[1][i])/width) ;
    }
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            DDW[i][j]=-DW[j]*(n*powr(rl,n-1)*(Kdelta[1][i])/width)
            -W*(n*(n-1)*powr(rl,n-2)*Kdelta[1][i]*Kdelta[1][j]/(width*width));
            }
    }
  }
  else /* W=1 */
  {
    ian_W=W=1.0;
    for(i=0;i<4;i++)
    {
      DW[i]=0.0;
      for(j=0;j<4;j++) DDW[i][j]=0.0;
    }
  }
}
