//Constants.c

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "Constants.h"

double ian_M;
double ian_Q=1.0;
double ian_R; 
double ian_Omega;
double ian_Omegad;

double ian_n;
double ian_width;

double Kdelta[4][4];

void set_parameters(double m, double q, double r, int N, double width)
{
    if(r<6.0*m){                                                     //effectice source only works for R>6M.
        perror("Orbit too close to the black hole. Pick R>6M. \n");
        abort();
    }

    ian_M = m;
    double M=ian_M;
    
    ian_R = r*M;   
    double R=ian_R;                               //radius of the orbit
    
    ian_Omega = sqrt(M/(R*R*R));                  //frequency for circular orbit
    ian_Omegad = ian_Omega*sqrt(1.0-3.0*M/R);

    ian_n = N;
    ian_width = width;
}

void eval_Kdelta()                       //Kronecker delta
{
    int i;
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            if(i==j) Kdelta[i][j]=1;
            else Kdelta[i][j]=0;
        }
    }
}
