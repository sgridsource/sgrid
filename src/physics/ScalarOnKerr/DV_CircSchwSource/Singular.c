//Singular.c

//Computes the regularizing singular field.

#include<stdio.h>
#include<math.h>
#include<assert.h>
#include "Constants.h"
#include "Schwarzschild.h"
#include "THZ.h"
#include "Singular.h"

    extern double ian_Q;

    extern double ian_x,ian_y,ian_z;
    extern double Dx[4], Dy[4], Dz[4];
    extern double DDx[4][4], DDy[4][4], DDz[4][4];

    double ian_rho;
    double ian_Psi;
    double DPsi[4];
    double DDPsi[4][4];
    
    void eval_singular(double ts, double r, double theta, double phi)
    {
        static double tsLast=0, rLast=0, thetaLast=0, phiLast=0;
        if (fabs(tsLast-ts) + fabs(rLast-r) + fabs(thetaLast-theta)+ fabs(phiLast-phi) < 1.e-14){
            return;
        }
        else{
            tsLast=ts; rLast=r; thetaLast=theta; phiLast=phi;
        }
        
        double Q=ian_Q;
        double x=ian_x;
        double y=ian_y;
        double z=ian_z;
        
        double rho,Psi;
        
        eval_g(ts,r,theta,phi);
        eval_THZ(ts,r,theta,phi);

        rho=sqrt(x*x+y*y+z*z); 
        ian_rho=rho;
        
        Psi=Q/rho;
        ian_Psi=Psi;
        
        int i;
        for(i=0;i<4;i++){
            DPsi[i]=-(Q/(rho*rho*rho))*(x*Dx[i]+y*Dy[i]+z*Dz[i]);
        }
        
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                DDPsi[i][j]=(3.0*Q/(rho*rho*rho*rho*rho))*(x*Dx[i]+y*Dy[i]+z*Dz[i])
                    *(x*Dx[j]+y*Dy[j]+z*Dz[j])-(Q/(rho*rho*rho))*(x*DDx[i][j]+Dx[i]*Dx[j]+
                    y*DDy[i][j]+Dy[i]*Dy[j]+z*DDz[i][j]+Dz[i]*Dz[j]);
            }
        }

        double Drho[4], DDrho[4][4];
   
        for(i=0;i<4;i++){
            Drho[i]=(x*Dx[i] + y*Dy[i] + z*Dz[i])/rho;
                int j;
                for(j=0;j<=i;j++){
                    DDrho[i][j] = DDrho[j][i] = (x*DDx[i][j]+Dx[i]*Dx[j]+ y*DDy[i][j]+Dy[i]*Dy[j]+z*DDz[i][j]+Dz[i]*Dz[j])/rho
                        - Drho[i]*Drho[j]/rho;
                }
        }
   
   } // end of void eval_singular()



