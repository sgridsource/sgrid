//Source.c

//Computes the regularized source for a scalar point charge in a circular
//orbit in Schwarzschild. 

//The source may be evaluated with Schwarzschild coordinates using
//SourceInSchwarzschild(double ts, double r, double theta, double phi), or 
//with Kerr-Schild coordinates with 
//SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS).

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include "Constants.h"
#include "Schwarzschild.h"
#include "THZ.h"
#include "Singular.h"
#include "Window_DONUT.h"
#include "Source.h"

extern double ian_M;
extern double ian_R;
extern double ian_Q;
extern double ian_Omega;

extern double g[4][4];
extern double detg;
extern double Ddetg[4];
extern double g_inv[4][4];
extern double div_g_inv[4];

extern double ian_t,ian_x,ian_y,ian_z;

extern double ian_Psi;
extern double DPsi[4];
extern double DDPsi[4][4];

extern double ian_W;
extern double DW[4];
extern double DDW[4][4];

double WPsi;
double Del2Psi;

void eval_WPsi(double ts, double r, double theta, double phi)
{
    double W=ian_W;
    double Psi=ian_Psi;

    eval_singular(ts,r,theta,phi);
    eval_window(ts,r,theta,phi);
    
    WPsi=W*Psi;
}

void eval_source(double ts, double r, double theta, double phi)
{
    double sum=0.0;
    double R=ian_R;
    double M=ian_M;
    double Pi=3.14159265358979323846264338;
    
    if(fabs(r-R)>3.5*M  && ian_W<1.0){
        Del2Psi=0.0;    
    }
    else{                                    
    
    eval_THZ(ts,r,theta,phi);
    
    double t=ian_t;
    double x=ian_x;
    double y=ian_y;
    double z=ian_z;
    
    if(sqrt(x*x+y*y+z*z)<0.005*M)                                       
    {                           
        r=R+0.006*M;
        
        eval_THZ(ts,r,theta,phi);
        t=ian_t;                           
        x=ian_x;
        y=ian_y;
        z=ian_z;
    }
    
    if(theta==0.0)                                       
    {                           
        theta+=1.0e-12;
        
        eval_THZ(ts,r,theta,phi);
        t=ian_t;                             
        x=ian_x;
        y=ian_y;
        z=ian_z;
    }
    
    if(fabs(theta-Pi)<1.0e-12)                                       
    {                           
        theta-=1.0e-12;
        
        eval_THZ(ts,r,theta,phi);
        t=ian_t;                             
        x=ian_x;
        y=ian_y;
        z=ian_z;
    }
    
    
    eval_g(ts,r,theta,phi);
    eval_singular(ts,r,theta,phi);
    eval_window(ts,r,theta,phi);
    
    double Psi=ian_Psi;
    double W=ian_W;
        
    int i;
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            sum+=Ddetg[i]*g_inv[i][j]*DPsi[j];
        }
    }
    sum=sum/detg;
    for(i=0;i<4;i++){
        sum+=div_g_inv[i]*DPsi[i];
    }
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            sum+=g_inv[i][j]*DDPsi[i][j];
        }
    }

                                        //window modifications.
    sum=sum*W;
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            sum+=(Psi/detg)*Ddetg[i]*g_inv[i][j]*DW[j];
        }
    }
    for(i=0;i<4;i++){
        sum+=Psi*div_g_inv[i]*DW[i];
    }
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            sum+=2.0*g_inv[i][j]*DPsi[i]*DW[j];
        }
    }
    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
            sum+=Psi*g_inv[i][j]*DDW[i][j];
        }
    }
    Del2Psi=sum;
    }
}
    
double SourceInSchwarzschild(double ts, double r, double theta, double phi)
{    
    eval_source(ts,r,theta,phi);
    return Del2Psi;
}

double SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS)
{
    double M=ian_M;  
    
    double Pi = 3.14159265358979323846264338;
    double r = sqrt(xKS*xKS+yKS*yKS+zKS*zKS);
    
    double theta; 
    if(zKS != 0.0){ 
        if(zKS>0.0) theta = atan(sqrt(xKS*xKS+yKS*yKS)/zKS);
        else theta = Pi-atan(sqrt(xKS*xKS+yKS*yKS)/fabs(zKS));
    }
    else theta = Pi/2.0;
    
    double phi;
    if(xKS == 0.0){
        if(yKS>= 0.0) phi = Pi/2.0;
        else phi = -Pi/2.0;
    } 
    else if(xKS > 0.0){ 
        if(yKS>=0.0) phi = atan(yKS/xKS);
        else phi = 2.0*Pi-atan(fabs(yKS)/xKS);
    }
    else{ 
        if(yKS>=0.0) phi = Pi-atan(yKS/fabs(xKS));
        else phi = Pi+atan(fabs(yKS/xKS));
    }
    
    double ts = tKS-2.0*M*log(r/(2.0*M)-1.0); 
        
    assert(fabs(xKS-r*sin(theta)*cos(phi))<1.e-12); 
    assert(fabs(yKS-r*sin(theta)*sin(phi))<1.e-12);
    assert(fabs(zKS-r*cos(theta))<1.e-12);
   
    eval_source(ts,r,theta,phi);
    return Del2Psi;
}

double SingularInKerrSchild(double tKS, double xKS, double yKS, double zKS)
{
    double M=ian_M;  
    
    double Pi = 3.14159265358979323846264338;
    double r = sqrt(xKS*xKS+yKS*yKS+zKS*zKS);
    
    double theta; 
    if(zKS != 0.0){ 
        if(zKS>0.0) theta = atan(sqrt(xKS*xKS+yKS*yKS)/zKS);
        else theta = Pi-atan(sqrt(xKS*xKS+yKS*yKS)/fabs(zKS));
    }
    else theta = Pi/2.0;
    
    double phi;
    if(xKS == 0.0){
        if(yKS>= 0.0) phi = Pi/2.0;
        else phi = -Pi/2.0;
    } 
    else if(xKS > 0.0){ 
        if(yKS>=0.0) phi = atan(yKS/xKS);
        else phi = 2.0*Pi-atan(fabs(yKS)/xKS);
    }
    else{ 
        if(yKS>=0.0) phi = Pi-atan(yKS/fabs(xKS));
        else phi = Pi+atan(fabs(yKS/xKS));
    }
    
    double ts = tKS-2.0*M*log(r/(2.0*M)-1.0); 
        
    assert(fabs(xKS-r*sin(theta)*cos(phi))<1.e-12); 
    assert(fabs(yKS-r*sin(theta)*sin(phi))<1.e-12);
    assert(fabs(zKS-r*cos(theta))<1.e-12);
   
    eval_THZ(ts,r,theta,phi);
    eval_g(ts,r,theta,phi);
    eval_singular(ts,r,theta,phi);
    return ian_Psi;
}

double current_x(double t_KS)
{
    double R=ian_R;
    double M=ian_M;
    double Omega=ian_Omega;
    double Pi=3.14159265358979323846264338;
    
    double ts = t_KS-2.0*M*log(R/(2.0*M)-1.0);
    
    return R*sin(Pi/2.0)*cos(Omega*ts);  
}

double current_y(double t_KS)
{
    double R=ian_R;
    double M=ian_M;
    double Omega=ian_Omega;
    double Pi=3.14159265358979323846264338;
 
    double ts = t_KS-2.0*M*log(R/(2.0*M)-1.0); 
    
    return R*sin(Pi/2.0)*sin(Omega*ts);  
}

double current_z(double t_KS)
{   
    return 0.0;  
}

