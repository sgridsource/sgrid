//THZ.c

//Implements the transformation from THZ to Schwarzchild coordinates.

#include<stdio.h>
#include<math.h>
#include "Constants.h"
#include "Schwarzschild.h"
#include "THZ.h"

    extern double ian_M;
    extern double ian_R; 
    extern double ian_Omega;
    extern double ian_Omegad;

    extern double Kdelta[4][4];

    double ian_X,ian_Y,ian_Z;                //All of these are functions of (ts,r,theta, phi).
    double DX[4],DY[4],DZ[4];                //Capital letters stand for the helping functions.
    double DDX[4][4],DDY[4][4],DDZ[4][4];    //Small letters are the THZ coordinates.

    double ian_t,ian_x,ian_y,ian_z;
    double Dx[4], Dy[4], Dz[4];
    double DDx[4][4], DDy[4][4], DDz[4][4];

    void eval_THZ(double ts, double r, double theta, double phi)
    {
        static double tsLast=0, rLast=0, thetaLast=0, phiLast=0;
        if (fabs(tsLast-ts) + fabs(rLast-r) + fabs(thetaLast-theta)+ fabs(phiLast-phi) < 1.e-14){
            return;
        }
        else{
            tsLast=ts; rLast=r; thetaLast=theta; phiLast=phi;
        }
        
        double M=ian_M;
        double R=ian_R;
        double Omega=ian_Omega;
        double Omegad=ian_Omegad;
        
        double X,Y,Z,t,x,y,z;
        
        eval_Kdelta();

        double cosTheta=cos(theta);
        double sinTheta=sin(theta);
        double cosPhi=cos(phi-Omega*ts);
        double sinPhi=sin(phi-Omega*ts);
        double cosTime=cos(Omegad*ts);
        double sinTime=sin(Omegad*ts);
                                                //Helping functions X,Y,Z, and their derivatives.
        X=(r-R)/(sqrt(1.0-2.0*M/R));                 
        Y=R*sinTheta*sinPhi*sqrt((R-2.0*M)/(R-3.0*M));
        Z=R*cosTheta;

        DX[0]=0.0;
        DX[1]=1.0/(sqrt(1.0-2.0*M/R));
        DX[2]=0.0;
        DX[3]=0.0;

        double c=R*sqrt((R-2.0*M)/(R-3.0*M));

        DY[0]=-c*Omega*sinTheta*cosPhi;
        DY[1]=0.0;
        DY[2]=c*cosTheta*sinPhi;
        DY[3]=c*sinTheta*cosPhi;

        DZ[0]=0.0;
        DZ[1]=0.0;
        DZ[2]=-R*sinTheta;
        DZ[3]=0.0;

        int i;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                DDX[i][j]=0.0;
                DDY[i][j]=0.0;
                DDZ[i][j]=0.0;
            }
        }
        
        DDY[0][0]=-c*Omega*Omega*sinTheta*sinPhi;
        DDY[0][2]=-c*Omega*cosTheta*cosPhi;
        DDY[2][0]=DDY[0][2];
        DDY[0][3]=c*Omega*sinTheta*sinPhi;
        DDY[3][0]=DDY[0][3];
        DDY[2][2]=-c*sinTheta*sinPhi;
        DDY[2][3]=c*cosTheta*cosPhi;
        DDY[3][2]=DDY[2][3];
        DDY[3][3]=-c*sinTheta*sinPhi;
        
        DDZ[2][2]=-R*cosTheta;

                                                //Evaluates THZ coordinates (t,x,y,z) and their derivatives.

        t=ts*sqrt(1.0-3.0*M/R)-(r*Omega*Y)/(sqrt(1-2.0*M/R))
            +((Omega*M*Y)/(R*R*sqrt(1-2.0*M/R)*(R-3.0*M)))
                 *(-0.5*X*X*(R-M)+M*Y*Y*(R-3.0*M)/(3.0*(R-2.0*M))+M*Z*Z)
            +(Omega*M*X*Y/(14.0*R*R*R*(R-2.0*M)*(R-3.0*M)))
            *(-X*X*(R*R-11.0*R*M+11.0*M*M)+Y*Y*(13.0*R*R-45.0*R*M+31.0*M*M)
                   +Z*Z*(13.0*R-5.0*M)*(R-2.0*M));

                                                //This is the part for z(ts,r,theta,phi).

        double A1,B1,B2,B3,C1,D1,D2,D3;

        A1=M/(2.0*R*R*R*(R-3.0*M));
        B1=-(2.0*R-3.0*M);
        B2=(R-3.0*M);
        B3=(R-2.0*M);
        C1=M/(14.0*R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        D1=M*(13.0*R-19.0*M);
        D2=14.0*R*R-36.0*R*M+9.0*M*M;
        D3=(R-2.0*M)*(14.0*R-15.0*M);

        z=r*cosTheta+A1*Z*(B1*X*X+B2*Y*Y+B3*Z*Z)+C1*X*Z*(D1*X*X+D2*Y*Y+D3*Z*Z);

        for(i=0;i<4;i++){

            Dz[i]=cosTheta*Kdelta[1][i]-r*sinTheta*Kdelta[2][i]
                +A1*DZ[i]*(B1*X*X+B2*Y*Y+B3*Z*Z)
                +2.0*A1*Z*(B1*X*DX[i]+B2*Y*DY[i]+B3*Z*DZ[i])
                +C1*DX[i]*Z*(D1*X*X+D2*Y*Y+D3*Z*Z)
                +C1*X*DZ[i]*(D1*X*X+D2*Y*Y+D3*Z*Z)+2.0*C1*X*Z*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDz[i][j]=-sinTheta*Kdelta[1][i]*Kdelta[2][j]
                        -sinTheta*Kdelta[1][j]*Kdelta[2][i]
                        -r*cosTheta*Kdelta[2][i]*Kdelta[2][j]
                        +A1*DDZ[i][j]*(B1*X*X+B2*Y*Y+B3*Z*Z)
                        +2.0*A1*DZ[i]*(B1*X*DX[j]+B2*Y*DY[j]+B3*Z*DZ[j])
                        +2.0*A1*DZ[j]*(B1*X*DX[i]+B2*Y*DY[i]+B3*Z*DZ[i])
                        +2.0*A1*Z*(B1*DX[i]*DX[j]+B1*X*DDX[i][j]+B2*DY[i]*DY[j]
                            +B2*Y*DDY[i][j]+B3*DZ[i]*DZ[j]+B3*Z*DDZ[i][j])
                        +C1*DDX[i][j]*Z*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +C1*DX[i]*DZ[j]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +2.0*C1*DX[i]*Z*(D1*X*DX[j]+D2*Y*DY[j]+D3*Z*DZ[j])
                        +C1*DX[j]*DZ[i]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +C1*X*DDZ[i][j]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +2.0*C1*X*DZ[i]*(D1*X*DX[j]+D2*Y*DY[j]+D3*Z*DZ[j])
                        +2.0*C1*DX[j]*Z*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i])
                        +2.0*C1*X*DZ[j]*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i])
                        +2.0*C1*X*Z*(D1*DX[i]*DX[j]+D1*X*DDX[i][j]+D2*DY[i]*DY[j]
                            +D2*Y*DDY[i][j]+D3*DZ[i]*DZ[j]+D3*Z*DDZ[i][j]);

            }
        }

                                //This is for x_tilde and y_tilde, functions which are used to evaluate x and y.

        double xt,yt;
        double  Dxt[4], Dyt[4];
        double DDxt[4][4], DDyt[4][4];

                                //y_tilde.

        double E1,E2,E3;

        A1=sqrt((R-2.0*M)/(R-3.0*M));
        B1=M/(2.0*R*R*R);
        C1=(R-3.0*M)/(R-2.0*M);
        D1=M/(14.0*R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        E1=2.0*M*(4.0*R-15.0*M);
        E2=14.0*R*R-69.0*M*R+89.0*M*M;
        E3=2.0*(R-2.0*M)*(7.0*R-24.0*M);

        yt=A1*r*sinTheta*sinPhi+B1*Y*(-2.0*X*X+C1*Y*Y+Z*Z)
            +D1*X*Y*(E1*X*X+E2*Y*Y+E3*Z*Z);
            
        for(i=0;i<4;i++){

            Dyt[i]=A1*sinTheta*sinPhi*Kdelta[1][i]
                  +A1*r*cosTheta*sinPhi*Kdelta[2][i]
                  +A1*r*sinTheta*cosPhi*Kdelta[3][i]
                  -A1*Omega*r*sinTheta*cosPhi*Kdelta[0][i]
                  +B1*DY[i]*(-2.0*X*X+C1*Y*Y+Z*Z)
                  +2.0*B1*Y*(-2.0*X*DX[i]+C1*Y*DY[i]+Z*DZ[i])
                  +D1*(Y*DX[i]+X*DY[i])*(E1*X*X+E2*Y*Y+E3*Z*Z)
                  +2.0*D1*X*Y*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDyt[i][j]=A1*cosTheta*sinPhi*Kdelta[1][i]*Kdelta[2][j]
                          +A1*sinTheta*cosPhi*Kdelta[1][i]*Kdelta[3][j]
                          -A1*Omega*sinTheta*cosPhi*Kdelta[1][i]*Kdelta[0][j]
                          +A1*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[1][j]
                          -A1*r*sinTheta*sinPhi*Kdelta[2][i]*Kdelta[2][j]
                          +A1*r*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[3][j]
                          -A1*Omega*r*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[0][j]
                          +A1*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[1][j]
                          +A1*r*cosTheta*cosPhi*Kdelta[3][i]*Kdelta[2][j]
                          -A1*r*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[3][j]
                          +A1*Omega*r*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[0][j]
                          -A1*Omega*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[1][j]
                          -A1*Omega*r*cosTheta*cosPhi*Kdelta[0][i]*Kdelta[2][j]
                          +A1*Omega*r*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[3][j]
                          -A1*Omega*Omega*r*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[0][j]
                          +B1*DDY[i][j]*(-2.0*X*X+C1*Y*Y+Z*Z)
                          +2.0*B1*DY[i]*(-2.0*X*DX[j]+C1*Y*DY[j]+Z*DZ[j])
                          +2.0*B1*DY[j]*(-2.0*X*DX[i]+C1*Y*DY[i]+Z*DZ[i])
                          +2.0*B1*Y*(-2.0*DX[j]*DX[i]-2.0*X*DDX[i][j]+C1*DY[i]*DY[j]+C1*Y*DDY[i][j]
                            +DZ[i]*DZ[j]+Z*DDZ[i][j])
                          +D1*DDX[i][j]*Y*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +D1*DX[i]*DY[j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*DX[i]*Y*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +D1*DX[j]*DY[i]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +D1*X*DDY[i][j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*X*DY[i]*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +2.0*D1*DX[j]*Y*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*DY[j]*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*Y*(E1*DX[i]*DX[j]+E1*X*DDX[i][j]+E2*DY[i]*DY[j]+E2*Y*DDY[i][j]
                            +E3*DZ[i]*DZ[j]+E3*Z*DDZ[i][j]);

            }
        }

                                                //x_tilde.

        double F1,G1,G2,G3,G4,G5,G6;

        A1=1.0/(sqrt(1.0-2.0*M/R));
        B1=M/(R*R*sqrt(1.0-2.0*M/R));
        C1=(R-3.0*M)/(R-2.0*M);
        D1=M/(2.0*R*R*R*(R-2.0*M)*(R-3.0*M));
        E1=-M*M;
        E2=(R-3.0*M)*(3.0*R-8.0*M);
        E3=3.0*(R-2.0*M)*(R-2.0*M);
        F1=M/(R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        G1=M*(R*R-R*M+3.0*M*M)/(8.0*(R-2.0*M));
        G2=(28.0*R*R-114.0*R*M+123.0*M*M)/28.0;
        G3=(14.0*R*R-48.0*R*M+33.0*M*M)/14.0;
        G4=M*(3.0*R*R*R-74.0*R*R*M+337.0*R*M*M-430.0*M*M*M)/(56.0*(R-2.0*M)*(R-2.0*M));
        G5=-M*M*(7.0*R-18.0*M)/(4.0*(R-2.0*M));
        G6=-M*(3.0*R+22.0*M)/56.0;

        xt=A1*(r*sinTheta*cosPhi-R)+B1*(-0.5*X*X+C1*Y*Y+Z*Z)
            +D1*X*(E1*X*X+E2*Y*Y+E3*Z*Z)
            +F1*(G1*X*X*X*X+G2*X*X*Y*Y+G3*X*X*Z*Z+G4*Y*Y*Y*Y+G5*Y*Y*Z*Z+G6*Z*Z*Z*Z);

        for(i=0;i<4;i++){

            Dxt[i]=A1*sinTheta*cosPhi*Kdelta[1][i]
                  +A1*r*cosTheta*cosPhi*Kdelta[2][i]
                  -A1*r*sinTheta*sinPhi*Kdelta[3][i]
                  +A1*Omega*r*sinTheta*sinPhi*Kdelta[0][i]
                  +B1*(-X*DX[i]+2.0*C1*Y*DY[i]+2.0*Z*DZ[i])
                  +D1*DX[i]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                  +2.0*D1*X*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                  +F1*(4.0*G1*X*X*X*DX[i]+2.0*G2*X*DX[i]*Y*Y+2.0*G2*X*X*Y*DY[i]
                    +2.0*G3*X*DX[i]*Z*Z+2.0*G3*X*X*Z*DZ[i]+4.0*G4*Y*Y*Y*DY[i]
                    +2.0*G5*Y*DY[i]*Z*Z+2.0*G5*Y*Y*Z*DZ[i]+4.0*G6*Z*Z*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDxt[i][j]=A1*cosTheta*cosPhi*Kdelta[1][i]*Kdelta[2][j]
                          -A1*sinTheta*sinPhi*Kdelta[1][i]*Kdelta[3][j]
                          +A1*Omega*sinTheta*sinPhi*Kdelta[1][i]*Kdelta[0][j]
                          +A1*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[1][j]
                          -A1*r*sinTheta*cosPhi*Kdelta[2][i]*Kdelta[2][j]
                          -A1*r*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[3][j]
                          +A1*Omega*r*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[0][j]
                          -A1*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[1][j]
                          -A1*r*cosTheta*sinPhi*Kdelta[3][i]*Kdelta[2][j]
                          -A1*r*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[3][j]
                          +A1*Omega*r*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[0][j]
                          +A1*Omega*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[1][j]
                          +A1*Omega*r*cosTheta*sinPhi*Kdelta[0][i]*Kdelta[2][j]
                          +A1*Omega*r*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[3][j]
                          -A1*Omega*Omega*r*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[0][j]
                          +B1*(-DX[i]*DX[j]-X*DDX[i][j]+2.0*C1*DY[j]*DY[i]+2.0*C1*Y*DDY[i][j]
                               +2.0*DZ[i]*DZ[j]+2.0*Z*DDZ[i][j])
                          +D1*DDX[i][j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*DX[i]*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +2.0*D1*DX[j]*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*(E1*DX[i]*DX[j]+E1*X*DDX[i][j]+E2*DY[i]*DY[j]+E2*Y*DDY[i][j]
                                     +E3*DZ[i]*DZ[j]+E3*Z*DDZ[i][j])
                          +F1*(12.0*G1*X*X*DX[j]*DX[i]+4.0*G1*X*X*X*DDX[i][j]+2.0*G2*DX[j]*DX[i]*Y*Y
                              +2.0*G2*X*DDX[i][j]*Y*Y+4.0*G2*X*DX[i]*Y*DY[j]+4.0*G2*X*DX[j]*Y*DY[i]
                              +2.0*G2*X*X*DY[j]*DY[i]+2.0*G2*X*X*Y*DDY[i][j]+2.0*G3*DX[j]*DX[i]*Z*Z
                              +2.0*G3*X*DDX[i][j]*Z*Z+4.0*G3*X*DX[i]*Z*DZ[j]+4.0*G3*X*DX[j]*Z*DZ[i]
                              +2.0*G3*X*X*DZ[j]*DZ[i]+2.0*G3*X*X*Z*DDZ[i][j]+12.0*G4*Y*Y*DY[i]*DY[j]
                              +4.0*G4*Y*Y*Y*DDY[i][j]+2.0*G5*DY[j]*DY[i]*Z*Z+2.0*G5*Y*DDY[i][j]*Z*Z
                              +4.0*G5*Y*DY[i]*Z*DZ[j]+4.0*G5*Y*DY[j]*Z*DZ[i]+2.0*G5*Y*Y*DZ[j]*DZ[i]
                              +2.0*G5*Y*Y*Z*DDZ[i][j]+12.0*G6*Z*Z*DZ[j]*DZ[i]+4.0*G6*Z*Z*Z*DDZ[i][j]);
            }
        }

                                            //x(ts,r,theta,phi) and derivatives.

        x=xt*cosTime-yt*sinTime;

        for(i=0;i<4;i++){

            Dx[i]=Dxt[i]*cosTime-xt*Omegad*sinTime*Kdelta[0][i]
                -Dyt[i]*sinTime-yt*Omegad*cosTime*Kdelta[0][i];

        }
        
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDx[i][j]=DDxt[i][j]*cosTime
                        -(Dxt[i]*Kdelta[0][j]+Dxt[j]*Kdelta[0][i])*Omegad*sinTime
                        -xt*Omegad*Omegad*cosTime*Kdelta[0][i]*Kdelta[0][j]
                        -DDyt[i][j]*sinTime
                        -(Dyt[i]*Kdelta[0][j]+Dyt[j]*Kdelta[0][i])*Omegad*cosTime
                        +yt*Omegad*Omegad*sinTime*Kdelta[0][i]*Kdelta[0][j];
            }
        }
                                            //y(ts,r,theta,phi) and derivatives.

        y=xt*sinTime+yt*cosTime;

        for(i=0;i<4;i++){
        
            Dy[i]=Dxt[i]*sinTime+xt*Omegad*cosTime*Kdelta[0][i]
                +Dyt[i]*cosTime-yt*Omegad*sinTime*Kdelta[0][i];

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDy[i][j]=DDxt[i][j]*sinTime
                        +(Dxt[i]*Kdelta[0][j]+Dxt[j]*Kdelta[0][i])*Omegad*cosTime
                        -xt*Omegad*Omegad*sinTime*Kdelta[0][i]*Kdelta[0][j]
                        +DDyt[i][j]*cosTime
                        -(Dyt[i]*Kdelta[0][j]+Dyt[j]*Kdelta[0][i])*Omegad*sinTime
                        -yt*Omegad*Omegad*cosTime*Kdelta[0][i]*Kdelta[0][j];
            }
        }   
    
        ian_X=X;
        ian_Y=Y;
        ian_Z=Z;
        ian_t=t;
        ian_x=x;
        ian_y=y;
        ian_z=z;
    }

