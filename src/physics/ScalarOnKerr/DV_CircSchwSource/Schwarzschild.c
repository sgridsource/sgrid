//Schwarzschild.c

//Contains elements specific to the Schwarzschild metric.

#include<stdio.h>
#include<math.h>
#include "Constants.h"
#include "Schwarzschild.h"

    extern double ian_M;
    
    double g[4][4];
    double detg;
    double Ddetg[4];
    double g_inv[4][4];
    double div_g_inv[4];

    double r2rstar(double r)
    {
        double M=ian_M;
        
        return r+2.*M*log(r/(2.*M)-1.);
    }

    double rstar2r(double r_star)
    {
        double M=ian_M;

        double x,xnew,xold,f,f_prime,c;
        int times=0;

        if(r_star>10.*M){
            c=(r_star-2.0*M)/(2.0*M);
            x=r_star;
            f=-1.;

            while(fabs(f)>1.e-15 && times<40){
                xold=x;
                f=x+log(x)-c;
                f_prime=(x+1.)/x;
                x=x-f/f_prime;
                xnew=x;
                if(fabs((xnew-xold)/xold)<1.e-14) break;
                times++;
            }
        } else{
            c=exp((r_star-2.*M)/(2.*M));
            x=2.5*M;
            f=-1.;

            while(fabs(f)>1.e-15 && times<40){
                xold=x;
                f=x*exp(x)-c;
                f_prime=(x+1.)*exp(x);
                x=x-f/f_prime;
                xnew=x;
                if(fabs((xnew-xold)/xold)<1.e-14) break;
                times++;
                }
        }

        x=2.*M*(x+1.);

        if(r_star>-10.0*M){
            double testRstar;
            testRstar=x+2.*M*log(x/(2.*M)-1.)-r_star;
            //assert(fabs(testRstar)<1.e-13);
        }

        return x;
    }

    void eval_g(double ts, double r, double theta, double phi)
    {
        double M=ian_M;
        
        detg=r*r*sin(theta);

        Ddetg[0]=0.0;
        Ddetg[1]=2.0*r*sin(theta);
        Ddetg[2]=r*r*cos(theta);
        Ddetg[3]=0.0;

        div_g_inv[0]=0.0;
        div_g_inv[1]=2.0*M/(r*r);
        div_g_inv[2]=0.0;
        div_g_inv[3]=0.0;

        int i;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                if(i!=j){
                    g[i][j]=0.0;
                    g_inv[i][j]=0.0;
                }
            }
        }

        g[0][0]=-(1.0-2.0*M/r);
        g[1][1]=1.0/(1.0-2.0*M/r);
        g[2][2]=r*r;
        g[3][3]=r*r*sin(theta)*sin(theta);

        g_inv[0][0]=1.0/g[0][0];
        g_inv[1][1]=1.0/g[1][1];
        g_inv[2][2]=1.0/g[2][2];
        g_inv[3][3]=1.0/g[3][3];
    }
