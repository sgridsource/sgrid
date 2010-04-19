/* Driver for routine odeint */

/* changed by WT */
#include "sgrid.h"
/*
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
*/
#define NRANSI

#define N 11

double dxsav,*xp,**yp;  /* defining declarations */
int kmax,kount;

int nrhs;   /* counts function evaluations */

//global constants used througth out all program
double nu, m, hi1, hi2, mu, deltam, m1, m2;
double c1, c2, c3, c4, c5, c6, c7, c8, c9, 
       c10, c11, c12, c13, c14, c15, c16,  
       c17, c18, c19, c20;

void compute_constants(void)
{
    double theta_cap, pi, gammaE;

    pi = 3.1415926535897932384626433;
    gammaE = 0.577215664901532860606512090082402431042;
    theta_cap = 1039.0/4620.0;
    m   = m1 + m2;
    mu  = m1*m2/m;
    nu  = mu/m;
    deltam = m1 - m2;

    c1  = 96.0*nu/5.0;
//    printf("%s %10.6e \n","c1:", c1);  
    c2  = (743.0 + 924.0*nu)/336.0;
//    printf("%s %10.6e \n","c2:", c2); 
    c3  = (113.0*m1*m1/m/m + 75*nu)/12.0/m1/m1;
    c4  = (113.0*m2*m2/m/m + 75*nu)/12.0/m2/m2;
    c5  = 4.0*pi;
    c6  = 34103.0/18144.0 + 13661.0/2016.0*nu + 59.0/18.0*nu*nu;
    c7  = 247.0*nu/m1/m1/m2/m2/48.0;
    c8  = 721.0*nu/m1/m1/m2/m2/48.0;
    c9  = (4159.0 + 15876.0*nu)*pi/672.0;
    c10 = 16447322263.0/139708800.0 -1712.0/105.0*gammaE + 16.0/3.0*pi*pi;
    c11 = (-273811877.0/1088640.0+451.0/48.0*pi*pi-88.0/3.0*theta_cap)*nu;
    c12 = 541.0/896.0*nu*nu - 5605.0/2592.0*nu*nu*nu;
    c13 = 856.0/105.0;
    c14 = (-4415.0/4032.0 + 358675.0/6048.0*nu + 91495.0/1512.0*nu*nu)*pi;
    c15 = 0.5/m;
    c16 = 4.0 + 3.0*m2/m1;
    c17 = 1.0/m/m;
    c18 = 4.0 + 3.0*m1/m2;
    c19 = 3.0/nu/m/m;
    c20 = 0.5/m/m/m;           
}
void derivs(double x,double y[],double dydx[])
{
//	dydx[1] = -y[2];
//	dydx[2]=y[1]-(1.0/x)*y[2];
//	dydx[3]=y[2]-(2.0/x)*y[3];
//	dydx[4]=y[3]-(3.0/x)*y[4];
        int i;
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	double *Ln_cap, *S1, *S2, *LNS1, *LNS2, *LNS, Ln_cap_dot_S1, Ln_cap_dot_S2, S1_dot_S2, Ln_cap_mod;

        Ln_cap = dvector(1,3); 
        S1 = dvector(1,3); 
        S2 = dvector(1,3); 
        LNS1 = dvector(1,3);
        LNS2 = dvector(1,3);
        LNS = dvector(1,3);

        t0  = m*y[1];            // mw^1
        t1  = pow(t0,(1.0/3.0)); //  ! mw^1/3
        t2  = t1*t1;          // mw^2/3
	t3  = t2*t2;          // mw^4/3 
	t4  = t3*t1;          // mw^5/3
        t5  = t0*t0;          // mw^2
        t6  = t5*t1;	      // mw^7/3
        t7  = 1.0/t1;         // mw^-1/3

        for(i=1;i<=3;i++){
            S1[i] = y[i+1];
            S2[i] = y[i+4];
            Ln_cap[i] = y[i+7];
        }
        Ln_cap_mod = sqrt(Ln_cap[1]*Ln_cap[1]+Ln_cap[2]*Ln_cap[2]+Ln_cap[3]*Ln_cap[3]);
        
        for(i=1;i<=3;i++)
            Ln_cap[i] = Ln_cap[i]/Ln_cap_mod;

        Ln_cap_dot_S1 = Ln_cap[1]*S1[1] + Ln_cap[2]*S1[2] + Ln_cap[3]*S1[3];
        Ln_cap_dot_S2 = Ln_cap[1]*S2[1] + Ln_cap[2]*S2[2] + Ln_cap[3]*S2[3];
        S1_dot_S2     = S1[1]*S2[1] + S1[2]*S2[2] + S1[3]*S2[3];

        for(i=1;i<=3;i++)
        {
            LNS1[i] = c15*t5*( nu*c16*t7*Ln_cap[i] + c17*(S2[i] - 3.0*Ln_cap_dot_S2*Ln_cap[i]) );
            LNS2[i] = c15*t5*( nu*c18*t7*Ln_cap[i] + c17*(S1[i] - 3.0*Ln_cap_dot_S1*Ln_cap[i]) );
            LNS[i]  = c16*S1[i] + c18*S2[i] - c19*t1*( S1[i]*Ln_cap_dot_S2 + S2[i]*Ln_cap_dot_S1 );      
        } 
// Omega
        dydx[1] = c1*y[1]*y[1]*t4*(1.0 
                                  - c2*t2 
                                  - (c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5)*t0 
                                  + c6*t3
                                  - c7*S1_dot_S2*t3 
                                  + c8*Ln_cap_dot_S1*Ln_cap_dot_S2*t3 
                                  - c9*t4     
                                  + (c10 + c11 + c12 - c13*log(16.0*t2))*t5 
                                  + c14*t6 );
// Spin1	 
	dydx[2] =  LNS1[2]*S1[3] - S1[2]*LNS1[3]; 
	dydx[3] = -LNS1[1]*S1[3] + S1[1]*LNS1[3];
	dydx[4] =  LNS1[1]*S1[2] - S1[1]*LNS1[2];
// Spin2
        dydx[5] =  LNS2[2]*S2[3] - S2[2]*LNS2[3];
	dydx[6] = -LNS2[1]*S2[3] + S2[1]*LNS2[3]; 
	dydx[7] =  LNS2[1]*S2[2] - S2[1]*LNS2[2]; 
// Ln_cap
        dydx[8] = c20*t5*( LNS[2]*Ln_cap[3] - Ln_cap[2]*LNS[3]);
	dydx[9] = c20*t5*(-LNS[1]*Ln_cap[3] + Ln_cap[1]*LNS[3]); 
	dydx[10] = c20*t5*( LNS[1]*Ln_cap[2] - Ln_cap[1]*LNS[2]);

        if((Ln_cap[1]==0.0)&&(Ln_cap[2]==0.0))
            dydx[11] = 0.0;
        else
            dydx[11] = y[1] - Ln_cap[3]*(Ln_cap[1]*dydx[9] - Ln_cap[2]*dydx[8])/(Ln_cap[1]*Ln_cap[1] + Ln_cap[2]*Ln_cap[2]);
        free_dvector(Ln_cap,1,3);
        free_dvector(S1,1,3);
        free_dvector(S2,1,3);
        free_dvector(LNS1,1,3);
        free_dvector(LNS2,1,3);
        free_dvector(LNS,1,3);  		 
}

double compute_r(double y[])
{
    double t1, t2, t3, t4, t5, t6, t7, r;
        t1 = m*y[1];
	t2 = pow(t1,(2.0/3.0));	 
	t3 = 1.0/t2;
	t4 = (y[8]*y[2] + y[9]*y[3] + y[10]*y[4]);
        t5 = (y[8]*y[5] + y[9]*y[6] + y[10]*y[7]);
	t6 = (y[2]*y[5] + y[3]*y[6] + y[4]*y[7]);
        t7 = t2*t2;
	r  = m*t3*( 1.0 
                  - (3.0 - nu)*t2/3.0 
                  - t1*(t4*(2.0*m1*m1/m/m+3.0*nu)/m1/m1 
                  + t5*(2.0*m1*m1/m/m+3.0*nu)/m2/m2)/3.0
                  + ( nu*(19.0/4.0 + nu/9.0) 
                      - 0.5*nu*(t6 - 3.0*t4*t5)/m1/m1/m2/m2 )*t7);
    return r;	    
}

void xodeint(double m1_in, double m2_in, double hi1_in, double hi2_in, double t1, double t2, double ystart_in[])
{
        int status; /* added by WT */
	int i,nbad,nok;
	double eps = 1.0e-6,
               h1 = 0.1,
               hmin = 0.0, 
               x1 = t1,
               x2 = t2,
               *ystart;

	ystart = dvector(1,N);
        for(i=1;i<=N;i++)
            ystart[i]=ystart_in[i];
	xp = dvector(1,200);
	yp = dmatrix(1,11,1,200);
        m1 = m1_in;
        m2 = m2_in;
        hi1 = hi1_in;
        hi2 = hi2_in;
        compute_constants();
//	ystart[1] = 0.01;
//	ystart[2] = 0.1;
//	ystart[3] = 0.2;
//	ystart[4] = 0.3;
//      ystart[5] =-0.1;
//        ystart[6] = 0.2;
//        ystart[7] = 0.3;
//        ystart[8] = 0.0;
//        ystart[9] = 0.0;
//        ystart[10] = 1.0;
//        ystart[11]= 0.0;
// printf("%s %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e\n","ystart:", 
//            ystart[1], ystart[2], ystart[3],ystart[4],ystart[5],ystart[6],ystart[7],ystart[8],ystart[9],ystart[10], ystart[11]);
	nrhs=0;
	kmax=100;
	dxsav=(x2-x1)/(kmax-1);
/*	odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);
*/	/* changed by WT */
	odeintegrate(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs,
	             kmax,&kount,xp,yp,dxsav, &status);
	             
//	printf("\n%s %13s %3d\n","successful steps:"," ",nok);
//	printf("%s %20s %3d\n","bad steps:"," ",nbad);
//	printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
//	printf("\n%s %3d\n","stored intermediate values:    ",kount);
//	printf("\n%8s %18s %15s\n","x","integral","bessj(3,x)");
//	for (i=1;i<=kount;i++)
//		printf("%10.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n",
//			xp[i],yp[1][i],yp[2][i],yp[5][i],yp[10][i],yp[11][i]);
//        printf("%10.4e %16.6e \n",x2, ystart[1]);
        for(i=1;i<=N;i++)
            ystart_in[i]=ystart[i];
        ystart_in[0] = compute_r(ystart);      
	free_dmatrix(yp,1,10,1,200);
	free_dvector(xp,1,200);
	free_dvector(ystart,1,N);
        return;
}
#undef NRANSI
