/* Driver for routine odeint */

#include "sgrid.h"

#define N 11


//int nrhs;   /* counts function evaluations */

/* global constants used throughout the program */
double nu, m, mu, deltam, m1, m2;
double c1, c2, c3, c4, c5, c6, c7, c8, c9, 
       c10, c11, c12, c13, c14, c15, c16,  
       c17, c18, c19, c20;
/* flags (0 or 1) that decide which order of v=om^{1/3} we include */
int f1, f2, f3, f4, f5, f6, f7;  


/* compute global constants that depend only on masses or pi,gammaE,theta_cap */
void PN_CircOrbit_compute_constants(double m1_in, double m2_in)
{
  double theta_cap, pi, gammaE;

  pi = 3.1415926535897932384626433;
  gammaE = 0.577215664901532860606512090082402431042;
  theta_cap = 1039.0/4620.0;
  m1  = m1_in;
  m2  = m2_in;
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

  /* flags (0 or 1) that decide which order of v=om^{1/3} we include */
  f1 = 1; /* order v */
  f2 = 1; /* order v^2 */
  f3 = 1; /* order v^3 */
  f4 = 1; /* order v^4 */
  f5 = 1; /* order v^5 */
  f6 = 1; /* order v^6 */
  f7 = 1; /* order v^7 */
}

/* function to be passed to odeintegrate */
void PN_CircOrbit_derivs(double x,double y[],double dydx[])
{
  int i;
  double om, v1, v2, v4, v5, v6, v7, oov1;
  double *Ln_cap, *S1, *S2, *LNS1, *LNS2, *LNS, Ln_cap_dot_S1, Ln_cap_dot_S2, S1_dot_S2, Ln_cap_mod;
  // nrhs++;

  LNS1 = dvector(1,3);
  LNS2 = dvector(1,3);
  LNS = dvector(1,3);

  om  = m*y[1];            // mw^1
  v1  = pow(om,(1.0/3.0)); // mw^1/3
  v2  = v1*v1;             // mw^2/3
  v4  = v2*v2;             // mw^4/3 
  v5  = v4*v1;             // mw^5/3
  v6  = om*om;             // mw^2
  v7  = v6*v1;	         // mw^7/3
  oov1 = 1.0/v1;           // mw^-1/3

  /* pointers S1, S2, Ln_cap to correct place inside y */
  S1 = y+1;
  S2 = y+4;
  Ln_cap = y+7;

  /* renormalize Ln_cap and thus (also inside y) */
  Ln_cap_mod = sqrt(Ln_cap[1]*Ln_cap[1]+Ln_cap[2]*Ln_cap[2]+Ln_cap[3]*Ln_cap[3]);
  
  for(i=1;i<=3;i++)
      Ln_cap[i] = Ln_cap[i]/Ln_cap_mod;

  Ln_cap_dot_S1 = Ln_cap[1]*S1[1] + Ln_cap[2]*S1[2] + Ln_cap[3]*S1[3];
  Ln_cap_dot_S2 = Ln_cap[1]*S2[1] + Ln_cap[2]*S2[2] + Ln_cap[3]*S2[3];
  S1_dot_S2     = S1[1]*S2[1] + S1[2]*S2[2] + S1[3]*S2[3];

  for(i=1;i<=3;i++)
  {
    LNS1[i] = c15*v6*( nu*c16*oov1*Ln_cap[i] + c17*(S2[i] - 3.0*Ln_cap_dot_S2*Ln_cap[i]) );
    LNS2[i] = c15*v6*( nu*c18*oov1*Ln_cap[i] + c17*(S1[i] - 3.0*Ln_cap_dot_S1*Ln_cap[i]) );
    LNS[i]  = c16*S1[i] + c18*S2[i] - c19*v1*( S1[i]*Ln_cap_dot_S2 + S2[i]*Ln_cap_dot_S1 );      
  } 
// Omega
  dydx[1] = c1*y[1]*y[1]*v5*(1.0 
                            - c2*v2 
                            - (c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5)*om 
                            + c6*v4
                            - c7*S1_dot_S2*v4 
                            + c8*Ln_cap_dot_S1*Ln_cap_dot_S2*v4 
                            - c9*v5     
                            + (c10 + c11 + c12 - c13*log(16.0*v2))*v6 
                            + c14*v7 );
// Spin1	 
  dydx[2] = LNS1[2]*S1[3] - S1[2]*LNS1[3]; 
  dydx[3] = -LNS1[1]*S1[3] + S1[1]*LNS1[3];
  dydx[4] = LNS1[1]*S1[2] - S1[1]*LNS1[2];
// Spin2
  dydx[5] = LNS2[2]*S2[3] - S2[2]*LNS2[3];
  dydx[6] = -LNS2[1]*S2[3] + S2[1]*LNS2[3]; 
  dydx[7] = LNS2[1]*S2[2] - S2[1]*LNS2[2]; 
// Ln_cap
  dydx[8] = c20*v6*( LNS[2]*Ln_cap[3] - Ln_cap[2]*LNS[3]);
  dydx[9] = c20*v6*(-LNS[1]*Ln_cap[3] + Ln_cap[1]*LNS[3]); 
  dydx[10] = c20*v6*( LNS[1]*Ln_cap[2] - Ln_cap[1]*LNS[2]);
 
  if((Ln_cap[1]==0.0)&&(Ln_cap[3]==0.0))
      dydx[11] = y[1];
  else
//        dydx(11) = y(1) - Ln_cap(2)*(dydx(8)*Ln_cap(3) - Ln_cap(1)*dydx(10))/(Ln_cap(1)*Ln_cap(1) + Ln_cap(3)*Ln_cap(3))
    dydx[11] = y[1] - Ln_cap[2]*(Ln_cap[3]*dydx[8] - Ln_cap[1]*dydx[10])/(Ln_cap[1]*Ln_cap[1] + Ln_cap[3]*Ln_cap[3]);
  free_dvector(LNS1,1,3);
  free_dvector(LNS2,1,3);
  free_dvector(LNS,1,3);  		 
}

/* get radius of orbit for a given y[], assuming PN_CircOrbit_compute_constants
   has been called before */
double PN_CircOrbit_compute_r(double y[])
{
  double om, v2, v4, oov2, LnS1, LnS2, SS, r;
  om = m*y[1];
  v2 = pow(om,(2.0/3.0));	 
  v4 = v2*v2;
  oov2 = 1.0/v2;
  LnS1 = (y[8]*y[2] + y[9]*y[3] + y[10]*y[4]);
  LnS2 = (y[8]*y[5] + y[9]*y[6] + y[10]*y[7]);
  SS = (y[2]*y[5] + y[3]*y[6] + y[4]*y[7]);
  r  = m*oov2*( 1.0 
               - (3.0 - nu)*v2/3.0 
               - om*( LnS1*(2.0*m1*m1/m/m+3.0*nu)/m1/m1 
                     +LnS2*(2.0*m1*m1/m/m+3.0*nu)/m2/m2 )/3.0
               + ( nu*(19.0/4.0 + nu/9.0) 
                   - 0.5*nu*(SS - 3.0*LnS1*LnS2)/m1/m1/m2/m2 )*v4);
  return r;
}

/* routine that calls odeintegrate, and integrates from t1 to t2 */
void PN_CircOrbit_xodeint(double m1_in, double m2_in, double t1, double t2, double ystart_in[])
{
  int status; /* added by WT */
  int i,nbad,nok;
  int kmax, kount;
  double dxsav, *xp, **yp;
  double eps = 1.0e-9,
         h1 = 0.01,
         hmin = 0.0, 
         x1 = t1,
         x2 = t2,
         *ystart;

  ystart = dvector(1,N);
  for(i=1;i<=N;i++) ystart[i]=ystart_in[i];
//         printf("%s %10.4e %16.6e \n","PHI before",x1, ystart[11]);
  kmax=10;
  dxsav=(x2-x1)/(kmax-1);
  xp = dvector(1,kmax*2);
  yp = dmatrix(1,11,1,kmax*2);

  /* set global consts */
  PN_CircOrbit_compute_constants(m1_in, m2_in);

  // nrhs=0;
  /* changed by WT */
  odeintegrate(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,PN_CircOrbit_derivs,rkqs,
                       kmax,&kount,xp,yp,dxsav, &status);
  if(status!=0)
    printf("odeintegrate: status=%d\n", status);
//	printf("\n%s %13s %3d\n","successful steps:"," ",nok);
//	printf("%s %20s %3d\n","bad steps:"," ",nbad);
//	printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
//	printf("\n%s %3d\n","stored intermediate values:    ",kount);
//	for (i=1;i<=kount;i++)
//		printf("%10.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n",
//			xp[i],yp[1][i],yp[2][i],yp[5][i],yp[10][i],yp[11][i]);
//        printf("%10.4e %16.6e \n",x2, ystart[1]);
  for(i=1;i<=N;i++) ystart_in[i]=ystart[i];
//         printf("%s %10.4e %16.6e \n","PHI after",x2, ystart_in[11]); 
  ystart_in[0] = PN_CircOrbit_compute_r(ystart);      
  free_dmatrix(yp,1,10,1,kmax*2);
  free_dvector(xp,1,kmax*2);
  free_dvector(ystart,1,N);
  return;
}
