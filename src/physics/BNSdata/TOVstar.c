/* TOV.c */
/* Wolfgang Tichy 2007 */


#include "sgrid.h"
#include "BNSdata.h"


/* actual Gamma amd K */
double Gamma, K;

/* TOV eqns */
void TOV_ODEs(double rf, double *y, double *dy);



int ode1()
{
  /* Variablen fuer odeint.c */
  int kmax=100;            /* max # of points outputed by odeint */
  int kount;               /* # of points outputed by odeint */
  double *rfp,**yp;          /* points outputed by odeint  */
  double drfsav;             /* approx. rf-distance between points */
  /* meine Variablen */
  double rf1,rf2;   /* intial and final rf-point */
  double *y;      /* The functions y1, y2, ...  */
  double *dy;     /* The functions' derivs dy1, dy2, ... */
  int   nvar=4;  /* The number of functions  */
  double eps, h1, hmin;   /* error, first step, minimum step  */
  int   nok,nbad;        /* # of ok steps, # of bad steps    */
  int i, stat;
  double rfe, ret;
  double mc, Pc, Phic, Psic, zeroP;

  y =vector(1,nvar);   /* The functions y1, y2, ... */
  dy=vector(1,nvar);   /* The functions' derivs dy1, dy2, ... */

  rfp=vector(1,kmax);         /* output of odeint */
  yp=matrix(1,nvar,1,kmax);

  /* set Gamm, K */
  Gamma = 1.6666666666666666;
//  K = 5.38e9;
  K = 1;
  printf("parameter Gamma=%g :\n",Gamma);

  /* initial values */
  mc = 0.0;
  Pc = K*pow(4, Gamma);
  Phic = 1.0;
  Psic =1.0;

  printf("#  ODE as an Initial Value Problem: \n\n");
  printf("integrate using rkqs:\n");


  rf1=0.0;   

  /* set initial values in y vec. */
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;

  printf("rf1=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
         rf1, y[1], y[2], y[3], y[4]);

  /* find rf2 where P is approx zero */
  rf2=rf1;
  hmin=1e-4;
  while(y[2]>=0.0)  /* y[2]=P */
  {
    zeroP = y[2]; /* save last val of P*/
    TOV_ODEs(rf2, y, dy);
    y[1] += dy[1]*hmin;
    y[2] += dy[2]*hmin;
    y[3] += dy[3]*hmin; 
    y[4] += dy[4]*hmin; 
    rf2 += hmin;
  }
  /*increase rf2 by 10% */
  rf2 *= 1.1;
  printf("rf2=%g zeroP=%g\n", rf2, zeroP);

  /* reset initial values in y vec. */
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;

  /* pars for odeintegrate */
  h1=1e-10;
  hmin=1e-10;
  eps=1e-12;
  drfsav=0.1;

  /* make one step to get away from rf=0 */
  TOV_ODEs(rf1, y, dy);
/*
  y[1] += dy[1]*(hmin*hmin);
  y[2] += dy[2]*(hmin*hmin);
  y[3] += dy[3]*(hmin*hmin); 
  y[4] += dy[4]*(hmin*hmin); 
  rf1 += (hmin*hmin);
*/
  y[1] += dy[1]*hmin;
  y[2] += dy[2]*hmin;
  y[3] += dy[3]*hmin; 
  y[4] += dy[4]*hmin; 
  rf1 += hmin;

  printf("rf=rf1=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
         rf1, y[1], y[2], y[3], y[4]);

  /* Here we use odeintegrate not ONLY to integrate, but also to find
     the rfe where P=0. The way we do this is odd, because we rely
     soly on the fact that odeintegrate will fail at P=0. If it does not
     fail we need a root finder to determine where P=0. */  
  rfe=rf2;
  for(;;)
  {
    ret=odeintegrate(y,nvar,rf1,rfe,eps,h1,hmin,&nok,&nbad,TOV_ODEs,rkqs,
                     kmax,&kount,rfp,yp,drfsav,&stat);  
    printf("ret=%g stat=%d\n", ret, stat);
    if(ret<rfe) rfe=ret;
    else break;
  }
  printf("rf=rfe=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n\n",
         rf1, y[1], y[2], y[3], y[4]);
  if( fabs((rfe-rf2)/rf2)> 0.3 || fabs(y[2])>1e-6*fabs(Pc))
    errorexit("we need a real root finder to find the rf where P=0");
 
  for(i=1; i<=kount; i++)
    printf("rf=%g:  m=%g  P=%g  Phi=%g  Psi=%g\n",
           rfp[i], yp[1][i], yp[2][i], yp[3][i], yp[4][i]);

  /* rescale Psi and rfe */
// ...
  /* add const to Phi */
// ...
   
  free_vector(dy, 1,nvar);
  free_vector(y,  1,nvar);
  free_vector(rfp, 1,nvar);
  free_matrix(yp, 1,nvar, 1,kmax);
  return 0;
}


/* TOV eqns */
/* r = area radial coord, rf = isotrop. like coord for conf. flatness */
void TOV_ODEs(double rf, double *y, double *dy)
{
  double r, dm_dr, dP_dr, dPhi_dr, dPsi_dr;
  double dr_drf, dm_drf, dP_drf, dPhi_drf, dPsi_drf;
  double m, P, Phi, Psi, rho;
  double A,B,C,D,E,F; /* Parameters in OV Eqn */
  double R,S;         /* Parameters for rho(P) */
  /* scaling as in gr1.c: GR HW 1
  const double rhob=1e15, Ms=1.99e33, rb=1e6, G=6.67e-8, c=3.00e10, Kb=5.38e9;
  const double Gammab=1.666666666666666666666666666666667;  */
  /* no scaling: just G=c=1 */
  const double rhob=1.0, Ms=1.0, rb=1.0, G=1.0, c=1.0, Kb=1.0;
  const double Gammab=1.0;
  double Pb; /* Pbar */

  /* scale as in gr1.c */
  Pb=Kb*pow(rhob,Gammab);          
  A=4*PI*rhob*rb*rb*rb/Ms;
  B=G*Ms*rhob/(Pb*rb);
  C=Pb/(rhob*c*c);
  D=4*PI*Pb*rb*rb*rb/(Ms*c*c);
  E=G*Ms/(rb*c*c);
  F=E;
  R=pow(Pb/K,1.0/Gamma)/rhob;
  S=Pb/((Gamma-1)*c*c*rhob);

  /* retrieve vars */
  m   = y[1];
  P   = y[2];
  Phi = y[3];
  Psi = y[4];
  rho = R*pow(P,1.0/Gamma)+S*P;
//printf("rho=%g P=%g\n", rho, P);

  /* r in terms of rf */
  r = rf*Psi*Psi;
  
  /* derivs with respect to r */
  dm_dr = A*rho*r*r;
  if(r>1e-4) // (1e-40)*E*rb)
  {
    dP_dr   =-B*rho*m/(r*r)*(1.0+C*P/rho)*(1.0+D*P*r*r*r/m)/(1.0-E*2.0*m/r);
    dPhi_dr = F*m/(r*r)*(1.0+D*P*r*r*r/m)/(1.0-E*2.0*m/r);
    dPsi_dr = 0.5*(1/r - 1/sqrt(r*r - 2.0*E*m*r))*Psi;
  }
  else
  {
    dP_dr   =-A*B/3.0*rho*rho*r*(1.0+C*P/rho)*(1.0+D*3.0/A*P/rho) / 
              (1.0-A/3.0*E*2.0*rho*r*r);
    dPhi_dr = A*F/3.0*rho*r*(1.0+D*3.0/A*P/rho)/(1.0-A/3.0*E*2.0*rho*r*r);
    dPsi_dr = 0.5*(-A*E*rho*r/3.0)*Psi;
  }

  /* derivs with respect to rf */
  dr_drf   = Psi*Psi/(1 - 2.0*rf*Psi*dPsi_dr);
  dm_drf   = dm_dr * dr_drf;
  dP_drf   = dP_dr * dr_drf;
  dPhi_drf = dPhi_dr * dr_drf;
  dPsi_drf = dPsi_dr * dr_drf;
  dy[1] = dm_drf;
  dy[2] = dP_drf;
  dy[3] = dPhi_drf;
  dy[4] = dPsi_drf;
/*
printf("r=%g m=%g ", r, m);
printf("P=%g ", P);
printf("Phi=%g ", Phi);
printf("Psi=%g ", Psi);
printf("\n");
printf("dr_drf=%g dm_dr=%g ", dr_drf, dm_dr);
printf("dP_dr=%g ", dP_dr);
printf("dPhi_dr=%g ", dPhi_dr);
printf("dPsi_dr=%g ", dPsi_dr);
printf("\n");
*/
}
