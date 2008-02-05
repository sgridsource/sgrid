/* TOV.c */
/* Wolfgang Tichy 2007 */


#include "sgrid.h"
#include "BNSdata.h"


/* actual Gamma and K */
double Gamma, K;

/* global scaling consts, some are initialized in TOV_init */
double A,B,C,D,E,F; /* Parameters in OV Eqn */
double R,S;         /* Parameters for rho(P) */
/* scaling as in gr1.c: GR HW 1
const double rhob=1e15, Ms=1.99e33, rb=1e6, G=6.67e-8, c=3.00e10, Kb=5.38e9;
const double Gammab=1.666666666666666666666666666666667; */
/* no scaling: just G=c=1 */
const double rhob=1.0, Ms=1.0, rb=1.0, G=1.0, c=1.0, Kb=1.0;
const double Gammab=1.0;



/* TOV eqns */
void TOV_ODEs(double rf, double *y, double *dy);




/* initialize the global vars A,B,C,D,E,F, R,S, Gamma, K
   compute *rf_surf,*m,*Phi_c,*Psi_c for a given Pc, K=kappa, Gamma=Gam */
int TOV_init(double Pc, double kappa, double Gam, 
             double *rf_surf, double *m, double *Phi_c, double *Psi_c)
{
  /* Variablen fuer odeint.c */
  int kmax=23;             /* max # of points outputed by odeint */
  int kount;               /* # of points outputed by odeint */
  double *rfp,**yp;        /* points outputed by odeint  */
  double drfsav;           /* approx. rf-distance between points */
  /* meine Variablen */
  double rf1,rf2;   /* intial and final rf-point */
  double *y;        /* The functions y1, y2, ...  */
  double *dy;       /* The functions' derivs dy1, dy2, ... */
  int   nvar=4;     /* The number of functions  */
  double eps, h1, hmin;   /* error, first step, minimum step  */
  int   nok,nbad;         /* # of ok steps, # of bad steps    */
  int i, stat;
  double rfe, ret;
  double mc, Phic, Psic, zeroP;
  double M, Psi_surf, r, Psi_surf_new, Phi_surf_new;
  double Pb; /* Pbar */

  /* */
  printf("TOV_init: \n");

  /* set Gamma, K */
  Gamma = Gam;
  K = kappa;

  /* initialize scaling as in gr1.c */
  Pb=Kb*pow(rhob,Gammab);          
  A=4*PI*rhob*rb*rb*rb/Ms;
  B=G*Ms*rhob/(Pb*rb);
  C=Pb/(rhob*c*c);
  D=4*PI*Pb*rb*rb*rb/(Ms*c*c);
  E=G*Ms/(rb*c*c);
  F=E;
  R=pow(Pb/Kb,1.0/Gamma)/rhob;
  S=Pb/((Gamma-1)*c*c*rhob);
  printf(" gr1 constants set to: G=%g c=%g\n", G, c);
  printf(" A=%g B=%g C=%g D=%g E=%g=F=%g\n K=%g Gamma=%g  R=%g S=%g\n",
         A,B,C,D,E,F, K,Gamma, R,S);

  /* allocate mem. */
  y =vector(1,nvar);   /* The functions y1, y2, ... */
  dy=vector(1,nvar);   /* The functions' derivs dy1, dy2, ... */

  rfp=vector(1,kmax);         /* output of odeint */
  yp=matrix(1,nvar,1,kmax);

  /* initial values */
  mc = 0.0;
  Psic = Phic = 1.0;
  rf1= 0.0;   

  /* set initial values in y vec. */
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;
//  printf("rf=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
//         rf1, y[1], y[2], y[3], y[4]);

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
//  printf("rf2=%g zeroP=%g\n", rf2, zeroP);
//  printf("rf=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
//         rf2, y[1], y[2], y[3], y[4]);

  /* reset initial values in y vec. */
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;

  /* pars for odeintegrate */
  h1=1e-10;
  hmin=1e-10;
  eps=1e-12;
  drfsav=0.0;

  /* make one step to get away from rf=0 */
  TOV_ODEs(rf1, y, dy);
  y[1] += dy[1]*hmin;
  y[2] += dy[2]*hmin;
  y[3] += dy[3]*hmin; 
  y[4] += dy[4]*hmin; 
  rf1 += hmin;

//  printf("rf=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
//         rf1, y[1], y[2], y[3], y[4]);


  /* Here we use odeintegrate not ONLY to integrate, but also to find
     the rfe where P=0. The way we do this is odd, because we rely
     soly on the fact that odeintegrate will fail at P=0. If it does not
     fail we need a root finder to determine where P=0. */  
  rfe=rf2;
  for(;;)
  {
    ret=odeintegrate(y,nvar,rf1,rfe,eps,h1,hmin,&nok,&nbad,TOV_ODEs,rkqs,
                     kmax,&kount,rfp,yp,drfsav,&stat);  
    // printf(" ret=%g stat=%d\n", ret, stat);
    if(ret<rfe) rfe=ret;
    else break;
  }
  if( fabs((rfe-rf2)/rf2)> 0.3 || fabs(y[2])>1e-6*fabs(Pc))
    errorexit("we need a real root finder to find the rf where P=0");

  /* we want Psi = 1 + M/(2*rf) = sqrt( (2*r)/(r-M +sqrt(r*r-2*M*r)) )
     outside the star */
  /* rescale Psi and rfe */
  M = y[1];         /* total mass */
  Psi_surf = y[4];   /* current Psi at surface */
  r = rfe*Psi_surf*Psi_surf;  /* Schw. r at surface, current rf at surface */
  Psi_surf_new = sqrt( (2*r)/(r-E*M +sqrt(r*r-2*E*M*r)) ); /* desired Psi at surf. */
  Psic = Psic * Psi_surf_new/Psi_surf; /* rescale Psi at center */ 
  rfe = r/(Psi_surf_new*Psi_surf_new); /* rescale surface radius rfe */
  rf2 = rfe;
  /* add const to Phi. We need e^{2Phi} = (1.0 - 2.0*M/r) outside */
  Phi_surf_new = 0.5*log(1.0 - 2.0*E*M/r);
  Phic = Phic + (Phi_surf_new-y[3]);
  /* make one step to get away from rf=0 */
  rf1 = 0.0;
  rf2 = rfe;
  drfsav = (rf2-rf1)/kmax;
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;
  TOV_ODEs(rf1, y, dy);
  y[1] += dy[1]*hmin;
  y[2] += dy[2]*hmin;
  y[3] += dy[3]*hmin; 
  y[4] += dy[4]*hmin; 
  rf1 += hmin;

  /* check rfe */
  for(;;)
  {
    ret=odeintegrate(y,nvar,rf1,rfe,eps,h1,hmin,&nok,&nbad,TOV_ODEs,rkqs,
                     kmax,&kount,rfp,yp,drfsav,&stat);  
    printf(" ret=%g stat=%d\n", ret, stat);
    if(ret<rfe) rfe=ret;
    else break;
  }
  //printf("rf=%g:  y[1]=%g  y[2]=%g  y[3]=%g  y[4]=%g\n",
  //       rfe, y[1], y[2], y[3], y[4]);
  if( fabs((rfe-rf2)/rf2)> 0.3 || fabs(y[2])>1e-6*fabs(Pc))
    errorexit("we need a real root finder to find the rf where P=0");


  /* output results */
  for(i=1; i<=kount; i++)
    printf(" rf=%g:  m=%g  P=%g  Phi=%g  Psi=%g\n",
           rfp[i], yp[1][i], yp[2][i], yp[3][i], yp[4][i]);

  /* compute *rf_surf,*m,*Phi_c,*Psi_c */
  *rf_surf=rfe;
  *m     = y[1];
  *Phi_c = Phic;
  *Psi_c = Psic;

  free_vector(dy, 1,nvar);
  free_vector(y,  1,nvar);
  free_vector(rfp, 1,nvar);
  free_matrix(yp, 1,nvar, 1,kmax);
  return 0;
}


/* find *m, *P, *Phi, *Psi at rf, 
   for a given rf_surf, kappa, Gam, Pc, Phic, Psic */
int TOV_m_P_Phi_Psi_OF_rf(double rf, double rf_surf,
                          double kappa, double Gam,
                          double Pc, double Phic, double Psic,
                          double *m, double *P, double *Phi, double *Psi)
{
  /* Variablen fuer odeint.c */
  int kmax=2  ;            /* max # of points outputed by odeint */
  int kount;               /* # of points outputed by odeint */
  double *rfp,**yp;        /* points outputed by odeint  */
  double drfsav;           /* approx. rf-distance between points */
  /* meine Variablen */
  double rf1,rf2;   /* intial and final rf-point */
  double *y;        /* The functions y1, y2, ...  */
  double *dy;       /* The functions' derivs dy1, dy2, ... */
  int nvar=4;       /* The number of functions  */
  double eps, h1, hmin;   /* error, first step, minimum step  */
  int nok,nbad;           /* # of ok steps, # of bad steps    */
  int stat;
  double rfe, ret;
  double mc;

  y =vector(1,nvar);          /* The functions y1, y2, ... */
  dy=vector(1,nvar);          /* The functions' derivs dy1, dy2, ... */
  rfp=vector(1,kmax);         /* output of odeint */
  yp=matrix(1,nvar,1,kmax);

  /* set Gamma, K */
  Gamma = Gam;
  K = kappa;

  /* set mc=0 */
  mc = 0.0;

  /* set rf1, rf2 for int. */
  rf1=0.0;
  if(rf>rf_surf) rf2 = rf_surf;
  else           rf2 = rf;

  /* set initial values in y vec. */
  y[1]=mc;
  y[2]=Pc;
  y[3]=Phic;
  y[4]=Psic;

  /* pars for odeintegrate */
  h1=1e-10;
  hmin=1e-10;
  eps=1e-12;
  drfsav=0;

  /* make one step to get away from rf=0 */
  TOV_ODEs(rf1, y, dy);
  y[1] += dy[1]*hmin;
  y[2] += dy[2]*hmin;
  y[3] += dy[3]*hmin; 
  y[4] += dy[4]*hmin; 
  rf1 += hmin;

  /* integrate up to rf2 */
  ret=odeintegrate(y,nvar,rf1,rf2,eps,h1,hmin,&nok,&nbad,TOV_ODEs,rkqs,
                   kmax,&kount,rfp,yp,drfsav,&stat);  
  *m   = y[1];
  *P   = y[2];
  *Phi = y[3];
  *Psi = y[4];
  
  if(rf>rf_surf)
  {
    double r;

    *P   = 0.0;
    *Psi = 1.0 + E*(*m)/(2*rf);
    r = rf*(*Psi)*(*Psi);  /* Schw. r at surface, current rf at surface */
    *Phi =0.5*log(1.0 - 2.0*E*(*m)/r);
  }
  
  free_vector(dy, 1,nvar);
  free_vector(y,  1,nvar);
  free_vector(rfp, 1,nvar);
  free_matrix(yp, 1,nvar, 1,kmax);
  return 0;
}







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
//  double A,B,C,D,E,F; /* Parameters in OV Eqn */
//  double R,S;         /* Parameters for rho(P) */
//  /* scaling as in gr1.c: GR HW 1
//  const double rhob=1e15, Ms=1.99e33, rb=1e6, G=6.67e-8, c=3.00e10, Kb=5.38e9;
//  const double Gammab=1.666666666666666666666666666666667;  */
//  /* no scaling: just G=c=1 */
//  const double rhob=1.0, Ms=1.0, rb=1.0, G=1.0, c=1.0, Kb=1.0;
//  const double Gammab=1.0;
//  double Pb; /* Pbar */
//
//  /* scale as in gr1.c */
//  Pb=Kb*pow(rhob,Gammab);          
//  A=4*PI*rhob*rb*rb*rb/Ms;
//  B=G*Ms*rhob/(Pb*rb);
//  C=Pb/(rhob*c*c);
//  D=4*PI*Pb*rb*rb*rb/(Ms*c*c);
//  E=G*Ms/(rb*c*c);
//  F=E;
//  R=pow(Pb/K,1.0/Gamma)/rhob;
//  S=Pb/((Gamma-1)*c*c*rhob);

  /* retrieve vars */
  m   = y[1];
  P   = y[2];
  Phi = y[3];
  Psi = y[4];
  rho = R*pow(P,1.0/Gamma)+S*P;
//printf("Pb=%g K=%g\n", Pb, K);
//printf("R=%g S=%g Gamma=%g\n", R, S, Gamma);
//printf("m=%g P=%g Phi=%g Psi=%g rho=%g\n", m, P, Phi, Psi, rho);

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
exit(22);
*/
}
