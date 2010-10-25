/* compute hplus and hcross up to 1.5 PN in Amplitude */

#include "sgrid.h"

void compute_P05Qc(double lambda[], double n[], double theta, double phi, double **P05Q, double deltam, double m);
void compute_P1Qc(double lambda[], double n[], double theta, double phi, double **P1Q, double deltam, double m, double Delta[], double nu);
void compute_P15Qc(double y[], double lambda[], double n[], double theta, double phi, double **P15Q, double deltam, double m1, double m2, double nu);
void compute_P2QSSc(double y[], double n[], double theta, double phi, double **P2QSSc, double m1, double m2);

void compute_hcross_hplus(double y[], double *hcross, double *hplus, double D, double theta, double phi, double m1, double m2)
{
  int i,j;
  double r, deltam, m, nu;   
  double
    *es1, 
    *es2, 
    *n, 
    *lambda, 
    **h,  
    *etheta, 
    *ephi,
    *Delta, 
    **Tcross, 
    **Tplus,
    **P05Q,
    **P1Q,
    **P15Q,
    **P2QSSc;
  int AOv1, AOv2, AOv3, AOv4; /* flags for amplitude */

  /* set flags for amplitude */
  AOv1 = Geti("PN_CircularOrbit_GWs_AmpOv1");
  AOv2 = Geti("PN_CircularOrbit_GWs_AmpOv2");
  AOv3 = Geti("PN_CircularOrbit_GWs_AmpOv3");
  AOv4 = Geti("PN_CircularOrbit_GWs_AmpOv4");

//    printf("%s %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e %10.6e %10.6e \n %10.6e\n","y in compute_GW:", 
//            y[0], y[1], y[2], y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10], y[11]);
  r = y[0];
  deltam = m1 - m2;
  m = m1 + m2;
  nu  = m1*m2/m/m;
  es1 = dvector(1,3); 
  es2 = dvector(1,3);
  n   = dvector(1,3);
  Delta = dvector(1,3);
  lambda = dvector(1,3); 
  h      = dmatrix(1,3,1,3);  
  etheta = dvector(1,3); 
  ephi   = dvector(1,3); 
  Tcross = dmatrix(1,3,1,3);  
  Tplus  = dmatrix(1,3,1,3);
  P05Q   = dmatrix(1,3,1,3);
  P1Q    = dmatrix(1,3,1,3);
  P15Q   = dmatrix(1,3,1,3);
  P2QSSc = dmatrix(1,3,1,3);

  /* compute Delta from spins */
  Delta[1] = m*(y[5]/m2 - y[2]/m1);
  Delta[2] = m*(y[6]/m2 - y[3]/m1);
  Delta[3] = m*(y[7]/m2 - y[4]/m1);

  /* let es1, es2 and L_hat for a right handed Dreibein */
  /* compute vector es1 from L_hat: es1 = ey \cross L_hat / (sin i)
     ey /cross L_hat = (L_hat_z, 0 , -L_hat_x)
     (sin i)^2 = L_hat_x^2 + L_hat_z^2  */
  if(y[8]==0.0&&y[10]==0.0)
  {
    es1[1] = 1.0;
    es1[2] = 0.0;
    es1[3] = 0.0;
  }
  else
  { 
    es1[1] =  y[10]/sqrt(y[8]*y[8] + y[10]*y[10]);
    es1[2] =  0.0;
    es1[3] = -y[8]/sqrt(y[8]*y[8] + y[10]*y[10]);
  }  
//printf("%s  %10.6e %10.6e %10.6e \n","es1",es1[1],es1[2],es1[3]);

  /* compute vector es2: es2 = L_hat \cross es1 = (ey - L_hat cos i) / (sin i)
     cos i = L_hat_y,  (sin i)^2 = L_hat_x^2 + L_hat_z^2   */
  if(y[8]==0.0&&y[10]==0.0)
  {
    es2[1] = 0.0;
    es2[2] = 0.0;
    es2[3] = -y[9]/fabs(y[9]);
  }
  else
  {
    es2[1] =       -y[8]*y[9]/sqrt(y[8]*y[8] + y[10]*y[10]);
    es2[2] = (1.0 - y[9]*y[9])/sqrt(y[8]*y[8] + y[10]*y[10]);
    es2[3] =       -y[9]*y[10]/sqrt(y[8]*y[8] + y[10]*y[10]);
  }
//printf("%s  %10.6e %10.6e %10.6e \n","es2",es2[1],es2[2],es2[3]);

  /* compute etheta */
  etheta[1] = cos(theta)*cos(phi); 
  etheta[2] = cos(theta)*sin(phi);
  etheta[3] = -sin(theta);
//printf("%s  %10.6e %10.6e %10.6e \n","etheta",etheta[1],etheta[2],etheta[3]);

  /* compute ephi */
  ephi[1] = -sin(phi); 
  ephi[2] =  cos(phi);
  ephi[3] =  0.0;
//printf("%s  %10.6e %10.6e %10.6e \n","ephi",ephi[1],ephi[2],ephi[3]);

  /* compute Tplus, Tcross tensors from etheta, ephi */
  for(i=1;i<=3;i++)
  {  
    for(j=1;j<=3;j++)
    {
//       printf("%s %3d \n","i:", i);
//       printf("%s%3d \n","j:",j);
      Tcross[i][j] = etheta[i]*ephi[j] + ephi[i]*etheta[j];
      Tplus[i][j] = etheta[i]*etheta[j] - ephi[i]*ephi[j];
    }  
  }
//printf("%s \n","after Tcross, Tplus");

  /* compute  n = es1 cos(Phi) + es2 sin(Phi) 
              lambda = -es1 sin(Phi) + es2 cos(Phi)  */
  for(i=1;i<=3;i++)
  {
//       printf("%s%3d\n","i:", i); 
      n[i] = es1[i]*cos(y[11]) + es2[i]*sin(y[11]);
      lambda[i] = -es1[i]*sin(y[11]) + es2[i]*cos(y[11]);
  }
//printf("%s\n","after n lambda");

  /* compute post 0.5, 1, 1.5 PN terms in h_ij */
  compute_P05Qc(lambda, n, theta, phi, P05Q, deltam, m);
  compute_P1Qc(lambda, n, theta, phi, P1Q, deltam, m, Delta, nu);
  compute_P15Qc(y, lambda, n, theta, phi, P15Q, deltam, m1, m2, nu);  
  compute_P2QSSc(y, n, theta,phi, P2QSSc, m1,m2);

  /* get h_ij / [2mu/D * m/r] */
  for(i=1;i<=3;i++)
  {
    for(j=1;j<=3;j++)
    {
       double Qij = 2.0*(lambda[i]*lambda[j]-n[i]*n[j]);  /* Eqn (4.9a) of Kidder, PRD 52, 821 (1995) */

       h[i][j] = Qij
                 + P05Q[i][j]*sqrt(m/r)    * AOv1
                 + P1Q[i][j]*(m/r)         * AOv2
                 + P15Q[i][j]*pow(m/r,1.5) * AOv3
                 + P2QSSc[i][j]*m*m/(r*r)  * AOv4;
    }
  }
//printf("%s\n", "after h[i][j]"); 
//   m   = m1 + m2;
//   mu  = m1*m2/m;

  /* multiply h[i][j] by [2mu/D * m/r] to obtain h_ij */
  for(i=1;i<=3;i++)
  {
    for(j=1;j<=3;j++)
    {
//         printf("%s %3d \n","i:", i);
//         printf("%s%3d \n","j:",j);          
        h[i][j] = 2.0*(m1*m2)/D/y[0]*h[i][j];
    }
  }
//printf("%s\n","after hij"); 

  /* compute h+ abd hx */
  *hcross=0.5*( h[1][1]*Tcross[1][1] + h[1][2]*Tcross[1][2] + h[1][3]*Tcross[1][3] 
               +h[2][1]*Tcross[2][1] + h[2][2]*Tcross[2][2] + h[2][3]*Tcross[2][3] 
               +h[3][1]*Tcross[3][1] + h[3][2]*Tcross[3][2] + h[3][3]*Tcross[3][3]);

  *hplus=0.5*( h[1][1]*Tplus[1][1] + h[1][2]*Tplus[1][2] + h[1][3]*Tplus[1][3] 
              +h[2][1]*Tplus[2][1] + h[2][2]*Tplus[2][2] + h[2][3]*Tplus[2][3] 
              +h[3][1]*Tplus[3][1] + h[3][2]*Tplus[3][2] + h[3][3]*Tplus[3][3]);

  /* free memory */
  free_dvector(es1,1,3); 
  free_dvector(es2,1,3);
  free_dvector(n,1,3);
  free_dvector(lambda,1,3);
  free_dvector(Delta,1,3); 
  free_dmatrix(h,1,3,1,3);  
  free_dvector(etheta,1,3); 
  free_dvector(ephi,1,3); 
  free_dmatrix(Tcross,1,3,1,3);  
  free_dmatrix(Tplus,1,3,1,3);
  free_dmatrix(P05Q,1,3,1,3);
  free_dmatrix(P1Q,1,3,1,3);
  free_dmatrix(P15Q,1,3,1,3);
  free_dmatrix(P2QSSc,1,3,1,3);
  return;
}

//**********************************************************
//* subprogram compute 0.5 post-Newtonian correction to hij
//**********************************************************
void compute_P05Qc(double lambda[], double n[], double theta, double phi, double **P05Q, double deltam, double m)
{
  int i,j;
  double NN_dot_n, NN_dot_lambda;
  double *NN;

  NN = dvector(1,3);
   
  /* get N_hat */
  NN[1] = sin(theta)*cos(phi);
  NN[2] = sin(theta)*sin(phi);
  NN[3] = cos(theta);
//printf("%10.6e %10.6e %10.6e \n", lambda[1], lambda[2], lambda[3]);
//printf("NN");
//printf("%10.6e %10.6e %10.6e \n", NN[1], NN[2], NN[3]);

  /* set scalar products */
  NN_dot_n = NN[1]*n[1]+NN[2]*n[2]+NN[3]*n[3];
  NN_dot_lambda = NN[1]*lambda[1]+NN[2]*lambda[2]+NN[3]*lambda[3];
//printf("NN_dot_n, NN_dot_lambda");
//printf("%10.6e %10.6e \n", NN_dot_n, NN_dot_lambda);

  /* Eqn (4.9b) of Kidder, PRD 52, 821 (1995) */
  for(i=1;i<=3;i++){
      for(j=1;j<=3;j++){
          P05Q[i][j] = deltam/m*(3.0*NN_dot_n*( n[i]*lambda[j]+n[j]*lambda[i])
                                 +NN_dot_lambda*( n[i]*n[j]-2.0*lambda[i]*lambda[j]));
  
      }
  }
//    printf("%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n", 
//                               P05Q[1][1], P05Q[1][2], P05Q[1][3],
//                               P05Q[2][1], P05Q[2][2], P05Q[2][3],
//                               P05Q[3][1], P05Q[3][2], P05Q[3][3]);
  free_dvector(NN,1,3);
  return; 
}

void compute_P1Qc(double lambda[], double n[], double theta, double phi, double **P1Q, double deltam, double m, double Delta[], double nu)
{
  int i,j;
  double *NN, *Delta_cross_NN;
  double NN_dot_n, 
         NN_dot_lambda;

  NN = dvector(1,3);
  Delta_cross_NN = dvector(1,3);

  /* get N_hat */
  NN[1] = sin(theta)*cos(phi);
  NN[2] = sin(theta)*sin(phi);
  NN[3] = cos(theta);

  /* Delta \cross N */
  Delta_cross_NN[1] = Delta[2]*NN[3] - NN[2]*Delta[3];
  Delta_cross_NN[2] = NN[1]*Delta[3] - Delta[1]*NN[3];
  Delta_cross_NN[3] = Delta[1]*NN[2] - NN[1]*Delta[2];

  /* set scalar products */
  NN_dot_n = NN[1]*n[1]+NN[2]*n[2]+NN[3]*n[3];
  NN_dot_lambda = NN[1]*lambda[1]+NN[2]*lambda[2]+NN[3]*lambda[3];

  /* Eqn (4.9c) of Kidder, PRD 52, 821 (1995) */
  for(i=1;i<=3;i++){
      for(j=1;j<=3;j++){
        P1Q[i][j] = 2.0/3.0*(1.0-3.0*nu) 
                    *(   pow(NN_dot_n,2)*( 5.0*n[i]*n[j]-7.0*lambda[i]*lambda[j] )
                       - 8.0*NN_dot_n*NN_dot_lambda*( n[i]*lambda[j] + n[j]*lambda[i] ) 
                       + pow(NN_dot_lambda,2)*( 3.0*lambda[i]*lambda[j] - n[i]*n[j] )              
                     )                                                                  
                    -(19.0-3.0*nu)/3.0*( n[i]*n[j]-lambda[i]*lambda[j] )                
                    +( n[i]*Delta_cross_NN[j]+n[j]*Delta_cross_NN[i])/m/m;   
      }
  } 
  free_dvector(NN,1,3);
  free_dvector(Delta_cross_NN,1,3);
  return;
}

void compute_P15Qc(double y[], double lambda[], double n[], double theta, double phi, double **P15Q, double deltam, double m1, double m2, double nu)
{
  int i,j;
  double m,
         NN_dot_n, 
         NN_dot_lambda,
         Ln_cap_dot_S2_p_Delta,
         Ln_cap_dot_S5_p_3Delta;
  double *NN,
         *S, 
         *Delta, 
         *Ln_cap, 
         *S5_p_3Delta,
         *S2_p_Delta,
         *S_p_Delta, 
         *S9_p_5Delta,
         *n_X_S_p_Delta,
         *lambda_X_S9_p_5Delta,
         *S_p_Delta_X_NN;  
   
  NN = dvector(1,3);
  S = dvector(1,3); 
  Delta = dvector(1,3);
  Ln_cap = dvector(1,3); 
  S5_p_3Delta = dvector(1,3);
  S2_p_Delta = dvector(1,3);
  S_p_Delta = dvector(1,3); 
  S9_p_5Delta = dvector(1,3);

  n_X_S_p_Delta = dvector(1,3);
  lambda_X_S9_p_5Delta = dvector(1,3);
  S_p_Delta_X_NN = dvector(1,3);

  /* get mass */ 
  m = m1 + m2;

  /* set Delta, S, Ln_hat */
  for(i=1;i<=3;i++){
    Delta[i] = m*(y[i+4]/m2 - y[i+1]/m1);
    S[i] = y[i+1] + y[i+4];
    Ln_cap[i] = y[i+7];
  }

  /* renormalize Ln_hat */
  for(i=1;i<=3;i++){
      Ln_cap[i] = Ln_cap[i]/sqrt(Ln_cap[1]*Ln_cap[1] + Ln_cap[2]*Ln_cap[2] + Ln_cap[3]*Ln_cap[3]);
  }

  /* get N_hat */
  NN[1] = sin(theta)*cos(phi);
  NN[2] = sin(theta)*sin(phi);
  NN[3] = cos(theta);

  /* set some sums of S and Delta */  
  for(i=1;i<=3;i++){
      S5_p_3Delta[i] = 5.0*S[i] + 3.0*deltam/m*Delta[i];
      S2_p_Delta[i]  = 2.0*S[i] +     deltam/m*Delta[i];
      S_p_Delta[i]   =     S[i] +     deltam/m*Delta[i];
      S9_p_5Delta[i] = 9.0*S[i] + 5.0*deltam/m*Delta[i];
  }

  /* set scalar products */
  NN_dot_n = NN[1]*n[1]+NN[2]*n[2]+NN[3]*n[3];
  NN_dot_lambda = NN[1]*lambda[1]+NN[2]*lambda[2]+NN[3]*lambda[3];

  /*  dot_product(Ln_cap, S5_p_3Delta) */
  Ln_cap_dot_S5_p_3Delta = Ln_cap[1]*S5_p_3Delta[1] 
                         + Ln_cap[2]*S5_p_3Delta[2] 
                         + Ln_cap[3]*S5_p_3Delta[3];

  /* dot_product(Ln_cap,S2_p_Delta) */
  Ln_cap_dot_S2_p_Delta = Ln_cap[1]*S2_p_Delta[1]
                        + Ln_cap[2]*S2_p_Delta[2]
                        + Ln_cap[3]*S2_p_Delta[3];

  /* cross_prod(n, S_p_Delta) =  n_X_S_p_Delta */
  n_X_S_p_Delta[1] = n[2]*S_p_Delta[3] - S_p_Delta[2]*n[3];
  n_X_S_p_Delta[2] = S_p_Delta[1]*n[3] - n[1]*S_p_Delta[3];
  n_X_S_p_Delta[3] = n[1]*S_p_Delta[2] - S_p_Delta[1]*n[2];    

  /*  cross_prod(lambda, S9_p_5Delta) =  lambda_X_S9_p_5Delta */
  lambda_X_S9_p_5Delta[1] = lambda[2]*S9_p_5Delta[3] - S9_p_5Delta[2]*lambda[3];
  lambda_X_S9_p_5Delta[2] = S9_p_5Delta[1]*lambda[3] - lambda[1]*S9_p_5Delta[3];
  lambda_X_S9_p_5Delta[3] = lambda[1]*S9_p_5Delta[2] - S9_p_5Delta[1]*lambda[2];

  /* cross_prod(S_p_Delta, NN) = S_p_Delta_X_NN */
  S_p_Delta_X_NN[1] = S_p_Delta[2]*NN[3] - NN[2]*S_p_Delta[3];
  S_p_Delta_X_NN[2] = NN[1]*S_p_Delta[3] - S_p_Delta[1]*NN[3];
  S_p_Delta_X_NN[3] = S_p_Delta[1]*NN[2] - NN[1]*S_p_Delta[2];

  /* Eqn (4.9d) of Kidder, PRD 52, 821 (1995) */
  for(i=1;i<=3;i++){
      for(j=1;j<=3;j++){
        P15Q[i][j] =  deltam/m*( 
                                (1.0-2.0*nu)*(
                                                 0.5*pow(NN_dot_lambda,3)*(n[i]*n[j]-4.0*lambda[i]*lambda[j])
                                               + 0.25*pow(NN_dot_n,2)*NN_dot_lambda*(58.0*lambda[i]*lambda[j]-37.0*n[i]*n[j])
                                               - 65.0/12.0*pow(NN_dot_n,3)*(n[i]*lambda[j]+n[j]*lambda[i])
                                               + 7.5*NN_dot_n*pow(NN_dot_lambda,2)*(n[i]*lambda[j]+n[j]*lambda[i])
                                             )
                                - NN_dot_lambda*(1.0/12.0*(101.0 -12.0*nu)*n[i]*n[j] - 0.5*(19.0-4.0*nu)*lambda[i]*lambda[j])
                                - 1.0/12.0*(149.0-36.0*nu)*(n[i]*lambda[j]+n[j]*lambda[i]) 
                              ) 
                    - 2.0/m/m*(  
                                lambda[i]*lambda[j]*Ln_cap_dot_S5_p_3Delta
                              - 6.0*n[i]*n[j]*Ln_cap_dot_S2_p_Delta
                              + (lambda[i]*n_X_S_p_Delta[j] + lambda[j]*n_X_S_p_Delta[i])
                              + 0.5*(n[i]*lambda_X_S9_p_5Delta[j] +n[j]*lambda_X_S9_p_5Delta[i])
                              + NN_dot_lambda*( S_p_Delta_X_NN[i]*n[j] + S_p_Delta_X_NN[j]*n[i] )
                              + NN_dot_n*( S_p_Delta_X_NN[i]*lambda[j] + S_p_Delta_X_NN[j]*lambda[i] )
                     );
      }
  }

  free_dvector(NN,1,3);
  free_dvector(S,1,3); 
  free_dvector(Delta,1,3);
  free_dvector(Ln_cap,1,3); 
  free_dvector(S5_p_3Delta,1,3);
  free_dvector(S2_p_Delta,1,3);
  free_dvector(S_p_Delta,1,3); 
  free_dvector(S9_p_5Delta,1,3); 

  free_dvector(n_X_S_p_Delta,1,3);
  free_dvector(lambda_X_S9_p_5Delta,1,3);
  free_dvector(S_p_Delta_X_NN,1,3);
}

/* Spin-Spin term (at 2PN order) */
void compute_P2QSSc(double y[], double n[], double theta, double phi, double **P2QSSc, double m1, double m2)
{
  int i,j;
  double m, mqube, mu,
         S1_dot_S2, n_dot_S1, n_dot_S2;
  double *S1, *S2;
  double rom3_P2QSS_ij;

  /* spin pointers */
  S1 = y+1;
  S2 = y+4;

  /* get mass and mu */
  m = m1 + m2;
  mu = m1*m2/m;
  mqube = m*m*m;

  /* set scalar products */
  S1_dot_S2 = S1[1]*S2[1]+S1[2]*S2[2]+S1[3]*S2[3];
  n_dot_S2  =  n[1]*S2[1]+ n[2]*S2[2]+ n[3]*S2[3];
  n_dot_S1  = S1[1]*n[1] +S1[2]*n[2] +S1[3]*n[3];

  /* Eqn (3.22b) of Kidder, PRD 52, 821 (1995) */
  for(i=1;i<=3;i++)
  {
    for(j=1;j<=3;j++)
    {
      /* (r/m)^3 P^2Q_{SS}^ij */
      rom3_P2QSS_ij = (-6.0/(mu*mqube))*
                      ( n[i]*n[j]*(S1_dot_S2-5.0*n_dot_S1*n_dot_S2) + 
                        (n[i]*S1[j]+n[j]*S1[i])*n_dot_S2 +
                        (n[i]*S2[j]+n[j]*S2[i])*n_dot_S1 
                      );
      P2QSSc[i][j] = rom3_P2QSS_ij;
    }
  }
}
