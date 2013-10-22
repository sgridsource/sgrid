/* Driver for routine odeint */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"

#define N 11  /* number of vars in odeint */
#define GammaE 0.577215664901532860606512090082402431042  /* Euler's gamma */

//int nrhs;   /* counts function evaluations */

/* global constants used throughout the program */
double nu, m, mu, deltam, m1, m2;
double c1, c2, c3, c4, c5, c6, c7, c8, c9, 
       c10, c11, c12, c13, c14, c15, c16,  
       c17, c18, c19, c20;
/* flags (0 or 1) that decide which order of v=om^{1/3} we include */
int f1, f2, f3, f4, f5, f6, f7;
double v2cut, cutoffP, v2trans, v6cut;
int OLS1, OLS2, OLS3, OSS1;
int EOM_type; /* types of EOM from PN_CircularOrbit_GWs.h */



/* compute global constants that depend only on masses or pi,gammaE,theta_cap */
void PN_CircOrbit_compute_constants(double m1_in, double m2_in)
{
  double theta_cap, pi, gammaE;

  pi = PI;
  gammaE = GammaE;
  theta_cap = 1039.0/4620.0;
  m1  = m1_in;
  m2  = m2_in;
  m   = m1 + m2;
  mu  = m1*m2/m;
  nu  = mu/m;
  deltam = m1 - m2;

  /* consts for omega evo */
  c1  = 96.0*nu/5.0;
  c2  = (743.0 + 924.0*nu)/336.0;
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

  /* consts for S1, S2, Ln evo */
  c15 = 0.5/m;
  c16 = 4.0 + 3.0*m2/m1;
  c17 = 1.0/m/m;
  c18 = 4.0 + 3.0*m1/m2;
  c19 = 3.0/nu/m/m;
  c20 = 0.5/m/m/m;

  /* set EOM_type from pars */
  if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","TaylorT4_bug")) 
    EOM_type=TaylorT4_bug;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","TaylorT4")) 
    EOM_type=TaylorT4;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","Kidder1995"))
    EOM_type=Kidder1995;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","Kidder1995_v2cut"))
    EOM_type=Kidder1995_v2cut;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","BuonannoEtAl2003_v2cut"))
    EOM_type=BuonannoEtAl2003_v2cut;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","Kidder1995_v6cut"))
    EOM_type=Kidder1995_v6cut;
  else if(Getv("PN_CircularOrbit_GWs_OrbitEOMtype","BuonannoEtAl2003_v6cut"))
    EOM_type=BuonannoEtAl2003_v6cut;
  else
    EOM_type=BuonannoEtAl2003;

  /* flags (0 or 1) that decide which order of v=om^{1/3} we include */
  f1 = Geti("PN_CircularOrbit_GWs_OrbitOv1"); /* order v */
  f2 = Geti("PN_CircularOrbit_GWs_OrbitOv2"); /* order v^2 */
  f3 = Geti("PN_CircularOrbit_GWs_OrbitOv3"); /* order v^3 */
  f4 = Geti("PN_CircularOrbit_GWs_OrbitOv4"); /* order v^4 */
  f5 = Geti("PN_CircularOrbit_GWs_OrbitOv5"); /* order v^5 */
  f6 = Geti("PN_CircularOrbit_GWs_OrbitOv6"); /* order v^6 */
  f7 = Geti("PN_CircularOrbit_GWs_OrbitOv7"); /* order v^7 */

  /* term we can put into EOM to stop omega from growing */
  v2cut = Getd("PN_CircularOrbit_GWs_Orbitv2cutoff");
  v2trans = Getd("PN_CircularOrbit_GWs_Orbitv2trans");
  cutoffP = Getd("PN_CircularOrbit_GWs_OrbittransP");
  v6cut = Getd("PN_CircularOrbit_GWs_Orbitv6cutoff");

  /* spin-orbit and spin-spin flags to decide which order in LS and SS
     we incluide in EOM */
  OLS1 = Geti("PN_CircularOrbit_GWs_OrbitOLS1");
  OLS2 = Geti("PN_CircularOrbit_GWs_OrbitOLS2");
  OLS3 = Geti("PN_CircularOrbit_GWs_OrbitOLS3");
  OSS1 = Geti("PN_CircularOrbit_GWs_OrbitOSS1");
}

/* function to be passed to odeintegrate */
void PN_CircOrbit_derivs(double x,double y[],double dydx[])
{
  int i;
  double om, v1, v2, v3, v4, v5, v6, v7, oov1;
  double V1, V2, V3, V4, V5, V6, V7;  /* powers of v multiplied by flags */
  double *Ln_cap, *S1, *S2, *LNS1, *LNS2, *LNS, Ln_cap_dot_S1, Ln_cap_dot_S2, S1_dot_S2, Ln_cap_mod;
  // nrhs++;

  LNS1 = dvector(1,3);
  LNS2 = dvector(1,3);
  LNS = dvector(1,3);

  om  = m*y[1];            // mw^1  /* here w is dimensionful orb. ang. vel. */
  v1  = pow(om,(1.0/3.0)); // mw^1/3
  v2  = v1*v1;             // mw^2/3
  v3  = om;                // mw
  v4  = v2*v2;             // mw^4/3 
  v5  = v4*v1;             // mw^5/3
  v6  = v3*v3;             // mw^2
  v7  = v6*v1;	           // mw^7/3
  oov1 = 1.0/v1;           // mw^-1/3

  /* set some powers of v to zero according to flags */
  V1 = v1*f1;
  V2 = v2*f2;
  V3 = v3*f3;
  V4 = v4*f4;
  V5 = v5*f5;
  V6 = v6*f6;
  V7 = v7*f7;

  /* pointers S1, S2, Ln_cap to correct place inside y */
  S1 = y+1;
  S2 = y+4;
  Ln_cap = y+7;

  /* renormalize Ln_cap and thus (also inside y) */
  Ln_cap_mod = sqrt(Ln_cap[1]*Ln_cap[1]+Ln_cap[2]*Ln_cap[2]+Ln_cap[3]*Ln_cap[3]);
  
  for(i=1;i<=3;i++)
      Ln_cap[i] = Ln_cap[i]/Ln_cap_mod;

  /* scalar products */
  Ln_cap_dot_S1 = Ln_cap[1]*S1[1] + Ln_cap[2]*S1[2] + Ln_cap[3]*S1[3];
  Ln_cap_dot_S2 = Ln_cap[1]*S2[1] + Ln_cap[2]*S2[2] + Ln_cap[3]*S2[3];
  S1_dot_S2     = S1[1]*S2[1] + S1[2]*S2[2] + S1[3]*S2[3];

  /* decide which eqn we use for freq evo */
  if(EOM_type==TaylorT4_bug)
  {
    /* stuff needed for Taylor T4 omega evo */
    /* products with 
       chis = 0.5*(chi1+chi2);  chia = 0.5*(chi1-chi2); */
    double Ln_cap_dot_chis, Ln_cap_dot_chia, Ln_cap_dot_chis2, 
           Ln_cap_dot_chis3, Ln_cap_dot_chia2, Ln_cap_dot_chia3, 
           S1_dot_S1, S2_dot_S2, chi1_dot_chi1, chi2_dot_chi2, 
           chis_dot_chis, chia_dot_chia, chis_dot_chia;

    Ln_cap_dot_chis = 0.5*(Ln_cap_dot_S1/(m1*m1) + Ln_cap_dot_S2/(m2*m2));
    Ln_cap_dot_chia = 0.5*(Ln_cap_dot_S1/(m1*m1) - Ln_cap_dot_S2/(m2*m2));
    Ln_cap_dot_chis2 = Ln_cap_dot_chis*Ln_cap_dot_chis;
    Ln_cap_dot_chis3 = Ln_cap_dot_chis2*Ln_cap_dot_chis;
    Ln_cap_dot_chia2 = Ln_cap_dot_chia*Ln_cap_dot_chia;
    Ln_cap_dot_chia3 = Ln_cap_dot_chia2*Ln_cap_dot_chia;
    S1_dot_S1     = S1[1]*S1[1] + S1[2]*S1[2] + S1[3]*S1[3];
    S2_dot_S2     = S2[1]*S2[1] + S2[2]*S2[2] + S2[3]*S2[3];
    chi1_dot_chi1 = S1_dot_S1/(m1*m1*m1*m1);
    chi2_dot_chi2 = S2_dot_S2/(m2*m2*m2*m2);
    chis_dot_chis = 0.25*( chi1_dot_chi1 
                          + 2*S1_dot_S2/(m1*m1*m2*m2) + chi2_dot_chi2 );
    chia_dot_chia = 0.25*( chi1_dot_chi1 
                          - 2*S1_dot_S2/(m1*m1*m2*m2) + chi2_dot_chi2 );
    chis_dot_chia = 0.25*( chi1_dot_chi1 - chi2_dot_chi2 );

    /* set some terms to zero (according to flags) */
    Ln_cap_dot_chis *= OLS1;
    Ln_cap_dot_chia *= OLS1;
    Ln_cap_dot_chis2 *= OLS2;
    Ln_cap_dot_chia2 *= OLS2;
    Ln_cap_dot_chis3 *= OLS3;
    Ln_cap_dot_chia3 *= OLS3;
    chis_dot_chis *= OSS1;
    chia_dot_chia *= OSS1;
    chis_dot_chia *= OSS1;

    /* Taylor T4 from Michael Boyle's SpEC old notebook with bug */
    dydx[1] = c1*y[1]*y[1]*v5*
        (1.0 + (-2.2113095238095237 - (11*nu)/4.)*V2 + 
          ((48*PI - 113*(deltam/m)*Ln_cap_dot_chia - 113*Ln_cap_dot_chis)/12. + 
             (19*nu*Ln_cap_dot_chis)/3.)*V3 + 
          ((59*nu*nu)/18. + 
             (34103 - 44037*chia_dot_chia + 
                189*(719*(Ln_cap_dot_chia2 + 
                      2*(deltam/m)*Ln_cap_dot_chia*Ln_cap_dot_chis + 
                      Ln_cap_dot_chis2) - 
                   466*(deltam/m)*chis_dot_chia) - 
                44037*chis_dot_chis)/18144. + 
             (nu*(13661 - 72618*Ln_cap_dot_chia2 + 
                  24486*chia_dot_chia + 12222*Ln_cap_dot_chis2 - 
                  4914*chis_dot_chis))/2016.)*V4 + 
          ((-12477*PI - 62638*(deltam/m)*Ln_cap_dot_chia - 62638*Ln_cap_dot_chis)/
              2016. - (79*nu*nu*Ln_cap_dot_chis)/3. + 
             (nu*(-11907*PI + 24339*(deltam/m)*Ln_cap_dot_chia + 
                  45950*Ln_cap_dot_chis))/504.)*V5 + 
          (117.72574285227559 - (1712*GammaE)/105. + (16*PI*PI)/3. - 
             (5605*nu*nu*nu)/2592. - (145*chia_dot_chia)/448. + 
             (5*((1035 + 50624*(deltam/m)*(deltam/m))*Ln_cap_dot_chia2 - 
                  21504*PI*Ln_cap_dot_chis + 51659*Ln_cap_dot_chis2 + 
                  2*(deltam/m)*Ln_cap_dot_chia*
                   (-10752*PI + 51659*Ln_cap_dot_chis) - 
                  522*(deltam/m)*chis_dot_chia))/4032. + 
             (nu*nu*(4869 + 350756*Ln_cap_dot_chia2 - 
                  116732*chia_dot_chia + 
                  178388*Ln_cap_dot_chis2 - 3276*chis_dot_chis))/
              8064. - (145*chis_dot_chis)/448. + 
             (nu*(-56198689 + 2045736*PI*PI + 
                  1888002*chia_dot_chia + 2903040*PI*Ln_cap_dot_chis - 
                  54*(102229*Ln_cap_dot_chia2 + 
                     386526*(deltam/m)*Ln_cap_dot_chia*Ln_cap_dot_chis + 
                     304997*Ln_cap_dot_chis2 - 
                     30002*(deltam/m)*chis_dot_chia) + 
                  13986*chis_dot_chis))/217728. - (856*log(16))/105. - 
             (1712*log(v1))/105.)*V6 + 
          ((2045*nu*nu*nu*Ln_cap_dot_chis)/216. + 
             (nu*nu*(365980*PI - 291109*(deltam/m)*Ln_cap_dot_chia - 
                  1415652*Ln_cap_dot_chia2*Ln_cap_dot_chis + 
                  3*Ln_cap_dot_chis*
                   (-398269 + 158228*chia_dot_chia + 
                     40740*Ln_cap_dot_chis2 - 16380*chis_dot_chis
                     )))/6048. + (nu*
                (3228075*PI + 10713953*Ln_cap_dot_chis + 
                  9*(2477496*(deltam/m)*Ln_cap_dot_chia3 + 
                     126*Ln_cap_dot_chia2*
                      (-3456*PI + 22229*Ln_cap_dot_chis) + 
                     42*(27*chia_dot_chia*
                         (128*PI - 827*Ln_cap_dot_chis) + 
                        Ln_cap_dot_chis*
                         (1879*Ln_cap_dot_chis2 - 
                           5066*(deltam/m)*chis_dot_chia - 
                           193*chis_dot_chis)) + 
                     (deltam/m)*Ln_cap_dot_chia*
                      (843811 - 831432*chia_dot_chia + 
                        402276*Ln_cap_dot_chis2 + 
                        98280*chis_dot_chis))))/54432. + 
             (-119205*PI - 10076804*Ln_cap_dot_chis + 
                4*(-2512188*(deltam/m)*Ln_cap_dot_chia3 + 
                   756*Ln_cap_dot_chia2*
                    (648*PI - 3323*(1 + 2*(deltam/m)*(deltam/m))*Ln_cap_dot_chis) + 
                   (deltam/m)*Ln_cap_dot_chia*
                    (-2519201 + 824796*chia_dot_chia + 
                      2268*(432*PI - 3323*Ln_cap_dot_chis)*Ln_cap_dot_chis + 
                      1649592*(deltam/m)*chis_dot_chia + 
                      824796*chis_dot_chis) - 
                   756*(chia_dot_chia*(216*PI - 1091*Ln_cap_dot_chis) - 
                      648*PI*Ln_cap_dot_chis2 + 
                      3323*Ln_cap_dot_chis3 + 
                      216*PI*(2*(deltam/m)*chis_dot_chia + 
                         chis_dot_chis) - 
                      1091*Ln_cap_dot_chis*
                       (2*(deltam/m)*chis_dot_chia + chis_dot_chis))
                   ))/108864.)*V7);
  }
  else if(EOM_type==TaylorT4)
  {
    /* stuff needed for Taylor T4 omega evo */
    /* products with 
       chis = 0.5*(chi1+chi2);  chia = 0.5*(chi1-chi2); */
    double Ln_cap_dot_chis, Ln_cap_dot_chia, Ln_cap_dot_chis2, 
           Ln_cap_dot_chis3, Ln_cap_dot_chia2, Ln_cap_dot_chia3, 
           S1_dot_S1, S2_dot_S2, chi1_dot_chi1, chi2_dot_chi2, 
           chis_dot_chis, chia_dot_chia, chis_dot_chia;

    Ln_cap_dot_chis = 0.5*(Ln_cap_dot_S1/(m1*m1) + Ln_cap_dot_S2/(m2*m2));
    Ln_cap_dot_chia = 0.5*(Ln_cap_dot_S1/(m1*m1) - Ln_cap_dot_S2/(m2*m2));
    Ln_cap_dot_chis2 = Ln_cap_dot_chis*Ln_cap_dot_chis;
    Ln_cap_dot_chis3 = Ln_cap_dot_chis2*Ln_cap_dot_chis;
    Ln_cap_dot_chia2 = Ln_cap_dot_chia*Ln_cap_dot_chia;
    Ln_cap_dot_chia3 = Ln_cap_dot_chia2*Ln_cap_dot_chia;
    S1_dot_S1     = S1[1]*S1[1] + S1[2]*S1[2] + S1[3]*S1[3];
    S2_dot_S2     = S2[1]*S2[1] + S2[2]*S2[2] + S2[3]*S2[3];
    chi1_dot_chi1 = S1_dot_S1/(m1*m1*m1*m1);
    chi2_dot_chi2 = S2_dot_S2/(m2*m2*m2*m2);
    chis_dot_chis = 0.25*( chi1_dot_chi1 
                          + 2*S1_dot_S2/(m1*m1*m2*m2) + chi2_dot_chi2 );
    chia_dot_chia = 0.25*( chi1_dot_chi1 
                          - 2*S1_dot_S2/(m1*m1*m2*m2) + chi2_dot_chi2 );
    chis_dot_chia = 0.25*( chi1_dot_chi1 - chi2_dot_chi2 );

    /* set some terms to zero (according to flags) */
    Ln_cap_dot_chis *= OLS1;
    Ln_cap_dot_chia *= OLS1;
    Ln_cap_dot_chis2 *= OLS2;
    Ln_cap_dot_chia2 *= OLS2;
    Ln_cap_dot_chis3 *= OLS3;
    Ln_cap_dot_chia3 *= OLS3;
    chis_dot_chis *= OSS1;
    chia_dot_chia *= OSS1;
    chis_dot_chia *= OSS1;

    /* Taylor T4 from Michael Boyle's SpEC notebook */
    dydx[1] = c1*y[1]*y[1]*v5*
      (1 + (-2.2113095238095237 - (11*nu)/4.)*V2 + 
        V3*((48*PI - 113*(deltam/m)*Ln_cap_dot_chia - 113*Ln_cap_dot_chis)/
            12. + (19*nu*Ln_cap_dot_chis)/3.) + 
        V4*((59*nu*nu)/18. + 
           nu*(10*chia_dot_chia - 30*Ln_cap_dot_chia2 + 
              (13661 - 588*chis_dot_chis + 84*Ln_cap_dot_chis2)/2016.)\
            + (34103 - 44037*chia_dot_chia - 44037*chis_dot_chis + 
              189*(-466*(deltam/m)*chis_dot_chia + 
                 719*(Ln_cap_dot_chia2 + 
                    2*(deltam/m)*Ln_cap_dot_chia*Ln_cap_dot_chis + Ln_cap_dot_chis2)
                 ))/18144.) +
        V5*((-12477*PI - 62638*(deltam/m)*Ln_cap_dot_chia - 
              62638*Ln_cap_dot_chis)/2016. - (79*nu*nu*Ln_cap_dot_chis)/3. + 
           (nu*(-11907*PI + 24339*(deltam/m)*Ln_cap_dot_chia + 45950*Ln_cap_dot_chis))/
            504.) +
        V6*(117.72574285227559 - (1712*GammaE)/105. - 
           (5605*nu*nu*nu)/2592. + (16*PI*PI)/3. - 
           (145*chia_dot_chia)/448. - (145*chis_dot_chis)/448. + 
           nu*nu*(0.6037946428571429 - (89*chia_dot_chia)/6. + 
              (89*Ln_cap_dot_chia2)/2. - (7*chis_dot_chis)/144. + 
              (3041*Ln_cap_dot_chis2)/144.) + 
           (5*((1035 + 50624*(deltam/m)*(deltam/m))*Ln_cap_dot_chia2 - 
                522*(deltam/m)*chis_dot_chia - 21504*PI*Ln_cap_dot_chis + 
                51659*Ln_cap_dot_chis2 + 
                2*(deltam/m)*Ln_cap_dot_chia*(-10752*PI + 51659*Ln_cap_dot_chis)))/4032.\
            + (nu*(-56198689 + 2045736*PI*PI + 1187190*chia_dot_chia + 
                714798*chis_dot_chis + 2903040*PI*Ln_cap_dot_chis - 
                54*(65815*Ln_cap_dot_chia2 - 
                   30002*(deltam/m)*chis_dot_chia + 
                   386526*(deltam/m)*Ln_cap_dot_chia*Ln_cap_dot_chis + 
                   341411*Ln_cap_dot_chis2)))/217728. - 
           (856*log(16))/105. - (1712*log(v1))/105.) +
        V7*((2045*nu*nu*nu*Ln_cap_dot_chis)/216. + 
           (nu*nu*(365980*PI - 291109*(deltam/m)*Ln_cap_dot_chia - 
                1294272*Ln_cap_dot_chia2*Ln_cap_dot_chis + 
                3*(-398269 + 143808*chia_dot_chia - 1960*chis_dot_chis)*
                 Ln_cap_dot_chis + 840*Ln_cap_dot_chis3))/6048. + 
           nu*((739*(deltam/m)*Ln_cap_dot_chia3)/2. + 
              chia_dot_chia*(24*PI - (20269*Ln_cap_dot_chis)/144.) + 
              Ln_cap_dot_chia2*(-72*PI + (60907*Ln_cap_dot_chis)/144.) + 
              ((deltam/m)*Ln_cap_dot_chia*
                 (843811 - 744912*chia_dot_chia + 11760*chis_dot_chis + 
                   645036*Ln_cap_dot_chis2))/6048. + 
              (3228075*PI + Ln_cap_dot_chis*
                  (10713953 - 1914948*(deltam/m)*chis_dot_chia - 
                    851634*chis_dot_chis + 2895102*Ln_cap_dot_chis2))/
               54432.) + (-119205*PI - 10076804*Ln_cap_dot_chis + 
              4*(-2512188*(deltam/m)*Ln_cap_dot_chia3 + 
                 756*Ln_cap_dot_chia2*
                  (648*PI - 3323*(1 + 2*(deltam/m)*(deltam/m))*Ln_cap_dot_chis) + 
                 (deltam/m)*Ln_cap_dot_chia*
                  (-2519201 + 824796*chia_dot_chia + 
                    1649592*(deltam/m)*chis_dot_chia + 824796*chis_dot_chis + 
                    2268*(432*PI - 3323*Ln_cap_dot_chis)*Ln_cap_dot_chis) - 
                 756*(216*PI*(2*(deltam/m)*chis_dot_chia + chis_dot_chis) + 
                    chia_dot_chia*(216*PI - 1091*Ln_cap_dot_chis) - 
                    1091*(2*(deltam/m)*chis_dot_chia + chis_dot_chis)*
                     Ln_cap_dot_chis - 648*PI*Ln_cap_dot_chis2 + 
                    3323*Ln_cap_dot_chis3)))/108864.) 
      );
  }
  else if(EOM_type==Kidder1995)
  {
    /* Evo of Omega */
    /* Eqn (4.14) Kidder, PRD 52, 821 (1995) */
    dydx[1] = c1*y[1]*y[1]*v5*(1.0 
                               - c2*V2 
                               - (c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5)*V3
                               + c6*V4 * 0
                               - c7*S1_dot_S2*V4 
                               + c8*Ln_cap_dot_S1*Ln_cap_dot_S2*V4 );
  }
  else if(EOM_type==Kidder1995_v2cut)
  {
    double B2 = 1.0/(v2cut*v2cut);
    double a2 = -c2 + B2;
    double a3 = -(c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5);
    double a4 = c6*0 - c7*S1_dot_S2 + c8*Ln_cap_dot_S1*Ln_cap_dot_S2 + a2*B2;
    /* Evo of Omega, Kidder1995 with v2cut */
    dydx[1] = c1*y[1]*y[1]*v5*(1 - B2*V2)*
              (1.0 + a2*V2 + a3*V3 + a4*V4);
  }
  else if(EOM_type==BuonannoEtAl2003_v2cut)
  {
    double B2 = 1.0/(v2cut*v2cut);
    double a2 = -c2 + B2;
    double a3 = -(c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5);
    double a4 = c6 - c7*S1_dot_S2 + c8*Ln_cap_dot_S1*Ln_cap_dot_S2 + a2*B2;
    double a5 = -c9 + a3*B2;
    double a6 = c10 + c11 + c12 - c13*log(16.0*v2) + a4*B2;
    double a7 = c14 + a5*B2;
    /* Note: a2*V2 + a3*V3 + a4*V4 + a5*V5 + a6*V6 + a7*V7 
              = tanh(e2*V2 + e3*V3 + e4*V4 + e5*V5 + e6*V6 + e7*V7)
       if: */
    double e2 = a2;
    double e3 = a3;
    double e4 = a4;
    double e5 = a5;
    double e6 = a6 + e2*e2*e2/3.0;
    double e7 = a7 + e2*e2*e3;
    double AA = cutoffP;
    double BB = 1.0/(v2trans*v2trans);
    double CC = BB*BB*BB*BB*BB;
    /* Evo of Omega, BuonannoEtAl2003_v2cut */
    dydx[1] = c1*y[1]*y[1]*v5*(1 - B2*V2)*
              (1.0 + tanh(e2*V2 + e3*V3 + e4*V4 + e5*V5 + e6*V6 + e7*V7) +
               AA*(1-BB*v2)*(1-BB*v2) * CC*v6*v2*V2);
  }
  else if(EOM_type==Kidder1995_v6cut)
  {
    double B6 = 1.0/pow(v6cut,6.0);
    double a2 = -c2;
    double a3 = -(c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5);
    double a4 = c6*0 - c7*S1_dot_S2 + c8*Ln_cap_dot_S1*Ln_cap_dot_S2;
    /* Evo of Omega, Kidder1995 with v6cut */
    dydx[1] = c1*y[1]*y[1]*v5*(1 - B6*V6)*
              (1.0 + a2*V2 + a3*V3 + a4*V4);
  }
  else if(EOM_type==BuonannoEtAl2003_v6cut)
  {
    double B6 = 1.0/pow(v6cut,6.0);
    double a2 = -c2;
    double a3 = -(c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5);
    double a4 = c6 - c7*S1_dot_S2 + c8*Ln_cap_dot_S1*Ln_cap_dot_S2;
    double a5 = -c9;
    double a6 = c10 + c11 + c12 - c13*log(16.0*v2) + B6;
    double a7 = c14;
    /* Evo of Omega, BuonannoEtAl2003_v6cut */
    dydx[1] = c1*y[1]*y[1]*v5*(1 - B6*V6)*
              (1.0 + a2*V2 + a3*V3 + a4*V4 + a5*V5 + a6*V6 + a7*V7);
  }
  else
  {
    /* Evo of Omega */
    /* Eqn (1) of Erratum in Buonanno et al, PRD 74, 029904(E) (2006) */
    dydx[1] = c1*y[1]*y[1]*v5*(1.0 
                               - c2*V2 
                               - (c3*Ln_cap_dot_S1 + c4*Ln_cap_dot_S2 - c5)*V3
                               + c6*V4
                               - c7*S1_dot_S2*V4 
                               + c8*Ln_cap_dot_S1*Ln_cap_dot_S2*V4 
                               - c9*V5     
                               + (c10 + c11 + c12 - c13*log(16.0*v2))*V6 
                               + c14*V7 );
  }

  /* terms in Eqn (2,3,9) of Buonanno et al, PRD 67, 104025 (2003) */
  for(i=1;i<=3;i++)
  {
    LNS1[i] = c15*( nu*c16*oov1*Ln_cap[i] + c17*(S2[i] - 3.0*Ln_cap_dot_S2*Ln_cap[i]) );
    LNS2[i] = c15*( nu*c18*oov1*Ln_cap[i] + c17*(S1[i] - 3.0*Ln_cap_dot_S1*Ln_cap[i]) );
    LNS[i]  = c16*S1[i] + c18*S2[i] - c19*v1*( S1[i]*Ln_cap_dot_S2 + S2[i]*Ln_cap_dot_S1 );
  } 

  /* Evo of Spin1, Eqn (2), Buonanno et al, PRD 67, 104025 (2003) */	 
  dydx[2] = v6*( LNS1[2]*S1[3] - S1[2]*LNS1[3]); 
  dydx[3] = v6*(-LNS1[1]*S1[3] + S1[1]*LNS1[3]);
  dydx[4] = v6*( LNS1[1]*S1[2] - S1[1]*LNS1[2]);

  /* Evo of Spin2, Eqn (3), Buonanno et al, PRD 67, 104025 (2003) */	 
  dydx[5] = v6*( LNS2[2]*S2[3] - S2[2]*LNS2[3]);
  dydx[6] = v6*(-LNS2[1]*S2[3] + S2[1]*LNS2[3]); 
  dydx[7] = v6*( LNS2[1]*S2[2] - S2[1]*LNS2[2]); 

  /* Evo of \hat{L}_N = Ln_cap, Eqn (9), Buonanno et al, PRD 67, 104025 (2003) */
  dydx[8]  = c20*v6*( LNS[2]*Ln_cap[3] - Ln_cap[2]*LNS[3]);
  dydx[9]  = c20*v6*(-LNS[1]*Ln_cap[3] + Ln_cap[1]*LNS[3]); 
  dydx[10] = c20*v6*( LNS[1]*Ln_cap[2] - Ln_cap[1]*LNS[2]);

  /* phase evo with respect to ascending node defined in zx-plane */ 
  if((Ln_cap[1]==0.0)&&(Ln_cap[3]==0.0))
      dydx[11] = y[1];
  else
    dydx[11] = y[1] - Ln_cap[2]*(Ln_cap[3]*dydx[8] - Ln_cap[1]*dydx[10])/(Ln_cap[1]*Ln_cap[1] + Ln_cap[3]*Ln_cap[3]);

  free_dvector(LNS1,1,3);
  free_dvector(LNS2,1,3);
  free_dvector(LNS,1,3);  		 
}

/* get radius of orbit for a given y[], assuming PN_CircOrbit_compute_constants
   has been called before */
double PN_CircOrbit_compute_r(double y[])
{
  double om, v2, v3, v4, oov2, LnS1, LnS2, SS, r;
  double V2, V3, V4;  /* powers of v multiplied by flags */
  om = m*y[1];
  v2 = pow(om,(2.0/3.0));
  v3 = om;
  v4 = v2*v2;
  oov2 = 1.0/v2;

  /* set some powers of v to zero according to flags */
  // V1 = v1*f1;
  V2 = v2*f2;
  V3 = v3*f3;
  V4 = v4*f4;
  // V5 = v5*f5;
  // V6 = v6*f6;
  // V7 = v7*f7;

  /* scalar products */
  LnS1 = (y[8]*y[2] + y[9]*y[3] + y[10]*y[4]);
  LnS2 = (y[8]*y[5] + y[9]*y[6] + y[10]*y[7]);
  SS = (y[2]*y[5] + y[3]*y[6] + y[4]*y[7]);

  /* from Kidder PRD 52, 821 (1995), Eq (4.13). I.e. only up to v^4 */
  r  = m*oov2*( 1.0 
               - (3.0 - nu)*V2/3.0 
               - V3*( LnS1*(2.0*m1*m1/m/m+3.0*nu)/m1/m1 
                     +LnS2*(2.0*m2*m2/m/m+3.0*nu)/m2/m2 )/3.0
               + ( nu*(19.0/4.0 + nu/9.0) 
                   - 0.5*nu*(SS - 3.0*LnS1*LnS2)/m1/m1/m2/m2 )*V4);
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
