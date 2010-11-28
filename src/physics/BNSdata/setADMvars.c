/* setADMvars.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 27.11.2010 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void setADMvars(tGrid *grid)
{
int corot1 = Getv("BNSdata_rotationstate1","corotation");
int corot2 = Getv("BNSdata_rotationstate2","corotation");
double n = Getd("BNSdata_n");
double kappa = Getd("BNSdata_kappa");
double Omega = Getd("BNSdata_Omega");
double xCM = Getd("BNSdata_x_CM");

int bi;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;



int index_BNSdata_Psi = Ind("BNSdata_Psi");
double *Psi = box->v[index_BNSdata_Psi + 0];
int index_BNSdata_Bx = Ind("BNSdata_Bx");
double *B1 = box->v[index_BNSdata_Bx + 0];
double *B2 = box->v[index_BNSdata_Bx + 1];
double *B3 = box->v[index_BNSdata_Bx + 2];
int index_BNSdata_alphaP = Ind("BNSdata_alphaP");
double *alphaP = box->v[index_BNSdata_alphaP + 0];
int index_BNSdata_Sigma = Ind("BNSdata_Sigma");
double *Sigma = box->v[index_BNSdata_Sigma + 0];
int index_BNSdata_q = Ind("BNSdata_q");
double *q = box->v[index_BNSdata_q + 0];
int index_BNSdata_wBx = Ind("BNSdata_wBx");
double *wB1 = box->v[index_BNSdata_wBx + 0];
double *wB2 = box->v[index_BNSdata_wBx + 1];
double *wB3 = box->v[index_BNSdata_wBx + 2];
int index_BNSdata_Bxx = Ind("BNSdata_Bxx");
double *dB11 = box->v[index_BNSdata_Bxx + 0];
double *dB12 = box->v[index_BNSdata_Bxx + 1];
double *dB13 = box->v[index_BNSdata_Bxx + 2];
double *dB21 = box->v[index_BNSdata_Bxx + 3];
double *dB22 = box->v[index_BNSdata_Bxx + 4];
double *dB23 = box->v[index_BNSdata_Bxx + 5];
double *dB31 = box->v[index_BNSdata_Bxx + 6];
double *dB32 = box->v[index_BNSdata_Bxx + 7];
double *dB33 = box->v[index_BNSdata_Bxx + 8];
int index_BNSdata_Sigmax = Ind("BNSdata_Sigmax");
double *dSigma1 = box->v[index_BNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_BNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_BNSdata_Sigmax + 2];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_gxx = Ind("gxx");
double *g11 = box->v[index_gxx + 0];
double *g12 = box->v[index_gxx + 1];
double *g13 = box->v[index_gxx + 2];
double *g22 = box->v[index_gxx + 3];
double *g23 = box->v[index_gxx + 4];
double *g33 = box->v[index_gxx + 5];
int index_psi = Ind("psi");
double *psi = box->v[index_psi + 0];
int index_Kxx = Ind("Kxx");
double *K11 = box->v[index_Kxx + 0];
double *K12 = box->v[index_Kxx + 1];
double *K13 = box->v[index_Kxx + 2];
double *K22 = box->v[index_Kxx + 3];
double *K23 = box->v[index_Kxx + 4];
double *K33 = box->v[index_Kxx + 5];
int index_alpha = Ind("alpha");
double *alpha = box->v[index_alpha + 0];
int index_betax = Ind("betax");
double *beta1 = box->v[index_betax + 0];
double *beta2 = box->v[index_betax + 1];
double *beta3 = box->v[index_betax + 2];
int index_rho = Ind("rho");
double *rho = box->v[index_rho + 0];
int index_jx = Ind("jx");
double *jdo1 = box->v[index_jx + 0];
double *jdo2 = box->v[index_jx + 1];
double *jdo3 = box->v[index_jx + 2];
int index_Sxx = Ind("Sxx");
double *Sdo11 = box->v[index_Sxx + 0];
double *Sdo12 = box->v[index_Sxx + 1];
double *Sdo13 = box->v[index_Sxx + 2];
double *Sdo22 = box->v[index_Sxx + 3];
double *Sdo23 = box->v[index_Sxx + 4];
double *Sdo33 = box->v[index_Sxx + 5];


double alpha2;
double dSigmaUp1;
double dSigmaUp2;
double dSigmaUp3;
double gdB;
double h;
double h2;
double jup1;
double jup2;
double jup3;
double LB11;
double LB12;
double LB13;
double LB22;
double LB23;
double LB33;
double LBdo11;
double LBdo12;
double LBdo13;
double LBdo21;
double LBdo22;
double LBdo23;
double LBdo31;
double LBdo32;
double LBdo33;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double P;
double Psi2;
double Psi4;
double Psim2;
double Psim4;
double Psim6;
double rho0;
double rhoE;
double uzerosqr;
double vR1;
double vR2;
double vR3;
double vRplusbetado1;
double vRplusbetado2;
double vRplusbetado3;
double w1;
double w2;
double w3;
double wBDown1;
double wBDown2;
double wBDown3;
double wDown1;
double wDown2;
double wDown3;



/* Jetzt geht's los! */

FirstDerivsOf_Sa(box, Ind("BNSdata_Bx"),       Ind("BNSdata_Bxx")); 


FirstDerivsOf_S(box, Ind("BNSdata_Sigma"),     Ind("BNSdata_Sigmax")); 


forallpoints(box, ijk) { 

alpha[ijk]
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha[ijk])
;

Psi2
=
pow2(Psi[ijk])
;

Psi4
=
pow2(Psi2)
;

OmegaCrossR1
=
-(Omega*y[ijk])
;

OmegaCrossR2
=
Omega*(-xCM + x[ijk])
;

OmegaCrossR3
=
0
;

beta1[ijk]
=
OmegaCrossR1 + B1[ijk]
;

beta2[ijk]
=
OmegaCrossR2 + B2[ijk]
;

beta3[ijk]
=
OmegaCrossR3 + B3[ijk]
;

gdB
=
dB11[ijk] + dB22[ijk] + dB33[ijk]
;

LB11
=
-0.66666666666666666667*gdB + 2.*dB11[ijk]
;

LB12
=
dB12[ijk] + dB21[ijk]
;

LB13
=
dB13[ijk] + dB31[ijk]
;

LB22
=
-0.66666666666666666667*gdB + 2.*dB22[ijk]
;

LB23
=
dB23[ijk] + dB32[ijk]
;

LB33
=
-0.66666666666666666667*gdB + 2.*dB33[ijk]
;

LBdo11
=
LB11
;

LBdo12
=
LB12
;

LBdo13
=
LB13
;

LBdo21
=
LB12
;

LBdo22
=
LB22
;

LBdo23
=
LB23
;

LBdo31
=
LB13
;

LBdo32
=
LB23
;

LBdo33
=
LB33
;

psi[ijk]
=
1.
;

g11[ijk]
=
Psi4
;

g12[ijk]
=
0
;

g13[ijk]
=
0
;

g22[ijk]
=
Psi4
;

g23[ijk]
=
0
;

g33[ijk]
=
Psi4
;

K11[ijk]
=
(0.5*LBdo11*Psi4)/alpha[ijk]
;

K12[ijk]
=
(0.5*LBdo12*Psi4)/alpha[ijk]
;

K12[ijk]
=
(0.5*LBdo21*Psi4)/alpha[ijk]
;

K13[ijk]
=
(0.5*LBdo13*Psi4)/alpha[ijk]
;

K13[ijk]
=
(0.5*LBdo31*Psi4)/alpha[ijk]
;

K22[ijk]
=
(0.5*LBdo22*Psi4)/alpha[ijk]
;

K23[ijk]
=
(0.5*LBdo23*Psi4)/alpha[ijk]
;

K23[ijk]
=
(0.5*LBdo32*Psi4)/alpha[ijk]
;

K33[ijk]
=
(0.5*LBdo33*Psi4)/alpha[ijk]
;



/* conditional */
if ((bi <= 1 || bi == 5) && corot1 || bi >= 2 && bi <= 4 && corot2) {

vR1
=
0
;

vR2
=
0
;

vR3
=
0
;

oouzerosqr
=
alpha2 - Psi4*(2.*(vR1*beta1[ijk] + vR2*beta2[ijk] + vR3*beta3[ijk]) + 
     pow2(vR1) + pow2(vR2) + pow2(vR3) + pow2(beta1[ijk]) + 
     pow2(beta2[ijk]) + pow2(beta3[ijk]))
;

uzerosqr
=
1./oouzerosqr
;


} else { /* if (!(bi <= 1 || bi == 5) && corot1 || bi >= 2 && bi <= 4 && corot2) */

Psim2
=
1/Psi2
;

Psim4
=
pow2(Psim2)
;

Psim6
=
Psim2*Psim4
;

dSigmaUp1
=
Psim4*dSigma1[ijk]
;

dSigmaUp2
=
Psim4*dSigma2[ijk]
;

dSigmaUp3
=
Psim4*dSigma3[ijk]
;

w1
=
Psim6*wB1[ijk]
;

w2
=
Psim6*wB2[ijk]
;

w3
=
Psim6*wB3[ijk]
;

wBDown1
=
wB1[ijk]
;

wBDown2
=
wB2[ijk]
;

wBDown3
=
wB3[ijk]
;

wDown1
=
Psim2*wBDown1
;

wDown2
=
Psim2*wBDown2
;

wDown3
=
Psim2*wBDown3
;

h
=
1. + q[ijk] + n*q[ijk]
;

h2
=
pow2(h)
;

uzerosqr
=
1/alpha2 + ((dSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
     (dSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
     (dSigmaUp3 + w3)*(wDown3 + dSigma3[ijk]))/(alpha2*h2)
;

}
/* if ((bi <= 1 || bi == 5) && corot1 || bi >= 2 && bi <= 4 && corot2) */


rho0
=
Power(q[ijk]/kappa,n)
;

P
=
rho0*q[ijk]
;

rhoE
=
rho0 + n*rho0*q[ijk]
;



/* conditional */
if (q[ijk] == 0) {

uzerosqr
=
0
;

}
/* if (q[ijk] == 0) */


rho[ijk]
=
alpha2*rhoE*uzerosqr + P*(-1. + alpha2*uzerosqr)
;

jup1
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR1 + beta1[ijk])
;

jup2
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR2 + beta2[ijk])
;

jup3
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR3 + beta3[ijk])
;

jdo1[ijk]
=
jup1*g11[ijk] + jup2*g12[ijk] + jup3*g13[ijk]
;

jdo2[ijk]
=
jup1*g12[ijk] + jup2*g22[ijk] + jup3*g23[ijk]
;

jdo3[ijk]
=
jup1*g13[ijk] + jup2*g23[ijk] + jup3*g33[ijk]
;

vRplusbetado1
=
(vR1 + beta1[ijk])*g11[ijk] + (vR2 + beta2[ijk])*g12[ijk] + 
  (vR3 + beta3[ijk])*g13[ijk]
;

vRplusbetado2
=
(vR1 + beta1[ijk])*g12[ijk] + (vR2 + beta2[ijk])*g22[ijk] + 
  (vR3 + beta3[ijk])*g23[ijk]
;

vRplusbetado3
=
(vR1 + beta1[ijk])*g13[ijk] + (vR2 + beta2[ijk])*g23[ijk] + 
  (vR3 + beta3[ijk])*g33[ijk]
;

Sdo11[ijk]
=
P*g11[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado1)
;

Sdo12[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado1*vRplusbetado2 + P*g12[ijk]
;

Sdo13[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado1*vRplusbetado3 + P*g13[ijk]
;

Sdo22[ijk]
=
P*g22[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado2)
;

Sdo23[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado2*vRplusbetado3 + P*g23[ijk]
;

Sdo33[ijk]
=
P*g33[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado3)
;


} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* setADMvars.c */
/* nvars = 51, n* = 180,  n/ = 40,  n+ = 184, n = 404, O = 1 */
