/* BNS_set_J_ADM_VolInt_integrand.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 27.11.2010 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BNS_set_J_ADM_VolInt_integrand(tGrid *grid, int iInteg)
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



int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
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
int index_BNSdata_Sigmax = Ind("BNSdata_Sigmax");
double *dSigma1 = box->v[index_BNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_BNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_BNSdata_Sigmax + 2];
double *Integ = box->v[iInteg+0];


double alpha;
double alpha2;
double beta1;
double beta2;
double beta3;
double dSigmaUp1;
double dSigmaUp2;
double dSigmaUp3;
double h;
double h2;
double jup1;
double jup2;
double jup3;
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

FirstDerivsOf_S(box, Ind("BNSdata_Sigma"),     Ind("BNSdata_Sigmax")); 


forallpoints(box, ijk) { 

alpha
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha)
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

beta1
=
OmegaCrossR1 + B1[ijk]
;

beta2
=
OmegaCrossR2 + B2[ijk]
;

beta3
=
OmegaCrossR3 + B3[ijk]
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
alpha2 - Psi4*(2.*(beta1*vR1 + beta2*vR2 + beta3*vR3) + pow2(beta1) + 
     pow2(beta2) + pow2(beta3) + pow2(vR1) + pow2(vR2) + pow2(vR3))
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


jup1
=
alpha*(P + rhoE)*uzerosqr*(beta1 + vR1)
;

jup2
=
alpha*(P + rhoE)*uzerosqr*(beta2 + vR2)
;

jup3
=
alpha*(P + rhoE)*uzerosqr*(beta3 + vR3)
;

Integ[ijk]
=
Psi2*(jup2*(-xCM + x[ijk]) - jup1*y[ijk])*pow2(Psi4)
;


} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* BNS_set_J_ADM_VolInt_integrand.c */
/* nvars = 16, n* = 82,  n/ = 31,  n+ = 70, n = 183, O = 1 */
