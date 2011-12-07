/* BNS_set_dlnIntegEuler.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 6.12.2011 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    (sqrt((double) (x)))
#define Abs(x)     (fabs((double) (x)))
#define Log(x)     (log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BNS_set_dlnIntegEuler(tGrid *grid, int ilnIntegEuler,                                  int idlnIntegEuler, double Om, double xcm)
{
int corot1 = Getv("BNSdata_rotationstate1","corotation");
int corot2 = Getv("BNSdata_rotationstate2","corotation");
double n = Getd("BNSdata_n");
double Omega = Getd("BNSdata_Omega");
double xCM = Getd("BNSdata_x_CM");

int bi;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;



double *lnIntegEuler = box->v[ilnIntegEuler+0];
double *dlnIntegEuler = box->v[idlnIntegEuler+0];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_BNSdata_q = Ind("BNSdata_q");
double *q = box->v[index_BNSdata_q + 0];
int index_BNSdata_Psi = Ind("BNSdata_Psi");
double *Psi = box->v[index_BNSdata_Psi + 0];
int index_BNSdata_alphaP = Ind("BNSdata_alphaP");
double *alphaP = box->v[index_BNSdata_alphaP + 0];
int index_BNSdata_Bx = Ind("BNSdata_Bx");
double *B1 = box->v[index_BNSdata_Bx + 0];
double *B2 = box->v[index_BNSdata_Bx + 1];
double *B3 = box->v[index_BNSdata_Bx + 2];
int index_BNSdata_wBx = Ind("BNSdata_wBx");
double *wB1 = box->v[index_BNSdata_wBx + 0];
double *wB2 = box->v[index_BNSdata_wBx + 1];
double *wB3 = box->v[index_BNSdata_wBx + 2];
int index_BNSdata_Sigma = Ind("BNSdata_Sigma");
double *Sigma = box->v[index_BNSdata_Sigma + 0];
int index_BNSdata_Sigmax = Ind("BNSdata_Sigmax");
double *dSigma1 = box->v[index_BNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_BNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_BNSdata_Sigmax + 2];


double alpha;
double alpha2;
double bet1;
double bet2;
double bet3;
double beta1;
double beta2;
double beta3;
double betxiw1;
double betxiw2;
double betxiw3;
double betxiwDown1;
double betxiwDown2;
double betxiwDown3;
double DSigmaUp1;
double DSigmaUp2;
double DSigmaUp3;
double Gamma;
double Gamma0;
double Gamman;
double h;
double h2;
double L2;
double OmCrossR1;
double OmCrossR2;
double OmCrossR3;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double Psi2;
double Psi4;
double Psim2;
double Psim4;
double Psim6;
double U01;
double U02;
double U03;
double U0Down1;
double U0Down2;
double U0Down3;
double U1;
double U2;
double U3;
double uzero;
double uzerosqr;
double w1;
double w2;
double w3;
double wBDown1;
double wBDown2;
double wBDown3;
double wDown1;
double wDown2;
double wDown3;
double xi1;
double xi2;
double xi3;



/* Jetzt geht's los! */

FirstDerivsOf_S(box,index_BNSdata_Sigma,                                    Ind("BNSdata_Sigmax")); 



/* conditional */
if (bi == 1) {


                                                                      
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmax"), 0,1);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmay"), 0,1);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmaz"), 0,1);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBx"), 0,1);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBy"), 0,1);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBz"), 0,1);}
/* if (bi == 1) */




/* conditional */
if (bi == 2) {


                                                                      
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmax"), 3,2);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmay"), 3,2);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_Sigmaz"), 3,2);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBx"), 3,2);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBy"), 3,2);
     copy_Var_at_i0_from_Box1_Box2(grid, Ind("BNSdata_wBz"), 3,2);}
/* if (bi == 2) */



forallpoints(box, ijk) { 

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

OmCrossR1
=
-(Om*y[ijk])
;

OmCrossR2
=
Om*(-xcm + x[ijk])
;

OmCrossR3
=
0
;

bet1
=
OmCrossR1 + B1[ijk]
;

bet2
=
OmCrossR2 + B2[ijk]
;

bet3
=
OmCrossR3 + B3[ijk]
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



/* conditional */
if ((bi <= 1 || bi == 5) && corot1 || bi >= 2 && bi <= 4 && corot2) {

oouzerosqr
=
alpha2 - Psi4*(pow2(bet1) + pow2(bet2) + pow2(bet3))
;



/* conditional */
if (oouzerosqr == 0) {

oouzerosqr
=
1.
;

}
/* if (oouzerosqr == 0) */


lnIntegEuler[ijk]
=
log(oouzerosqr)
;


} else { /* if (!oouzerosqr == 0) */

Psim4
=
1/Psi4
;

Psim2
=
Psi2*Psim4
;

Psim6
=
Psim2*Psim4
;

DSigmaUp1
=
Psim4*dSigma1[ijk]
;

DSigmaUp2
=
Psim4*dSigma2[ijk]
;

DSigmaUp3
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

L2
=
h2 + (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
  (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
  (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])
;

uzerosqr
=
L2/(alpha2*h2)
;

uzero
=
sqrt(uzerosqr)
;

xi1
=
0
;

xi2
=
0
;

xi3
=
0
;

U01
=
(beta1 + w1/(h*uzero) + xi1)/alpha
;

U02
=
(beta2 + w2/(h*uzero) + xi2)/alpha
;

U03
=
(beta3 + w3/(h*uzero) + xi3)/alpha
;

U0Down1
=
Psi4*U01
;

U0Down2
=
Psi4*U02
;

U0Down3
=
Psi4*U03
;

U1
=
DSigmaUp1/(alpha*h*uzero)
;

U2
=
DSigmaUp2/(alpha*h*uzero)
;

U3
=
DSigmaUp3/(alpha*h*uzero)
;

Gamman
=
alpha*uzero
;

Gamma0
=
1/sqrt(1. - U01*U0Down1 - U02*U0Down2 - U03*U0Down3)
;

Gamma
=
Gamma0*(Gamman - Gamman*(U0Down1*U1 + U0Down2*U2 + U0Down3*U3 + 
       (w1*wDown1 + w2*wDown2 + w3*wDown3)/(alpha2*h2*uzerosqr)))
;

betxiw1
=
bet1 + w1/(h*uzero) + xi1
;

betxiw2
=
bet2 + w2/(h*uzero) + xi2
;

betxiw3
=
bet3 + w3/(h*uzero) + xi3
;

betxiwDown1
=
betxiw1*Psi4
;

betxiwDown2
=
betxiw2*Psi4
;

betxiwDown3
=
betxiw3*Psi4
;

lnIntegEuler[ijk]
=
log(alpha2 - betxiw1*betxiwDown1 - betxiw2*betxiwDown2 - 
    betxiw3*betxiwDown3) + 2.*log(Gamma)
;

}
/* if (oouzerosqr == 0) */



} /* end of points loop */ 


FirstDerivsOf_S(box, ilnIntegEuler, idlnIntegEuler); 

} /* end of boxes */


}  /* end of function */

/* BNS_set_dlnIntegEuler.c */
/* nvars = 19, n* = 110,  n/ = 50,  n+ = 87, n = 247, O = 1 */
