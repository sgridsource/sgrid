/* BNS_compute_new_q.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 25.3.2011 */
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




void BNS_compute_new_q(tGrid *grid, int iq)
{
int corot1 = Getv("BNSdata_rotationstate1","corotation");
int corot2 = Getv("BNSdata_rotationstate2","corotation");
double n = Getd("BNSdata_n");
double C1 = Getd("BNSdata_C1");
double C2 = Getd("BNSdata_C2");
double kappa = Getd("BNSdata_kappa");
double Omega = Getd("BNSdata_Omega");
double xCM = Getd("BNSdata_x_CM");

int bi;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;



double *q = box->v[iq+0];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
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
double bb;
double beta1;
double beta2;
double beta3;
double betadSigmaMinusCC;
double CC;
double DSigmaUp1;
double DSigmaUp2;
double DSigmaUp3;
double F;
double h;
double L2;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double Psi2;
double Psi4;
double Psim2;
double Psim4;
double Psim6;
double twoalpha2wdSigmapw;
double uzero;
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
double xi1;
double xi2;
double xi3;



/* Jetzt geht's los! */

FirstDerivsOf_S(box,index_BNSdata_Sigma,                   
                                   Ind("BNSdata_Sigmax"));


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
if (((bi <= 1 || bi == 5) && corot1) || (bi >= 2 && bi <= 4 && corot2)) {

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



/* conditional */
if (oouzerosqr == 0) {

oouzerosqr
=
-1.
;

}
/* if (oouzerosqr == 0) */




/* conditional */
if (oouzerosqr < 0) {

uzero
=
-1.
;


} else { /* if (!oouzerosqr < 0) */

uzero
=
Sqrt(1/oouzerosqr)
;

}
/* if (oouzerosqr < 0) */


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

F
=
uzero*(-alpha2 + Psi4*(vR1*xi1 + beta1*(vR1 + xi1) + vR2*xi2 + 
       beta2*(vR2 + xi2) + vR3*xi3 + beta3*(vR3 + xi3) + pow2(beta1) + 
       pow2(beta2) + pow2(beta3)))
;



/* conditional */
if (bi == 0 || bi == 1 || bi == 5) {

q[ijk]
=
(-1. + C1/F)/(1. + n)
;


} else { /* if (!bi == 0 || bi == 1 || bi == 5) */

q[ijk]
=
(-1. + C2/F)/(1. + n)
;

}
/* if (bi == 0 || bi == 1 || bi == 5) */



} else { /* if (!bi == 0 || bi == 1 || bi == 5) */



/* conditional */
if (bi <= 1 || bi == 5) {

CC
=
C1
;


} else { /* if (!bi <= 1 || bi == 5) */

CC
=
C2
;

}
/* if (bi <= 1 || bi == 5) */


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

twoalpha2wdSigmapw
=
2.*alpha2*(w1*(wDown1 + dSigma1[ijk]) + w2*(wDown2 + dSigma2[ijk]) + 
    w3*(wDown3 + dSigma3[ijk]))
;

betadSigmaMinusCC
=
-CC + beta1*dSigma1[ijk] + beta2*dSigma2[ijk] + beta3*dSigma3[ijk]
;

bb
=
twoalpha2wdSigmapw + pow2(betadSigmaMinusCC)
;

L2
=
(0.5*(bb + Sqrt(Abs(pow2(bb) + pow2(twoalpha2wdSigmapw)))))/alpha2
;

h
=
Sqrt(Abs(L2 - (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) - 
    (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) - 
    (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])))
;

q[ijk]
=
(-1. + h)/(1. + n)
;

}
/* if (bi <= 1 || bi == 5) */



} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* BNS_compute_new_q.c */
/* nvars = 15, n* = 107,  n/ = 60,  n+ = 92, n = 259, O = 1 */
