/* BNS_compute_new_q_atXYZ.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 27.11.2010 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    (sqrt((double) (x)))
#define Log(x)     (log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




double BNS_compute_new_q_atXYZ(tGrid *grid, int bi, double X, double Y, double Z)
{
double n = Getd("BNSdata_n");
double C1 = Getd("BNSdata_C1");
double C2 = Getd("BNSdata_C2");
double kappa = Getd("BNSdata_kappa");
double Omega = Getd("BNSdata_Omega");
double xCM = Getd("BNSdata_x_CM");

tBox *box = grid->box[bi];
int ijk;



int index_BNSdata_Psi = Ind("BNSdata_Psi");
double *BNSPsi = box->v[index_BNSdata_Psi + 0];
int index_BNSdata_alphaP = Ind("BNSdata_alphaP");
double *BNSalphaP = box->v[index_BNSdata_alphaP + 0];
int index_BNSdata_Bx = Ind("BNSdata_Bx");
double *BNSB1 = box->v[index_BNSdata_Bx + 0];
double *BNSB2 = box->v[index_BNSdata_Bx + 1];
double *BNSB3 = box->v[index_BNSdata_Bx + 2];
int index_BNSdata_q = Ind("BNSdata_q");
double *BNSq = box->v[index_BNSdata_q + 0];
int index_BNSdata_wBx = Ind("BNSdata_wBx");
double *BNSwB1 = box->v[index_BNSdata_wBx + 0];
double *BNSwB2 = box->v[index_BNSdata_wBx + 1];
double *BNSwB3 = box->v[index_BNSdata_wBx + 2];
int index_BNSdata_Sigma = Ind("BNSdata_Sigma");
double *BNSSigma = box->v[index_BNSdata_Sigma + 0];
int index_BNSdata_Sigmax = Ind("BNSdata_Sigmax");
double *BNSdSigma1 = box->v[index_BNSdata_Sigmax + 0];
double *BNSdSigma2 = box->v[index_BNSdata_Sigmax + 1];
double *BNSdSigma3 = box->v[index_BNSdata_Sigmax + 2];
int index_BNSdata_temp4 = Ind("BNSdata_temp4");
double *temp4 = box->v[index_BNSdata_temp4 + 0];
double Psi, B1,B2,B3, alphaP;
double Sigma, dSigma1,dSigma2,dSigma3, wB1,wB2,wB3, x,y;


double alpha;
double alpha2;
double beta1;
double beta2;
double beta3;
double F;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double Psi2;
double Psi4;
double q;
double uzero;
double vR1;
double vR2;
double vR3;
double vRI1;
double vRI2;
double vRI3;
double xi1;
double xi2;
double xi3;



/* Jetzt geht's los! */

FirstDerivsOf_S(box,index_BNSdata_Sigma,                                    Ind("BNSdata_Sigmax")); 


if(box->x_of_X[1] != NULL) { 


x = box->x_of_X[1]((void *) box, -1, X,Y,Z); 


y = box->x_of_X[2]((void *) box, -1, X,Y,Z); 


} else { 


x = X; 


y = Y; 


} 


spec_Coeffs(box, BNSPsi, temp4); 


Psi = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSB1, temp4); 


B1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSB2, temp4); 


B2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSB3, temp4); 


B3 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSalphaP, temp4); 


alphaP = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSSigma, temp4); 


Sigma = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSdSigma1, temp4); 


dSigma1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSdSigma2, temp4); 


dSigma2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSdSigma3, temp4); 


dSigma3 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSwB1, temp4); 


wB1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSwB2, temp4); 


wB2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, BNSwB3, temp4); 


wB3 = spec_interpolate(box, temp4, X,Y,Z); 

OmegaCrossR1
=
-(Omega*y)
;

OmegaCrossR2
=
Omega*(x - xCM)
;

OmegaCrossR3
=
0
;

beta1
=
B1 + OmegaCrossR1
;

beta2
=
B2 + OmegaCrossR2
;

beta3
=
B3 + OmegaCrossR3
;

alpha
=
alphaP/Psi
;

alpha2
=
pow2(alpha)
;

Psi2
=
pow2(Psi)
;

Psi4
=
pow2(Psi2)
;

vRI1
=
dSigma1
;

vRI2
=
dSigma2
;

vRI3
=
dSigma3
;

vR1
=
vRI1 + wB1
;

vR2
=
vRI2 + wB2
;

vR3
=
vRI3 + wB3
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

q
=
(-1. + C1/F)/(1. + n)
;


} else { /* if (!bi == 0 || bi == 1 || bi == 5) */

q
=
(-1. + C2/F)/(1. + n)
;

}
/* if (bi == 0 || bi == 1 || bi == 5) */



/* end of computation */ 

return q;
}  /* end of function */

/* BNS_compute_new_q_atXYZ.c */
/* nvars = 13, n* = 63,  n/ = 35,  n+ = 73, n = 171, O = 1 */
