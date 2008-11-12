/* BNS_set_restmassintegrand.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 11.11.2008 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Abs(x)     (fabs((double) (x)))
#define Sqrt(x)    (sqrt((double) (x)))
#define Log(x)     (log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BNS_set_restmassintegrand(tGrid *grid, int iInteg)
{
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
int index_BNSdata_q = Ind("BNSdata_q");
double *q = box->v[index_BNSdata_q + 0];
int index_BNSdata_vRSx = Ind("BNSdata_vRSx");
double *vRS1 = box->v[index_BNSdata_vRSx + 0];
double *vRS2 = box->v[index_BNSdata_vRSx + 1];
double *vRS3 = box->v[index_BNSdata_vRSx + 2];
int index_BNSdata_Sigma = Ind("BNSdata_Sigma");
double *Sigma = box->v[index_BNSdata_Sigma + 0];
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
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double Psi2;
double Psi4;
double rho0;
double uzero;
double vR1;
double vR2;
double vR3;
double vRI1;
double vRI2;
double vRI3;



/* Jetzt geht's los! */

FirstDerivsOf_S(box,index_BNSdata_Sigma,                                    Ind("BNSdata_Sigmax")); 


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

vRI1
=
dSigma1[ijk]
;

vRI2
=
dSigma2[ijk]
;

vRI3
=
dSigma3[ijk]
;

vR1
=
vRI1 + vRS1[ijk]
;

vR2
=
vRI2 + vRS2[ijk]
;

vR3
=
vRI3 + vRS3[ijk]
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




/* conditional */
if (q[ijk] >= 0.) {

rho0
=
Power(q[ijk]/kappa,n)
;


} else { /* if (!q[ijk] >= 0.) */

rho0
=
-Power(Abs(q[ijk]/kappa),n)
;

}
/* if (q[ijk] >= 0.) */


Integ[ijk]
=
alpha*Psi2*Psi4*rho0*uzero
;


} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* BNS_set_restmassintegrand.c */
/* nvars = 16, n* = 61,  n/ = 35,  n+ = 54, n = 150, O = 1 */
