/* SingleBHKS.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 9.5.2005 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "SingleBH.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void SingleBHKS(tGrid *grid, int i_x, int i_gb, int i_K,                             int i_psi, int i_dpsiopsi, int i_ddpsiopsi,                            int i_alpha, int i_beta)
{
int bi;

double M = Getd("BHmass1");
double Sx = Getd("BHsx1");
double Sy = Getd("BHsy1");
double Sz = Getd("BHsz1");
int ConformalFactor = Getv("SingleBH_ConformalFactor","yes");

for(bi = 0; bi < grid->nboxes; bi++)
{
tBox *box = grid->box[bi];
int ijk;

forallpoints(box, ijk)
{
double *x1 = box->v[i_x+0];
double *x2 = box->v[i_x+1];
double *x3 = box->v[i_x+2];
double *gb11 = box->v[i_gb+0];
double *gb12 = box->v[i_gb+1];
double *gb13 = box->v[i_gb+2];
double *gb22 = box->v[i_gb+3];
double *gb23 = box->v[i_gb+4];
double *gb33 = box->v[i_gb+5];
double *K11 = box->v[i_K+0];
double *K12 = box->v[i_K+1];
double *K13 = box->v[i_K+2];
double *K22 = box->v[i_K+3];
double *K23 = box->v[i_K+4];
double *K33 = box->v[i_K+5];
double *psi = box->v[i_psi+0];
double *dpsiopsi1 = box->v[i_dpsiopsi+0];
double *dpsiopsi2 = box->v[i_dpsiopsi+1];
double *dpsiopsi3 = box->v[i_dpsiopsi+2];
double *ddpsiopsi11 = box->v[i_ddpsiopsi+0];
double *ddpsiopsi12 = box->v[i_ddpsiopsi+1];
double *ddpsiopsi13 = box->v[i_ddpsiopsi+2];
double *ddpsiopsi22 = box->v[i_ddpsiopsi+3];
double *ddpsiopsi23 = box->v[i_ddpsiopsi+4];
double *ddpsiopsi33 = box->v[i_ddpsiopsi+5];
double *alpha = box->v[i_alpha+0];
double *beta1 = box->v[i_beta+0];
double *beta2 = box->v[i_beta+1];
double *beta3 = box->v[i_beta+2];


double a2;
double adotx;
double detgup;
double H;
double Omega;
double r;
double rho;
double acrossx1;
double acrossx2;
double acrossx3;
double ahat12;
double ahat13;
double ahat23;
double aS1;
double aS2;
double aS3;
double dH1;
double dH2;
double dH3;
double dl11;
double dl12;
double dl13;
double dl21;
double dl22;
double dl23;
double dl31;
double dl32;
double dl33;
double dr1;
double dr2;
double dr3;
double g11;
double g12;
double g13;
double g22;
double g23;
double g33;
double gup11;
double gup12;
double gup13;
double gup22;
double gup23;
double gup33;
double l1;
double l2;
double l3;



/* Jetzt geht's los! */
rho
=
sqrt(pow2(x1[ijk]) + pow2(x2[ijk]) + pow2(x3[ijk]))
;

aS1
=
Sx
;

aS2
=
Sy
;

aS3
=
Sz
;

adotx
=
aS1*x1[ijk] + aS2*x2[ijk] + aS3*x3[ijk]
;

a2
=
pow2(aS1) + pow2(aS2) + pow2(aS3)
;

Omega
=
sqrt(Power(rho,4) + pow2(a2) + 4.*pow2(adotx) - 2.*a2*pow2(rho))
;

r
=
sqrt(0.5*(-a2 + Omega + pow2(rho)))
;

dr1
=
((adotx*aS1 - 0.5*a2*x1[ijk])/Omega + 
    0.5*(x1[ijk] + (x1[ijk]*pow2(rho))/Omega))/r
;

dr2
=
((adotx*aS2 - 0.5*a2*x2[ijk])/Omega + 
    0.5*(x2[ijk] + (x2[ijk]*pow2(rho))/Omega))/r
;

dr3
=
((adotx*aS3 - 0.5*a2*x3[ijk])/Omega + 
    0.5*(x3[ijk] + (x3[ijk]*pow2(rho))/Omega))/r
;

acrossx1
=
-(aS3*x2[ijk]) + aS2*x3[ijk]
;

acrossx2
=
aS3*x1[ijk] - aS1*x3[ijk]
;

acrossx3
=
-(aS2*x1[ijk]) + aS1*x2[ijk]
;

ahat12
=
-aS3
;

ahat13
=
aS2
;

ahat23
=
-aS1
;

H
=
(M*Power(r,3))/(Power(r,4) + pow2(adotx))
;

dH1
=
(H*(-2.*adotx*aS1*r + dr1*(-Power(r,4) + 3.*pow2(adotx))))/
  (Power(r,5) + r*pow2(adotx))
;

dH2
=
(H*(-2.*adotx*aS2*r + dr2*(-Power(r,4) + 3.*pow2(adotx))))/
  (Power(r,5) + r*pow2(adotx))
;

dH3
=
(H*(-2.*adotx*aS3*r + dr3*(-Power(r,4) + 3.*pow2(adotx))))/
  (Power(r,5) + r*pow2(adotx))
;

l1
=
(adotx*aS1 - acrossx1*r + x1[ijk]*pow2(r))/(a2*r + Power(r,3))
;

l2
=
(adotx*aS2 - acrossx2*r + x2[ijk]*pow2(r))/(a2*r + Power(r,3))
;

l3
=
(adotx*aS3 - acrossx3*r + x3[ijk]*pow2(r))/(a2*r + Power(r,3))
;

dl11
=
(pow2(aS1) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl12
=
(aS1*aS2 - ahat12*r)/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl13
=
(aS1*aS3 - ahat13*r)/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl21
=
(aS1*aS2 + ahat12*r)/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl22
=
(pow2(aS2) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl23
=
(aS2*aS3 - ahat23*r)/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl31
=
(aS1*aS3 + ahat13*r)/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl32
=
(aS2*aS3 + ahat23*r)/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl33
=
(pow2(aS3) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3[ijk]*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

g11
=
1. + 2.*H*pow2(l1)
;

g12
=
2.*H*l1*l2
;

g13
=
2.*H*l1*l3
;

g22
=
1. + 2.*H*pow2(l2)
;

g23
=
2.*H*l2*l3
;

g33
=
1. + 2.*H*pow2(l3)
;

detgup
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

gup11
=
detgup*(g22*g33 - pow2(g23))
;

gup12
=
detgup*(g13*g23 - g12*g33)
;

gup13
=
detgup*(-(g13*g22) + g12*g23)
;

gup22
=
detgup*(g11*g33 - pow2(g13))
;

gup23
=
detgup*(g12*g13 - g11*g23)
;

gup33
=
detgup*(g11*g22 - pow2(g12))
;

alpha[ijk]
=
1/sqrt(1. + 2.*H)
;

beta1[ijk]
=
(2.*H*l1)/(1. + 2.*H)
;

beta2[ijk]
=
(2.*H*l2)/(1. + 2.*H)
;

beta3[ijk]
=
(2.*H*l3)/(1. + 2.*H)
;

K11[ijk]
=
alpha[ijk]*(4.*pow2(H)*(l1*(dl12*l2 + dl13*l3) + dl11*pow2(l1)) + 
    2.*(dH1*l1 + H*(dl11 + dH1*Power(l1,3) + (dH2*l2 + dH3*l3)*pow2(l1))))
;

K12[ijk]
=
alpha[ijk]*((dl12 + dl21)*H + dH1*l2 + dH2*(l1 + 2.*H*l1*pow2(l2)) + 
    2.*(l2*(dl13*l3*pow2(H) + H*(dH3*l1*l3 + dH1*pow2(l1))) + 
       pow2(H)*(l1*((dl11 + dl22)*l2 + dl23*l3) + dl21*pow2(l1) + 
          dl12*pow2(l2))))
;

K13[ijk]
=
alpha[ijk]*(dH1*l3 + l1*(dH3 + 2.*dl32*l2*pow2(H)) + 
    H*(dl13 + dl31 + 2.*dH3*l1*pow2(l3)) + 
    2.*(l3*(dl12*l2*pow2(H) + H*(dH2*l1*l2 + dH1*pow2(l1))) + 
       pow2(H)*((dl11 + dl33)*l1*l3 + dl31*pow2(l1) + dl13*pow2(l3))))
;

K22[ijk]
=
alpha[ijk]*(4.*pow2(H)*(dl21*l1*l2 + dl22*pow2(l2)) + 
    l3*(4.*dl23*l2*pow2(H) + 2.*dH3*H*pow2(l2)) + 
    2.*(dH2*l2 + H*(dl22 + dH2*Power(l2,3) + dH1*l1*pow2(l2))))
;

K23[ijk]
=
alpha[ijk]*((dl23 + dl32)*H + dH2*l3 + l2*(dH3 + 2.*dl31*l1*pow2(H)) + 
    2.*(l3*((dl21*l1 + dl22*l2)*pow2(H) + l2*(dH1*H*l1 + dl33*pow2(H))) + 
       (dH2*H*l3 + dl32*pow2(H))*pow2(l2) + 
       (dH3*H*l2 + dl23*pow2(H))*pow2(l3)))
;

K33[ijk]
=
alpha[ijk]*(4.*(dl31*l1 + dl32*l2)*l3*pow2(H) + 
    (2.*dH1*H*l1 + 4.*dl33*pow2(H))*pow2(l3) + 
    2.*(dH3*l3 + H*(dl33 + dH3*Power(l3,3) + dH2*l2*pow2(l3))))
;



/* conditional */
if (ConformalFactor) {

psi[ijk]
=
1. + (2.*M)/r
;

dpsiopsi1[ijk]
=
(-2.*M*x1[ijk])/(Power(r,3)*psi[ijk])
;

dpsiopsi2[ijk]
=
(-2.*M*x2[ijk])/(Power(r,3)*psi[ijk])
;

dpsiopsi3[ijk]
=
(-2.*M*x3[ijk])/(Power(r,3)*psi[ijk])
;

ddpsiopsi11[ijk]
=
(M*(-2./Power(r,3) + (6.*pow2(x1[ijk]))/Power(r,5)))/psi[ijk]
;

ddpsiopsi12[ijk]
=
(6.*M*x1[ijk]*x2[ijk])/(Power(r,5)*psi[ijk])
;

ddpsiopsi13[ijk]
=
(6.*M*x1[ijk]*x3[ijk])/(Power(r,5)*psi[ijk])
;

ddpsiopsi22[ijk]
=
(M*(-2./Power(r,3) + (6.*pow2(x2[ijk]))/Power(r,5)))/psi[ijk]
;

ddpsiopsi23[ijk]
=
(6.*M*x2[ijk]*x3[ijk])/(Power(r,5)*psi[ijk])
;

ddpsiopsi33[ijk]
=
(M*(-2./Power(r,3) + (6.*pow2(x3[ijk]))/Power(r,5)))/psi[ijk]
;


} else { /* if (!ConformalFactor) */

psi[ijk]
=
1.
;

dpsiopsi1[ijk]
=
0
;

dpsiopsi2[ijk]
=
0
;

dpsiopsi3[ijk]
=
0
;

ddpsiopsi11[ijk]
=
0
;

ddpsiopsi12[ijk]
=
0
;

ddpsiopsi13[ijk]
=
0
;

ddpsiopsi22[ijk]
=
0
;

ddpsiopsi23[ijk]
=
0
;

ddpsiopsi33[ijk]
=
0
;

}
/* if (ConformalFactor) */


gb11[ijk]
=
g11/Power(psi[ijk],4)
;

gb12[ijk]
=
g12/Power(psi[ijk],4)
;

gb13[ijk]
=
g13/Power(psi[ijk],4)
;

gb22[ijk]
=
g22/Power(psi[ijk],4)
;

gb23[ijk]
=
g23/Power(psi[ijk],4)
;

gb33[ijk]
=
g33/Power(psi[ijk],4)
;

} /* end of points */
} /* end of boxes */


}  /* end of function */

/* SingleBHKS.c */
/* nvars = 29, n* = 429,  n/ = 82,  n+ = 290, n = 801, O = 1 */
