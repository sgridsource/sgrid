/* Kerr3d.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 14.9.2007 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "ScalarOnKerr.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void Kerr3d(tGrid *grid, int i_x, int i_alpha, int i_beta, int i_g,    int i_K, int i_TrK, int i_gup, int i_Gam, int i_dalpha)
{
int bi;

double M = Getd("BHmass");
double Sx = Getd("BHsx");
double Sy = Getd("BHsy");
double Sz = Getd("BHsz");

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;

forallpoints(box, ijk)
{
double *x1 = box->v[i_x+0];
double *x2 = box->v[i_x+1];
double *x3 = box->v[i_x+2];
double *alpha = box->v[i_alpha+0];
double *beta1 = box->v[i_beta+0];
double *beta2 = box->v[i_beta+1];
double *beta3 = box->v[i_beta+2];
double *g11 = box->v[i_g+0];
double *g12 = box->v[i_g+1];
double *g13 = box->v[i_g+2];
double *g22 = box->v[i_g+3];
double *g23 = box->v[i_g+4];
double *g33 = box->v[i_g+5];
double *K11 = box->v[i_K+0];
double *K12 = box->v[i_K+1];
double *K13 = box->v[i_K+2];
double *K21 = box->v[i_K+3];
double *K22 = box->v[i_K+4];
double *K23 = box->v[i_K+5];
double *K31 = box->v[i_K+6];
double *K32 = box->v[i_K+7];
double *K33 = box->v[i_K+8];
double *TrK = box->v[i_TrK+0];
double *gup11 = box->v[i_gup+0];
double *gup12 = box->v[i_gup+1];
double *gup13 = box->v[i_gup+2];
double *gup22 = box->v[i_gup+3];
double *gup23 = box->v[i_gup+4];
double *gup33 = box->v[i_gup+5];
double *Gam111 = box->v[i_Gam+0];
double *Gam112 = box->v[i_Gam+1];
double *Gam113 = box->v[i_Gam+2];
double *Gam122 = box->v[i_Gam+3];
double *Gam123 = box->v[i_Gam+4];
double *Gam133 = box->v[i_Gam+5];
double *Gam211 = box->v[i_Gam+6];
double *Gam212 = box->v[i_Gam+7];
double *Gam213 = box->v[i_Gam+8];
double *Gam222 = box->v[i_Gam+9];
double *Gam223 = box->v[i_Gam+10];
double *Gam233 = box->v[i_Gam+11];
double *Gam311 = box->v[i_Gam+12];
double *Gam312 = box->v[i_Gam+13];
double *Gam313 = box->v[i_Gam+14];
double *Gam322 = box->v[i_Gam+15];
double *Gam323 = box->v[i_Gam+16];
double *Gam333 = box->v[i_Gam+17];
double *dalpha1 = box->v[i_dalpha+0];
double *dalpha2 = box->v[i_dalpha+1];
double *dalpha3 = box->v[i_dalpha+2];


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
double dg111;
double dg112;
double dg113;
double dg121;
double dg122;
double dg123;
double dg131;
double dg132;
double dg133;
double dg221;
double dg222;
double dg223;
double dg231;
double dg232;
double dg233;
double dg331;
double dg332;
double dg333;
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
double Gamdo111;
double Gamdo112;
double Gamdo113;
double Gamdo122;
double Gamdo123;
double Gamdo133;
double Gamdo211;
double Gamdo212;
double Gamdo213;
double Gamdo222;
double Gamdo223;
double Gamdo233;
double Gamdo311;
double Gamdo312;
double Gamdo313;
double Gamdo322;
double Gamdo323;
double Gamdo333;
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



/* conditional */
if (M != 0.) {

aS1
=
aS1/M
;

aS2
=
aS2/M
;

aS3
=
aS3/M
;

}
/* if (M != 0.) */


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

g11[ijk]
=
1. + 2.*H*pow2(l1)
;

g12[ijk]
=
2.*H*l1*l2
;

g13[ijk]
=
2.*H*l1*l3
;

g22[ijk]
=
1. + 2.*H*pow2(l2)
;

g23[ijk]
=
2.*H*l2*l3
;

g33[ijk]
=
1. + 2.*H*pow2(l3)
;

detgup
=
1/(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
    g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
    g11[ijk]*pow2(g23[ijk]))
;

gup11[ijk]
=
detgup*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

gup12[ijk]
=
detgup*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

gup13[ijk]
=
detgup*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

gup22[ijk]
=
detgup*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

gup23[ijk]
=
detgup*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

gup33[ijk]
=
detgup*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
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

dalpha1[ijk]
=
-(dH1/(sqrt(1. + 2.*H) + 2.*H*sqrt(1. + 2.*H)))
;

dalpha2[ijk]
=
-(dH2/(sqrt(1. + 2.*H) + 2.*H*sqrt(1. + 2.*H)))
;

dalpha3[ijk]
=
-(dH3/(sqrt(1. + 2.*H) + 2.*H*sqrt(1. + 2.*H)))
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

K21[ijk]
=
alpha[ijk]*((dl12 + dl21)*H + dH1*l2 + dH2*(l1 + 2.*H*l1*pow2(l2)) + 
    2.*(l2*(dl13*l3*pow2(H) + H*(dH3*l1*l3 + dH1*pow2(l1))) + 
       pow2(H)*(l1*((dl11 + dl22)*l2 + dl23*l3) + dl21*pow2(l1) + 
          dl12*pow2(l2))))
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

K31[ijk]
=
alpha[ijk]*(dH1*l3 + l1*(dH3 + 2.*dl32*l2*pow2(H)) + 
    H*(dl13 + dl31 + 2.*dH3*l1*pow2(l3)) + 
    2.*(l3*(dl12*l2*pow2(H) + H*(dH2*l1*l2 + dH1*pow2(l1))) + 
       pow2(H)*((dl11 + dl33)*l1*l3 + dl31*pow2(l1) + dl13*pow2(l3))))
;

K32[ijk]
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

TrK[ijk]
=
gup11[ijk]*K11[ijk] + gup12[ijk]*(K12[ijk] + K21[ijk]) + 
  gup22[ijk]*K22[ijk] + gup13[ijk]*(K13[ijk] + K31[ijk]) + 
  gup23[ijk]*(K23[ijk] + K32[ijk]) + gup33[ijk]*K33[ijk]
;

dg111
=
4.*dl11*H*l1 + 2.*dH1*pow2(l1)
;

dg112
=
4.*dl12*H*l1 + 2.*dH2*pow2(l1)
;

dg113
=
4.*dl13*H*l1 + 2.*dH3*pow2(l1)
;

dg121
=
2.*(dH1*l1*l2 + H*(dl21*l1 + dl11*l2))
;

dg122
=
2.*(dH2*l1*l2 + H*(dl22*l1 + dl12*l2))
;

dg123
=
2.*(dH3*l1*l2 + H*(dl23*l1 + dl13*l2))
;

dg131
=
2.*(dH1*l1*l3 + H*(dl31*l1 + dl11*l3))
;

dg132
=
2.*(dH2*l1*l3 + H*(dl32*l1 + dl12*l3))
;

dg133
=
2.*(dH3*l1*l3 + H*(dl33*l1 + dl13*l3))
;

dg221
=
4.*dl21*H*l2 + 2.*dH1*pow2(l2)
;

dg222
=
4.*dl22*H*l2 + 2.*dH2*pow2(l2)
;

dg223
=
4.*dl23*H*l2 + 2.*dH3*pow2(l2)
;

dg231
=
2.*(dH1*l2*l3 + H*(dl31*l2 + dl21*l3))
;

dg232
=
2.*(dH2*l2*l3 + H*(dl32*l2 + dl22*l3))
;

dg233
=
2.*(dH3*l2*l3 + H*(dl33*l2 + dl23*l3))
;

dg331
=
4.*dl31*H*l3 + 2.*dH1*pow2(l3)
;

dg332
=
4.*dl32*H*l3 + 2.*dH2*pow2(l3)
;

dg333
=
4.*dl33*H*l3 + 2.*dH3*pow2(l3)
;

Gamdo111
=
0.5*dg111
;

Gamdo112
=
0.5*dg112
;

Gamdo113
=
0.5*dg113
;

Gamdo122
=
dg122 - 0.5*dg221
;

Gamdo123
=
0.5*(dg123 + dg132 - dg231)
;

Gamdo133
=
dg133 - 0.5*dg331
;

Gamdo211
=
-0.5*dg112 + dg121
;

Gamdo212
=
0.5*dg221
;

Gamdo213
=
0.5*(dg123 - dg132 + dg231)
;

Gamdo222
=
0.5*dg222
;

Gamdo223
=
0.5*dg223
;

Gamdo233
=
dg233 - 0.5*dg332
;

Gamdo311
=
-0.5*dg113 + dg131
;

Gamdo312
=
0.5*(-dg123 + dg132 + dg231)
;

Gamdo313
=
0.5*dg331
;

Gamdo322
=
-0.5*dg223 + dg232
;

Gamdo323
=
0.5*dg332
;

Gamdo333
=
0.5*dg333
;

Gam111[ijk]
=
Gamdo111*gup11[ijk] + Gamdo211*gup12[ijk] + Gamdo311*gup13[ijk]
;

Gam112[ijk]
=
Gamdo112*gup11[ijk] + Gamdo212*gup12[ijk] + Gamdo312*gup13[ijk]
;

Gam113[ijk]
=
Gamdo113*gup11[ijk] + Gamdo213*gup12[ijk] + Gamdo313*gup13[ijk]
;

Gam122[ijk]
=
Gamdo122*gup11[ijk] + Gamdo222*gup12[ijk] + Gamdo322*gup13[ijk]
;

Gam123[ijk]
=
Gamdo123*gup11[ijk] + Gamdo223*gup12[ijk] + Gamdo323*gup13[ijk]
;

Gam133[ijk]
=
Gamdo133*gup11[ijk] + Gamdo233*gup12[ijk] + Gamdo333*gup13[ijk]
;

Gam211[ijk]
=
Gamdo111*gup12[ijk] + Gamdo211*gup22[ijk] + Gamdo311*gup23[ijk]
;

Gam212[ijk]
=
Gamdo112*gup12[ijk] + Gamdo212*gup22[ijk] + Gamdo312*gup23[ijk]
;

Gam213[ijk]
=
Gamdo113*gup12[ijk] + Gamdo213*gup22[ijk] + Gamdo313*gup23[ijk]
;

Gam222[ijk]
=
Gamdo122*gup12[ijk] + Gamdo222*gup22[ijk] + Gamdo322*gup23[ijk]
;

Gam223[ijk]
=
Gamdo123*gup12[ijk] + Gamdo223*gup22[ijk] + Gamdo323*gup23[ijk]
;

Gam233[ijk]
=
Gamdo133*gup12[ijk] + Gamdo233*gup22[ijk] + Gamdo333*gup23[ijk]
;

Gam311[ijk]
=
Gamdo111*gup13[ijk] + Gamdo211*gup23[ijk] + Gamdo311*gup33[ijk]
;

Gam312[ijk]
=
Gamdo112*gup13[ijk] + Gamdo212*gup23[ijk] + Gamdo312*gup33[ijk]
;

Gam313[ijk]
=
Gamdo113*gup13[ijk] + Gamdo213*gup23[ijk] + Gamdo313*gup33[ijk]
;

Gam322[ijk]
=
Gamdo122*gup13[ijk] + Gamdo222*gup23[ijk] + Gamdo322*gup33[ijk]
;

Gam323[ijk]
=
Gamdo123*gup13[ijk] + Gamdo223*gup23[ijk] + Gamdo323*gup33[ijk]
;

Gam333[ijk]
=
Gamdo133*gup13[ijk] + Gamdo233*gup23[ijk] + Gamdo333*gup33[ijk]
;

} /* end of points */
} /* end of boxes */


}  /* end of function */

/* Kerr3d.c */
/* nvars = 50, n* = 678,  n/ = 64,  n+ = 454, n = 1196, O = 1 */
