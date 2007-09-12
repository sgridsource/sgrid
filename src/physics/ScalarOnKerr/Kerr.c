/* Kerr.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 12.9.2007 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "ScalarOnKerr.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void Kerr(tGrid *grid, int i_x, int i_g, int i_gup, int i_Gam)
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
double *xp1 = box->v[i_x+0];
double *xp2 = box->v[i_x+1];
double *xp3 = box->v[i_x+2];
double *g00 = box->v[i_g+0];
double *g01 = box->v[i_g+1];
double *g02 = box->v[i_g+2];
double *g03 = box->v[i_g+3];
double *g11 = box->v[i_g+4];
double *g12 = box->v[i_g+5];
double *g13 = box->v[i_g+6];
double *g22 = box->v[i_g+7];
double *g23 = box->v[i_g+8];
double *g33 = box->v[i_g+9];
double *gup00 = box->v[i_gup+0];
double *gup01 = box->v[i_gup+1];
double *gup02 = box->v[i_gup+2];
double *gup03 = box->v[i_gup+3];
double *gup11 = box->v[i_gup+4];
double *gup12 = box->v[i_gup+5];
double *gup13 = box->v[i_gup+6];
double *gup22 = box->v[i_gup+7];
double *gup23 = box->v[i_gup+8];
double *gup33 = box->v[i_gup+9];
double *Gam000 = box->v[i_Gam+0];
double *Gam001 = box->v[i_Gam+1];
double *Gam002 = box->v[i_Gam+2];
double *Gam003 = box->v[i_Gam+3];
double *Gam011 = box->v[i_Gam+4];
double *Gam012 = box->v[i_Gam+5];
double *Gam013 = box->v[i_Gam+6];
double *Gam022 = box->v[i_Gam+7];
double *Gam023 = box->v[i_Gam+8];
double *Gam033 = box->v[i_Gam+9];
double *Gam100 = box->v[i_Gam+10];
double *Gam101 = box->v[i_Gam+11];
double *Gam102 = box->v[i_Gam+12];
double *Gam103 = box->v[i_Gam+13];
double *Gam111 = box->v[i_Gam+14];
double *Gam112 = box->v[i_Gam+15];
double *Gam113 = box->v[i_Gam+16];
double *Gam122 = box->v[i_Gam+17];
double *Gam123 = box->v[i_Gam+18];
double *Gam133 = box->v[i_Gam+19];
double *Gam200 = box->v[i_Gam+20];
double *Gam201 = box->v[i_Gam+21];
double *Gam202 = box->v[i_Gam+22];
double *Gam203 = box->v[i_Gam+23];
double *Gam211 = box->v[i_Gam+24];
double *Gam212 = box->v[i_Gam+25];
double *Gam213 = box->v[i_Gam+26];
double *Gam222 = box->v[i_Gam+27];
double *Gam223 = box->v[i_Gam+28];
double *Gam233 = box->v[i_Gam+29];
double *Gam300 = box->v[i_Gam+30];
double *Gam301 = box->v[i_Gam+31];
double *Gam302 = box->v[i_Gam+32];
double *Gam303 = box->v[i_Gam+33];
double *Gam311 = box->v[i_Gam+34];
double *Gam312 = box->v[i_Gam+35];
double *Gam313 = box->v[i_Gam+36];
double *Gam322 = box->v[i_Gam+37];
double *Gam323 = box->v[i_Gam+38];
double *Gam333 = box->v[i_Gam+39];


double a2;
double acrossx0;
double acrossx1;
double acrossx2;
double acrossx3;
double adotx;
double ahat01;
double ahat02;
double ahat03;
double ahat12;
double ahat13;
double ahat23;
double aS0;
double aS1;
double aS2;
double aS3;
double dg000;
double dg001;
double dg002;
double dg003;
double dg010;
double dg011;
double dg012;
double dg013;
double dg020;
double dg021;
double dg022;
double dg023;
double dg030;
double dg031;
double dg032;
double dg033;
double dg110;
double dg111;
double dg112;
double dg113;
double dg120;
double dg121;
double dg122;
double dg123;
double dg130;
double dg131;
double dg132;
double dg133;
double dg220;
double dg221;
double dg222;
double dg223;
double dg230;
double dg231;
double dg232;
double dg233;
double dg330;
double dg331;
double dg332;
double dg333;
double dH0;
double dH1;
double dH2;
double dH3;
double dl00;
double dl01;
double dl02;
double dl03;
double dl10;
double dl11;
double dl12;
double dl13;
double dl20;
double dl21;
double dl22;
double dl23;
double dl30;
double dl31;
double dl32;
double dl33;
double dr0;
double dr1;
double dr2;
double dr3;
double Gamdo000;
double Gamdo001;
double Gamdo002;
double Gamdo003;
double Gamdo011;
double Gamdo012;
double Gamdo013;
double Gamdo022;
double Gamdo023;
double Gamdo033;
double Gamdo100;
double Gamdo101;
double Gamdo102;
double Gamdo103;
double Gamdo111;
double Gamdo112;
double Gamdo113;
double Gamdo122;
double Gamdo123;
double Gamdo133;
double Gamdo200;
double Gamdo201;
double Gamdo202;
double Gamdo203;
double Gamdo211;
double Gamdo212;
double Gamdo213;
double Gamdo222;
double Gamdo223;
double Gamdo233;
double Gamdo300;
double Gamdo301;
double Gamdo302;
double Gamdo303;
double Gamdo311;
double Gamdo312;
double Gamdo313;
double Gamdo322;
double Gamdo323;
double Gamdo333;
double H;
double l0;
double l1;
double l2;
double l3;
double Omega;
double r;
double rho;
double x0;
double x1;
double x2;
double x3;



/* Jetzt geht's los! */
x0
=
0
;

x1
=
xp1[ijk]
;

x2
=
xp2[ijk]
;

x3
=
xp3[ijk]
;

rho
=
sqrt(pow2(x1) + pow2(x2) + pow2(x3))
;

aS0
=
0
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
if (M != 0) {

aS0
=
aS0/M
;

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
/* if (M != 0) */


adotx
=
aS1*x1 + aS2*x2 + aS3*x3
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

dr0
=
((adotx*aS0 - 0.5*a2*x0)/Omega + 0.5*(x0 + (x0*pow2(rho))/Omega))/r
;

dr1
=
((adotx*aS1 - 0.5*a2*x1)/Omega + 0.5*(x1 + (x1*pow2(rho))/Omega))/r
;

dr2
=
((adotx*aS2 - 0.5*a2*x2)/Omega + 0.5*(x2 + (x2*pow2(rho))/Omega))/r
;

dr3
=
((adotx*aS3 - 0.5*a2*x3)/Omega + 0.5*(x3 + (x3*pow2(rho))/Omega))/r
;

dr0
=
0
;

acrossx0
=
0
;

acrossx1
=
-(aS3*x2) + aS2*x3
;

acrossx2
=
aS3*x1 - aS1*x3
;

acrossx3
=
-(aS2*x1) + aS1*x2
;

ahat01
=
0
;

ahat02
=
0
;

ahat03
=
0
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

dH0
=
(H*(-2.*adotx*aS0*r + dr0*(-Power(r,4) + 3.*pow2(adotx))))/
  (Power(r,5) + r*pow2(adotx))
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

dH0
=
0
;

l0
=
(adotx*aS0 - acrossx0*r + x0*pow2(r))/(a2*r + Power(r,3))
;

l1
=
(adotx*aS1 - acrossx1*r + x1*pow2(r))/(a2*r + Power(r,3))
;

l2
=
(adotx*aS2 - acrossx2*r + x2*pow2(r))/(a2*r + Power(r,3))
;

l3
=
(adotx*aS3 - acrossx3*r + x3*pow2(r))/(a2*r + Power(r,3))
;

l0
=
-1.
;

dl00
=
(pow2(aS0) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr0*(2.*acrossx0*Power(r,3) - adotx*aS0*(a2 + 3.*pow2(r)) + 
       x0*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl01
=
(aS0*aS1 - ahat01*r)/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx0*Power(r,3) - adotx*aS0*(a2 + 3.*pow2(r)) + 
       x0*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl02
=
(aS0*aS2 - ahat02*r)/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx0*Power(r,3) - adotx*aS0*(a2 + 3.*pow2(r)) + 
       x0*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl03
=
(aS0*aS3 - ahat03*r)/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx0*Power(r,3) - adotx*aS0*(a2 + 3.*pow2(r)) + 
       x0*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl10
=
(aS0*aS1 + ahat01*r)/(a2*r + Power(r,3)) + 
  (dr0*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl11
=
(pow2(aS1) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl12
=
(aS1*aS2 - ahat12*r)/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl13
=
(aS1*aS3 - ahat13*r)/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx1*Power(r,3) - adotx*aS1*(a2 + 3.*pow2(r)) + 
       x1*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl20
=
(aS0*aS2 + ahat02*r)/(a2*r + Power(r,3)) + 
  (dr0*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl21
=
(aS1*aS2 + ahat12*r)/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl22
=
(pow2(aS2) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl23
=
(aS2*aS3 - ahat23*r)/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx2*Power(r,3) - adotx*aS2*(a2 + 3.*pow2(r)) + 
       x2*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl30
=
(aS0*aS3 + ahat03*r)/(a2*r + Power(r,3)) + 
  (dr0*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl31
=
(aS1*aS3 + ahat13*r)/(a2*r + Power(r,3)) + 
  (dr1*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl32
=
(aS2*aS3 + ahat23*r)/(a2*r + Power(r,3)) + 
  (dr2*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl33
=
(pow2(aS3) + pow2(r))/(a2*r + Power(r,3)) + 
  (dr3*(2.*acrossx3*Power(r,3) - adotx*aS3*(a2 + 3.*pow2(r)) + 
       x3*(-Power(r,4) + a2*pow2(r))))/
   (2.*a2*Power(r,4) + Power(r,6) + pow2(a2)*pow2(r))
;

dl00
=
0
;

dl10
=
0
;

dl20
=
0
;

dl30
=
0
;

g00[ijk]
=
-1. + 2.*H*pow2(l0)
;

g01[ijk]
=
2.*H*l0*l1
;

g02[ijk]
=
2.*H*l0*l2
;

g03[ijk]
=
2.*H*l0*l3
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

gup00[ijk]
=
-1. - 2.*H*pow2(l0)
;

gup01[ijk]
=
-2.*H*l0*l1
;

gup02[ijk]
=
-2.*H*l0*l2
;

gup03[ijk]
=
-2.*H*l0*l3
;

gup11[ijk]
=
1. - 2.*H*pow2(l1)
;

gup12[ijk]
=
-2.*H*l1*l2
;

gup13[ijk]
=
-2.*H*l1*l3
;

gup22[ijk]
=
1. - 2.*H*pow2(l2)
;

gup23[ijk]
=
-2.*H*l2*l3
;

gup33[ijk]
=
1. - 2.*H*pow2(l3)
;

dg000
=
4.*dl00*H*l0 + 2.*dH0*pow2(l0)
;

dg001
=
4.*dl01*H*l0 + 2.*dH1*pow2(l0)
;

dg002
=
4.*dl02*H*l0 + 2.*dH2*pow2(l0)
;

dg003
=
4.*dl03*H*l0 + 2.*dH3*pow2(l0)
;

dg010
=
2.*(dH0*l0*l1 + H*(dl10*l0 + dl00*l1))
;

dg011
=
2.*(dH1*l0*l1 + H*(dl11*l0 + dl01*l1))
;

dg012
=
2.*(dH2*l0*l1 + H*(dl12*l0 + dl02*l1))
;

dg013
=
2.*(dH3*l0*l1 + H*(dl13*l0 + dl03*l1))
;

dg020
=
2.*(dH0*l0*l2 + H*(dl20*l0 + dl00*l2))
;

dg021
=
2.*(dH1*l0*l2 + H*(dl21*l0 + dl01*l2))
;

dg022
=
2.*(dH2*l0*l2 + H*(dl22*l0 + dl02*l2))
;

dg023
=
2.*(dH3*l0*l2 + H*(dl23*l0 + dl03*l2))
;

dg030
=
2.*(dH0*l0*l3 + H*(dl30*l0 + dl00*l3))
;

dg031
=
2.*(dH1*l0*l3 + H*(dl31*l0 + dl01*l3))
;

dg032
=
2.*(dH2*l0*l3 + H*(dl32*l0 + dl02*l3))
;

dg033
=
2.*(dH3*l0*l3 + H*(dl33*l0 + dl03*l3))
;

dg110
=
4.*dl10*H*l1 + 2.*dH0*pow2(l1)
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

dg120
=
2.*(dH0*l1*l2 + H*(dl20*l1 + dl10*l2))
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

dg130
=
2.*(dH0*l1*l3 + H*(dl30*l1 + dl10*l3))
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

dg220
=
4.*dl20*H*l2 + 2.*dH0*pow2(l2)
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

dg230
=
2.*(dH0*l2*l3 + H*(dl30*l2 + dl20*l3))
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

dg330
=
4.*dl30*H*l3 + 2.*dH0*pow2(l3)
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

Gamdo000
=
0.5*dg000
;

Gamdo001
=
0.5*dg001
;

Gamdo002
=
0.5*dg002
;

Gamdo003
=
0.5*dg003
;

Gamdo011
=
dg011 - 0.5*dg110
;

Gamdo012
=
0.5*(dg012 + dg021 - dg120)
;

Gamdo013
=
0.5*(dg013 + dg031 - dg130)
;

Gamdo022
=
dg022 - 0.5*dg220
;

Gamdo023
=
0.5*(dg023 + dg032 - dg230)
;

Gamdo033
=
dg033 - 0.5*dg330
;

Gamdo100
=
-0.5*dg001 + dg010
;

Gamdo101
=
0.5*dg110
;

Gamdo102
=
0.5*(dg012 - dg021 + dg120)
;

Gamdo103
=
0.5*(dg013 - dg031 + dg130)
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

Gamdo200
=
-0.5*dg002 + dg020
;

Gamdo201
=
0.5*(-dg012 + dg021 + dg120)
;

Gamdo202
=
0.5*dg220
;

Gamdo203
=
0.5*(dg023 - dg032 + dg230)
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

Gamdo300
=
-0.5*dg003 + dg030
;

Gamdo301
=
0.5*(-dg013 + dg031 + dg130)
;

Gamdo302
=
0.5*(-dg023 + dg032 + dg230)
;

Gamdo303
=
0.5*dg330
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

Gam000[ijk]
=
Gamdo000*gup00[ijk] + Gamdo100*gup01[ijk] + Gamdo200*gup02[ijk] + 
  Gamdo300*gup03[ijk]
;

Gam001[ijk]
=
Gamdo001*gup00[ijk] + Gamdo101*gup01[ijk] + Gamdo201*gup02[ijk] + 
  Gamdo301*gup03[ijk]
;

Gam002[ijk]
=
Gamdo002*gup00[ijk] + Gamdo102*gup01[ijk] + Gamdo202*gup02[ijk] + 
  Gamdo302*gup03[ijk]
;

Gam003[ijk]
=
Gamdo003*gup00[ijk] + Gamdo103*gup01[ijk] + Gamdo203*gup02[ijk] + 
  Gamdo303*gup03[ijk]
;

Gam011[ijk]
=
Gamdo011*gup00[ijk] + Gamdo111*gup01[ijk] + Gamdo211*gup02[ijk] + 
  Gamdo311*gup03[ijk]
;

Gam012[ijk]
=
Gamdo012*gup00[ijk] + Gamdo112*gup01[ijk] + Gamdo212*gup02[ijk] + 
  Gamdo312*gup03[ijk]
;

Gam013[ijk]
=
Gamdo013*gup00[ijk] + Gamdo113*gup01[ijk] + Gamdo213*gup02[ijk] + 
  Gamdo313*gup03[ijk]
;

Gam022[ijk]
=
Gamdo022*gup00[ijk] + Gamdo122*gup01[ijk] + Gamdo222*gup02[ijk] + 
  Gamdo322*gup03[ijk]
;

Gam023[ijk]
=
Gamdo023*gup00[ijk] + Gamdo123*gup01[ijk] + Gamdo223*gup02[ijk] + 
  Gamdo323*gup03[ijk]
;

Gam033[ijk]
=
Gamdo033*gup00[ijk] + Gamdo133*gup01[ijk] + Gamdo233*gup02[ijk] + 
  Gamdo333*gup03[ijk]
;

Gam100[ijk]
=
Gamdo000*gup01[ijk] + Gamdo100*gup11[ijk] + Gamdo200*gup12[ijk] + 
  Gamdo300*gup13[ijk]
;

Gam101[ijk]
=
Gamdo001*gup01[ijk] + Gamdo101*gup11[ijk] + Gamdo201*gup12[ijk] + 
  Gamdo301*gup13[ijk]
;

Gam102[ijk]
=
Gamdo002*gup01[ijk] + Gamdo102*gup11[ijk] + Gamdo202*gup12[ijk] + 
  Gamdo302*gup13[ijk]
;

Gam103[ijk]
=
Gamdo003*gup01[ijk] + Gamdo103*gup11[ijk] + Gamdo203*gup12[ijk] + 
  Gamdo303*gup13[ijk]
;

Gam111[ijk]
=
Gamdo011*gup01[ijk] + Gamdo111*gup11[ijk] + Gamdo211*gup12[ijk] + 
  Gamdo311*gup13[ijk]
;

Gam112[ijk]
=
Gamdo012*gup01[ijk] + Gamdo112*gup11[ijk] + Gamdo212*gup12[ijk] + 
  Gamdo312*gup13[ijk]
;

Gam113[ijk]
=
Gamdo013*gup01[ijk] + Gamdo113*gup11[ijk] + Gamdo213*gup12[ijk] + 
  Gamdo313*gup13[ijk]
;

Gam122[ijk]
=
Gamdo022*gup01[ijk] + Gamdo122*gup11[ijk] + Gamdo222*gup12[ijk] + 
  Gamdo322*gup13[ijk]
;

Gam123[ijk]
=
Gamdo023*gup01[ijk] + Gamdo123*gup11[ijk] + Gamdo223*gup12[ijk] + 
  Gamdo323*gup13[ijk]
;

Gam133[ijk]
=
Gamdo033*gup01[ijk] + Gamdo133*gup11[ijk] + Gamdo233*gup12[ijk] + 
  Gamdo333*gup13[ijk]
;

Gam200[ijk]
=
Gamdo000*gup02[ijk] + Gamdo100*gup12[ijk] + Gamdo200*gup22[ijk] + 
  Gamdo300*gup23[ijk]
;

Gam201[ijk]
=
Gamdo001*gup02[ijk] + Gamdo101*gup12[ijk] + Gamdo201*gup22[ijk] + 
  Gamdo301*gup23[ijk]
;

Gam202[ijk]
=
Gamdo002*gup02[ijk] + Gamdo102*gup12[ijk] + Gamdo202*gup22[ijk] + 
  Gamdo302*gup23[ijk]
;

Gam203[ijk]
=
Gamdo003*gup02[ijk] + Gamdo103*gup12[ijk] + Gamdo203*gup22[ijk] + 
  Gamdo303*gup23[ijk]
;

Gam211[ijk]
=
Gamdo011*gup02[ijk] + Gamdo111*gup12[ijk] + Gamdo211*gup22[ijk] + 
  Gamdo311*gup23[ijk]
;

Gam212[ijk]
=
Gamdo012*gup02[ijk] + Gamdo112*gup12[ijk] + Gamdo212*gup22[ijk] + 
  Gamdo312*gup23[ijk]
;

Gam213[ijk]
=
Gamdo013*gup02[ijk] + Gamdo113*gup12[ijk] + Gamdo213*gup22[ijk] + 
  Gamdo313*gup23[ijk]
;

Gam222[ijk]
=
Gamdo022*gup02[ijk] + Gamdo122*gup12[ijk] + Gamdo222*gup22[ijk] + 
  Gamdo322*gup23[ijk]
;

Gam223[ijk]
=
Gamdo023*gup02[ijk] + Gamdo123*gup12[ijk] + Gamdo223*gup22[ijk] + 
  Gamdo323*gup23[ijk]
;

Gam233[ijk]
=
Gamdo033*gup02[ijk] + Gamdo133*gup12[ijk] + Gamdo233*gup22[ijk] + 
  Gamdo333*gup23[ijk]
;

Gam300[ijk]
=
Gamdo000*gup03[ijk] + Gamdo100*gup13[ijk] + Gamdo200*gup23[ijk] + 
  Gamdo300*gup33[ijk]
;

Gam301[ijk]
=
Gamdo001*gup03[ijk] + Gamdo101*gup13[ijk] + Gamdo201*gup23[ijk] + 
  Gamdo301*gup33[ijk]
;

Gam302[ijk]
=
Gamdo002*gup03[ijk] + Gamdo102*gup13[ijk] + Gamdo202*gup23[ijk] + 
  Gamdo302*gup33[ijk]
;

Gam303[ijk]
=
Gamdo003*gup03[ijk] + Gamdo103*gup13[ijk] + Gamdo203*gup23[ijk] + 
  Gamdo303*gup33[ijk]
;

Gam311[ijk]
=
Gamdo011*gup03[ijk] + Gamdo111*gup13[ijk] + Gamdo211*gup23[ijk] + 
  Gamdo311*gup33[ijk]
;

Gam312[ijk]
=
Gamdo012*gup03[ijk] + Gamdo112*gup13[ijk] + Gamdo212*gup23[ijk] + 
  Gamdo312*gup33[ijk]
;

Gam313[ijk]
=
Gamdo013*gup03[ijk] + Gamdo113*gup13[ijk] + Gamdo213*gup23[ijk] + 
  Gamdo313*gup33[ijk]
;

Gam322[ijk]
=
Gamdo022*gup03[ijk] + Gamdo122*gup13[ijk] + Gamdo222*gup23[ijk] + 
  Gamdo322*gup33[ijk]
;

Gam323[ijk]
=
Gamdo023*gup03[ijk] + Gamdo123*gup13[ijk] + Gamdo223*gup23[ijk] + 
  Gamdo323*gup33[ijk]
;

Gam333[ijk]
=
Gamdo033*gup03[ijk] + Gamdo133*gup13[ijk] + Gamdo233*gup23[ijk] + 
  Gamdo333*gup33[ijk]
;

} /* end of points */
} /* end of boxes */


}  /* end of function */

/* Kerr.c */
/* nvars = 64, n* = 855,  n/ = 76,  n+ = 598, n = 1529, O = 1 */
