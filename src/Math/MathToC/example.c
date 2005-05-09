/* example.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 9.5.2005 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "RealisticBBH.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void example(tBox *box, tVarList *vlMetric, tVarList *vlCurv)
{
int ijk=0;

double *psi = box->v[Ind("psi")+0];
double *g11 = box->v[Ind("gxx")+0];
double *g12 = box->v[Ind("gxx")+1];
double *g13 = box->v[Ind("gxx")+2];
double *g22 = box->v[Ind("gxx")+3];
double *g23 = box->v[Ind("gxx")+4];
double *g33 = box->v[Ind("gxx")+5];
double *K11 = box->v[Ind("Kxx")+0];
double *K12 = box->v[Ind("Kxx")+1];
double *K13 = box->v[Ind("Kxx")+2];
double *K22 = box->v[Ind("Kxx")+3];
double *K23 = box->v[Ind("Kxx")+4];
double *K33 = box->v[Ind("Kxx")+5];
double *AR[1, 1] = box->v[Ind("RealisticBBH_ARxx")+0];
double *AR[1, 2] = box->v[Ind("RealisticBBH_ARxx")+1];
double *AR[1, 3] = box->v[Ind("RealisticBBH_ARxx")+2];
double *AR[2, 1] = box->v[Ind("RealisticBBH_ARxx")+3];
double *AR[2, 2] = box->v[Ind("RealisticBBH_ARxx")+4];
double *AR[2, 3] = box->v[Ind("RealisticBBH_ARxx")+5];
double *AR[3, 1] = box->v[Ind("RealisticBBH_ARxx")+6];
double *AR[3, 2] = box->v[Ind("RealisticBBH_ARxx")+7];
double *AR[3, 3] = box->v[Ind("RealisticBBH_ARxx")+8];
double *ATT[1, 1] = box->v[Ind("RealisticBBH_ATTxx")+0];
double *ATT[1, 2] = box->v[Ind("RealisticBBH_ATTxx")+1];
double *ATT[1, 3] = box->v[Ind("RealisticBBH_ATTxx")+2];
double *ATT[2, 1] = box->v[Ind("RealisticBBH_ATTxx")+3];
double *ATT[2, 2] = box->v[Ind("RealisticBBH_ATTxx")+4];
double *ATT[2, 3] = box->v[Ind("RealisticBBH_ATTxx")+5];
double *ATT[3, 1] = box->v[Ind("RealisticBBH_ATTxx")+6];
double *ATT[3, 2] = box->v[Ind("RealisticBBH_ATTxx")+7];
double *ATT[3, 3] = box->v[Ind("RealisticBBH_ATTxx")+8];
double *tau = box->v[Ind("RealisticBBH_tau")+0];
double *U = box->v[Ind("RealisticBBH_U")+0];
double *LW11 = box->v[Ind("RealisticBBH_LWxx")+0];
double *LW12 = box->v[Ind("RealisticBBH_LWxx")+1];
double *LW13 = box->v[Ind("RealisticBBH_LWxx")+2];
double *LW22 = box->v[Ind("RealisticBBH_LWxx")+3];
double *LW23 = box->v[Ind("RealisticBBH_LWxx")+4];
double *LW33 = box->v[Ind("RealisticBBH_LWxx")+5];


double Psi;
double Psi2;
double Psi4;
double psqrPN1;
double Bdown11;
double Bdown12;
double Bdown13;
double Bdown22;
double Bdown23;
double Bdown33;
double Bup11;
double Bup12;
double Bup13;
double Bup22;
double Bup23;
double Bup33;
double V111;
double V112;
double V113;
double V121;
double V122;
double V123;
double V131;
double V132;
double V133;
double V211;
double V212;
double V213;
double V221;
double V222;
double V223;
double V231;
double V232;
double V233;
double V311;
double V312;
double V313;
double V321;
double V322;
double V323;
double V331;
double V332;
double V333;
double W111;
double W112;
double W113;
double W121;
double W122;
double W123;
double W131;
double W132;
double W133;
double W211;
double W212;
double W213;
double W221;
double W222;
double W223;
double W231;
double W232;
double W233;
double W311;
double W312;
double W313;
double W321;
double W322;
double W323;
double W331;
double W332;
double W333;
double Z1111;
double Z1112;
double Z1113;
double Z1121;
double Z1122;
double Z1123;
double Z1131;
double Z1132;
double Z1133;
double Z1211;
double Z1212;
double Z1213;
double Z1221;
double Z1222;
double Z1223;
double Z1231;
double Z1232;
double Z1233;
double Z1311;
double Z1312;
double Z1313;
double Z1321;
double Z1322;
double Z1323;
double Z1331;
double Z1332;
double Z1333;
double Z2111;
double Z2112;
double Z2113;
double Z2121;
double Z2122;
double Z2123;
double Z2131;
double Z2132;
double Z2133;
double Z2211;
double Z2212;
double Z2213;
double Z2221;
double Z2222;
double Z2223;
double Z2231;
double Z2232;
double Z2233;
double Z2311;
double Z2312;
double Z2313;
double Z2321;
double Z2322;
double Z2323;
double Z2331;
double Z2332;
double Z2333;
double Z3111;
double Z3112;
double Z3113;
double Z3121;
double Z3122;
double Z3123;
double Z3131;
double Z3132;
double Z3133;
double Z3211;
double Z3212;
double Z3213;
double Z3221;
double Z3222;
double Z3223;
double Z3231;
double Z3232;
double Z3233;
double Z3311;
double Z3312;
double Z3313;
double Z3321;
double Z3322;
double Z3323;
double Z3331;
double Z3332;
double Z3333;



/* Jetzt geht's los! */

forallpoints(box,ijk) { 

Psi
=
U + psi[ijk]
;

Psi2
=
pow2(Psi)
;

Psi4
=
pow2(Psi2)
;

g11[ijk]
=
(Power(Psi,4)*g11[ijk])/Power(psi[ijk],4)
;

g12[ijk]
=
(Power(Psi,4)*g12[ijk])/Power(psi[ijk],4)
;

g13[ijk]
=
(Power(Psi,4)*g13[ijk])/Power(psi[ijk],4)
;

g22[ijk]
=
(Power(Psi,4)*g22[ijk])/Power(psi[ijk],4)
;

g23[ijk]
=
(Power(Psi,4)*g23[ijk])/Power(psi[ijk],4)
;

g33[ijk]
=
(Power(Psi,4)*g33[ijk])/Power(psi[ijk],4)
;

psqrPN1
=
pow2(mu)*(1/r12m + 4.*pow2inv(r12m))
;

W111
=
dg111[ijk]
;

W112
=
dg121[ijk]
;

W113
=
dg131[ijk]
;

W121
=
dg121[ijk]
;

W122
=
dg221[ijk]
;

W123
=
dg231[ijk]
;

W131
=
dg131[ijk]
;

W132
=
dg231[ijk]
;

W133
=
dg331[ijk]
;

W211
=
dg112[ijk]
;

W212
=
dg122[ijk]
;

W213
=
dg132[ijk]
;

W221
=
dg122[ijk]
;

W222
=
dg222[ijk]
;

W223
=
dg232[ijk]
;

W231
=
dg132[ijk]
;

W232
=
dg232[ijk]
;

W233
=
dg332[ijk]
;

W311
=
dg113[ijk]
;

W312
=
dg123[ijk]
;

W313
=
dg133[ijk]
;

W321
=
dg123[ijk]
;

W322
=
dg223[ijk]
;

W323
=
dg233[ijk]
;

W331
=
dg133[ijk]
;

W332
=
dg233[ijk]
;

W333
=
dg333[ijk]
;

Z1111
=
ddg1111[ijk]
;

Z1112
=
ddg1211[ijk]
;

Z1113
=
ddg1311[ijk]
;

Z1121
=
ddg1211[ijk]
;

Z1122
=
ddg2211[ijk]
;

Z1123
=
ddg2311[ijk]
;

Z1131
=
ddg1311[ijk]
;

Z1132
=
ddg2311[ijk]
;

Z1133
=
ddg3311[ijk]
;

Z1211
=
ddg1112[ijk]
;

Z1212
=
ddg1212[ijk]
;

Z1213
=
ddg1312[ijk]
;

Z1221
=
ddg1212[ijk]
;

Z1222
=
ddg2212[ijk]
;

Z1223
=
ddg2312[ijk]
;

Z1231
=
ddg1312[ijk]
;

Z1232
=
ddg2312[ijk]
;

Z1233
=
ddg3312[ijk]
;

Z1311
=
ddg1113[ijk]
;

Z1312
=
ddg1213[ijk]
;

Z1313
=
ddg1313[ijk]
;

Z1321
=
ddg1213[ijk]
;

Z1322
=
ddg2213[ijk]
;

Z1323
=
ddg2313[ijk]
;

Z1331
=
ddg1313[ijk]
;

Z1332
=
ddg2313[ijk]
;

Z1333
=
ddg3313[ijk]
;

Z2111
=
ddg1112[ijk]
;

Z2112
=
ddg1212[ijk]
;

Z2113
=
ddg1312[ijk]
;

Z2121
=
ddg1212[ijk]
;

Z2122
=
ddg2212[ijk]
;

Z2123
=
ddg2312[ijk]
;

Z2131
=
ddg1312[ijk]
;

Z2132
=
ddg2312[ijk]
;

Z2133
=
ddg3312[ijk]
;

Z2211
=
ddg1122[ijk]
;

Z2212
=
ddg1222[ijk]
;

Z2213
=
ddg1322[ijk]
;

Z2221
=
ddg1222[ijk]
;

Z2222
=
ddg2222[ijk]
;

Z2223
=
ddg2322[ijk]
;

Z2231
=
ddg1322[ijk]
;

Z2232
=
ddg2322[ijk]
;

Z2233
=
ddg3322[ijk]
;

Z2311
=
ddg1123[ijk]
;

Z2312
=
ddg1223[ijk]
;

Z2313
=
ddg1323[ijk]
;

Z2321
=
ddg1223[ijk]
;

Z2322
=
ddg2223[ijk]
;

Z2323
=
ddg2323[ijk]
;

Z2331
=
ddg1323[ijk]
;

Z2332
=
ddg2323[ijk]
;

Z2333
=
ddg3323[ijk]
;

Z3111
=
ddg1113[ijk]
;

Z3112
=
ddg1213[ijk]
;

Z3113
=
ddg1313[ijk]
;

Z3121
=
ddg1213[ijk]
;

Z3122
=
ddg2213[ijk]
;

Z3123
=
ddg2313[ijk]
;

Z3131
=
ddg1313[ijk]
;

Z3132
=
ddg2313[ijk]
;

Z3133
=
ddg3313[ijk]
;

Z3211
=
ddg1123[ijk]
;

Z3212
=
ddg1223[ijk]
;

Z3213
=
ddg1323[ijk]
;

Z3221
=
ddg1223[ijk]
;

Z3222
=
ddg2223[ijk]
;

Z3223
=
ddg2323[ijk]
;

Z3231
=
ddg1323[ijk]
;

Z3232
=
ddg2323[ijk]
;

Z3233
=
ddg3323[ijk]
;

Z3311
=
ddg1133[ijk]
;

Z3312
=
ddg1233[ijk]
;

Z3313
=
ddg1333[ijk]
;

Z3321
=
ddg1233[ijk]
;

Z3322
=
ddg2233[ijk]
;

Z3323
=
ddg2333[ijk]
;

Z3331
=
ddg1333[ijk]
;

Z3332
=
ddg2333[ijk]
;

Z3333
=
ddg3333[ijk]
;

W111
=
d1g11[ijk]
;

W112
=
d1g12[ijk]
;

W113
=
d1g13[ijk]
;

W121
=
d1g12[ijk]
;

W122
=
d1g22[ijk]
;

W123
=
d1g23[ijk]
;

W131
=
d1g13[ijk]
;

W132
=
d1g23[ijk]
;

W133
=
d1g33[ijk]
;

W211
=
d2g11[ijk]
;

W212
=
d2g12[ijk]
;

W213
=
d2g13[ijk]
;

W221
=
d2g12[ijk]
;

W222
=
d2g22[ijk]
;

W223
=
d2g23[ijk]
;

W231
=
d2g13[ijk]
;

W232
=
d2g23[ijk]
;

W233
=
d2g33[ijk]
;

W311
=
d3g11[ijk]
;

W312
=
d3g12[ijk]
;

W313
=
d3g13[ijk]
;

W321
=
d3g12[ijk]
;

W322
=
d3g22[ijk]
;

W323
=
d3g23[ijk]
;

W331
=
d3g13[ijk]
;

W332
=
d3g23[ijk]
;

W333
=
d3g33[ijk]
;

Z1111
=
d1d1g11[ijk]
;

Z1112
=
d1d1g12[ijk]
;

Z1113
=
d1d1g13[ijk]
;

Z1121
=
d1d1g12[ijk]
;

Z1122
=
d1d1g22[ijk]
;

Z1123
=
d1d1g23[ijk]
;

Z1131
=
d1d1g13[ijk]
;

Z1132
=
d1d1g23[ijk]
;

Z1133
=
d1d1g33[ijk]
;

Z1211
=
d1d2g11[ijk]
;

Z1212
=
d1d2g12[ijk]
;

Z1213
=
d1d2g13[ijk]
;

Z1221
=
d1d2g12[ijk]
;

Z1222
=
d1d2g22[ijk]
;

Z1223
=
d1d2g23[ijk]
;

Z1231
=
d1d2g13[ijk]
;

Z1232
=
d1d2g23[ijk]
;

Z1233
=
d1d2g33[ijk]
;

Z1311
=
d1d3g11[ijk]
;

Z1312
=
d1d3g12[ijk]
;

Z1313
=
d1d3g13[ijk]
;

Z1321
=
d1d3g12[ijk]
;

Z1322
=
d1d3g22[ijk]
;

Z1323
=
d1d3g23[ijk]
;

Z1331
=
d1d3g13[ijk]
;

Z1332
=
d1d3g23[ijk]
;

Z1333
=
d1d3g33[ijk]
;

Z2111
=
d1d2g11[ijk]
;

Z2112
=
d1d2g12[ijk]
;

Z2113
=
d1d2g13[ijk]
;

Z2121
=
d1d2g12[ijk]
;

Z2122
=
d1d2g22[ijk]
;

Z2123
=
d1d2g23[ijk]
;

Z2131
=
d1d2g13[ijk]
;

Z2132
=
d1d2g23[ijk]
;

Z2133
=
d1d2g33[ijk]
;

Z2211
=
d2d2g11[ijk]
;

Z2212
=
d2d2g12[ijk]
;

Z2213
=
d2d2g13[ijk]
;

Z2221
=
d2d2g12[ijk]
;

Z2222
=
d2d2g22[ijk]
;

Z2223
=
d2d2g23[ijk]
;

Z2231
=
d2d2g13[ijk]
;

Z2232
=
d2d2g23[ijk]
;

Z2233
=
d2d2g33[ijk]
;

Z2311
=
d2d3g11[ijk]
;

Z2312
=
d2d3g12[ijk]
;

Z2313
=
d2d3g13[ijk]
;

Z2321
=
d2d3g12[ijk]
;

Z2322
=
d2d3g22[ijk]
;

Z2323
=
d2d3g23[ijk]
;

Z2331
=
d2d3g13[ijk]
;

Z2332
=
d2d3g23[ijk]
;

Z2333
=
d2d3g33[ijk]
;

Z3111
=
d1d3g11[ijk]
;

Z3112
=
d1d3g12[ijk]
;

Z3113
=
d1d3g13[ijk]
;

Z3121
=
d1d3g12[ijk]
;

Z3122
=
d1d3g22[ijk]
;

Z3123
=
d1d3g23[ijk]
;

Z3131
=
d1d3g13[ijk]
;

Z3132
=
d1d3g23[ijk]
;

Z3133
=
d1d3g33[ijk]
;

Z3211
=
d2d3g11[ijk]
;

Z3212
=
d2d3g12[ijk]
;

Z3213
=
d2d3g13[ijk]
;

Z3221
=
d2d3g12[ijk]
;

Z3222
=
d2d3g22[ijk]
;

Z3223
=
d2d3g23[ijk]
;

Z3231
=
d2d3g13[ijk]
;

Z3232
=
d2d3g23[ijk]
;

Z3233
=
d2d3g33[ijk]
;

Z3311
=
d3d3g11[ijk]
;

Z3312
=
d3d3g12[ijk]
;

Z3313
=
d3d3g13[ijk]
;

Z3321
=
d3d3g12[ijk]
;

Z3322
=
d3d3g22[ijk]
;

Z3323
=
d3d3g23[ijk]
;

Z3331
=
d3d3g13[ijk]
;

Z3332
=
d3d3g23[ijk]
;

Z3333
=
d3d3g33[ijk]
;

V111
=
func(1.,g11[ijk])
;

V112
=
func(1.,g12[ijk])
;

V113
=
func(1.,g13[ijk])
;

V121
=
func(1.,g12[ijk])
;

V122
=
func(1.,g22[ijk])
;

V123
=
func(1.,g23[ijk])
;

V131
=
func(1.,g13[ijk])
;

V132
=
func(1.,g23[ijk])
;

V133
=
func(1.,g33[ijk])
;

V211
=
func(2.,g11[ijk])
;

V212
=
func(2.,g12[ijk])
;

V213
=
func(2.,g13[ijk])
;

V221
=
func(2.,g12[ijk])
;

V222
=
func(2.,g22[ijk])
;

V223
=
func(2.,g23[ijk])
;

V231
=
func(2.,g13[ijk])
;

V232
=
func(2.,g23[ijk])
;

V233
=
func(2.,g33[ijk])
;

V311
=
func(3.,g11[ijk])
;

V312
=
func(3.,g12[ijk])
;

V313
=
func(3.,g13[ijk])
;

V321
=
func(3.,g12[ijk])
;

V322
=
func(3.,g22[ijk])
;

V323
=
func(3.,g23[ijk])
;

V331
=
func(3.,g13[ijk])
;

V332
=
func(3.,g23[ijk])
;

V333
=
func(3.,g33[ijk])
;

Bup11
=
LW11[ijk]
;

Bup12
=
LW12[ijk]
;

Bup13
=
LW13[ijk]
;

Bup22
=
LW22[ijk]
;

Bup23
=
LW23[ijk]
;

Bup33
=
LW33[ijk]
;

Bdown11
=
2.*(Bup23*g12[ijk]*g13[ijk] + g11[ijk]*(Bup12*g12[ijk] + Bup13*g13[ijk])) + 
  Bup11*pow2(g11[ijk]) + Bup22*pow2(g12[ijk]) + Bup33*pow2(g13[ijk])
;

Bdown12
=
g12[ijk]*(Bup11*g11[ijk] + Bup13*g13[ijk] + Bup22*g22[ijk]) + 
  (Bup13*g11[ijk] + Bup33*g13[ijk])*g23[ijk] + 
  Bup23*(g13[ijk]*g22[ijk] + g12[ijk]*g23[ijk]) + 
  Bup12*(g11[ijk]*g22[ijk] + pow2(g12[ijk]))
;

Bdown13
=
(Bup12*g11[ijk] + Bup22*g12[ijk])*g23[ijk] + Bup23*g12[ijk]*g33[ijk] + 
  g13[ijk]*(Bup11*g11[ijk] + Bup12*g12[ijk] + Bup23*g23[ijk] + 
     Bup33*g33[ijk]) + Bup13*(g11[ijk]*g33[ijk] + pow2(g13[ijk]))
;

Bdown22
=
2.*(Bup23*g22[ijk]*g23[ijk] + g12[ijk]*(Bup12*g22[ijk] + Bup13*g23[ijk])) + 
  Bup11*pow2(g12[ijk]) + Bup22*pow2(g22[ijk]) + Bup33*pow2(g23[ijk])
;

Bdown23
=
g13[ijk]*(Bup11*g12[ijk] + Bup12*g22[ijk] + Bup13*g23[ijk]) + 
  Bup13*g12[ijk]*g33[ijk] + g23[ijk]*
   (Bup12*g12[ijk] + Bup22*g22[ijk] + Bup33*g33[ijk]) + 
  Bup23*(g22[ijk]*g33[ijk] + pow2(g23[ijk]))
;

Bdown33
=
2.*(Bup23*g23[ijk]*g33[ijk] + g13[ijk]*(Bup12*g23[ijk] + Bup13*g33[ijk])) + 
  Bup11*pow2(g13[ijk]) + Bup22*pow2(g23[ijk]) + Bup33*pow2(g33[ijk])
;

K11[ijk]
=
Bdown11
;

K12[ijk]
=
Bdown12
;

K13[ijk]
=
Bdown13
;

K22[ijk]
=
Bdown22
;

K23[ijk]
=
Bdown23
;

K33[ijk]
=
Bdown33
;


}  /* end loop */ 



}  /* end of function */

/* example.c */
/* nvars = 19, n* = 127,  n/ = 20,  n+ = 119, n = 266, O = 1 */
