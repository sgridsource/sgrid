/* ADMconstraints.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 18.12.2007 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "ADMvars.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void ADMconstraints(tVarList *u)
{
tGrid *grid = u->grid;
int bi;

int normConstr = Getv("ADMvars_normalizedConstraints", "yes");
int TermByTerm = GetvLax("ADMvars_ConstraintNorm", "TermByTerm");

for(bi = 0; bi < grid->nboxes; bi++)
{
tBox *box = grid->box[bi];
int ijk;

allDerivsOf_Sab(box, Ind("gxx"), Ind("ADMvars_dgxxx"),                       Ind("ADMvars_ddgxxxx"));
FirstDerivsOf_Sab(box, Ind("Kxx"), Ind("ADMvars_dKxxx"));

forallpoints(box, ijk)
{
double *psi = vlldataptr(psiandderivs, box, 0);
double *dpop1 = vlldataptr(psiandderivs, box, 1);
double *dpop2 = vlldataptr(psiandderivs, box, 2);
double *dpop3 = vlldataptr(psiandderivs, box, 3);
double *ddpop11 = vlldataptr(psiandderivs, box, 4);
double *ddpop12 = vlldataptr(psiandderivs, box, 5);
double *ddpop13 = vlldataptr(psiandderivs, box, 6);
double *ddpop22 = vlldataptr(psiandderivs, box, 7);
double *ddpop23 = vlldataptr(psiandderivs, box, 8);
double *ddpop33 = vlldataptr(psiandderivs, box, 9);
double *g11 = vlldataptr(u, box, 0);
double *g12 = vlldataptr(u, box, 1);
double *g13 = vlldataptr(u, box, 2);
double *g22 = vlldataptr(u, box, 3);
double *g23 = vlldataptr(u, box, 4);
double *g33 = vlldataptr(u, box, 5);
double *K11 = vlldataptr(u, box, 6);
double *K12 = vlldataptr(u, box, 7);
double *K13 = vlldataptr(u, box, 8);
double *K22 = vlldataptr(u, box, 9);
double *K23 = vlldataptr(u, box, 10);
double *K33 = vlldataptr(u, box, 11);
double *rho = vlldataptr(u, box, 12);
double *j1 = vlldataptr(u, box, 13);
double *j2 = vlldataptr(u, box, 14);
double *j3 = vlldataptr(u, box, 15);
double *ham = vlldataptr(u, box, 16);
double *mom1 = vlldataptr(u, box, 17);
double *mom2 = vlldataptr(u, box, 18);
double *mom3 = vlldataptr(u, box, 19);
double *trK = vlldataptr(u, box, 20);
double *normham = vlldataptr(u, box, 21);
double *normmom1 = vlldataptr(u, box, 22);
double *normmom2 = vlldataptr(u, box, 23);
double *normmom3 = vlldataptr(u, box, 24);
int index_ADMvars_dgxxx = Ind("ADMvars_dgxxx");
double *dg111 = box->v[index_ADMvars_dgxxx + 0];
double *dg112 = box->v[index_ADMvars_dgxxx + 1];
double *dg113 = box->v[index_ADMvars_dgxxx + 2];
double *dg121 = box->v[index_ADMvars_dgxxx + 3];
double *dg122 = box->v[index_ADMvars_dgxxx + 4];
double *dg123 = box->v[index_ADMvars_dgxxx + 5];
double *dg131 = box->v[index_ADMvars_dgxxx + 6];
double *dg132 = box->v[index_ADMvars_dgxxx + 7];
double *dg133 = box->v[index_ADMvars_dgxxx + 8];
double *dg221 = box->v[index_ADMvars_dgxxx + 9];
double *dg222 = box->v[index_ADMvars_dgxxx + 10];
double *dg223 = box->v[index_ADMvars_dgxxx + 11];
double *dg231 = box->v[index_ADMvars_dgxxx + 12];
double *dg232 = box->v[index_ADMvars_dgxxx + 13];
double *dg233 = box->v[index_ADMvars_dgxxx + 14];
double *dg331 = box->v[index_ADMvars_dgxxx + 15];
double *dg332 = box->v[index_ADMvars_dgxxx + 16];
double *dg333 = box->v[index_ADMvars_dgxxx + 17];
int index_ADMvars_ddgxxxx = Ind("ADMvars_ddgxxxx");
double *ddg1111 = box->v[index_ADMvars_ddgxxxx + 0];
double *ddg1112 = box->v[index_ADMvars_ddgxxxx + 1];
double *ddg1113 = box->v[index_ADMvars_ddgxxxx + 2];
double *ddg1122 = box->v[index_ADMvars_ddgxxxx + 3];
double *ddg1123 = box->v[index_ADMvars_ddgxxxx + 4];
double *ddg1133 = box->v[index_ADMvars_ddgxxxx + 5];
double *ddg1211 = box->v[index_ADMvars_ddgxxxx + 6];
double *ddg1212 = box->v[index_ADMvars_ddgxxxx + 7];
double *ddg1213 = box->v[index_ADMvars_ddgxxxx + 8];
double *ddg1222 = box->v[index_ADMvars_ddgxxxx + 9];
double *ddg1223 = box->v[index_ADMvars_ddgxxxx + 10];
double *ddg1233 = box->v[index_ADMvars_ddgxxxx + 11];
double *ddg1311 = box->v[index_ADMvars_ddgxxxx + 12];
double *ddg1312 = box->v[index_ADMvars_ddgxxxx + 13];
double *ddg1313 = box->v[index_ADMvars_ddgxxxx + 14];
double *ddg1322 = box->v[index_ADMvars_ddgxxxx + 15];
double *ddg1323 = box->v[index_ADMvars_ddgxxxx + 16];
double *ddg1333 = box->v[index_ADMvars_ddgxxxx + 17];
double *ddg2211 = box->v[index_ADMvars_ddgxxxx + 18];
double *ddg2212 = box->v[index_ADMvars_ddgxxxx + 19];
double *ddg2213 = box->v[index_ADMvars_ddgxxxx + 20];
double *ddg2222 = box->v[index_ADMvars_ddgxxxx + 21];
double *ddg2223 = box->v[index_ADMvars_ddgxxxx + 22];
double *ddg2233 = box->v[index_ADMvars_ddgxxxx + 23];
double *ddg2311 = box->v[index_ADMvars_ddgxxxx + 24];
double *ddg2312 = box->v[index_ADMvars_ddgxxxx + 25];
double *ddg2313 = box->v[index_ADMvars_ddgxxxx + 26];
double *ddg2322 = box->v[index_ADMvars_ddgxxxx + 27];
double *ddg2323 = box->v[index_ADMvars_ddgxxxx + 28];
double *ddg2333 = box->v[index_ADMvars_ddgxxxx + 29];
double *ddg3311 = box->v[index_ADMvars_ddgxxxx + 30];
double *ddg3312 = box->v[index_ADMvars_ddgxxxx + 31];
double *ddg3313 = box->v[index_ADMvars_ddgxxxx + 32];
double *ddg3322 = box->v[index_ADMvars_ddgxxxx + 33];
double *ddg3323 = box->v[index_ADMvars_ddgxxxx + 34];
double *ddg3333 = box->v[index_ADMvars_ddgxxxx + 35];
int index_ADMvars_dKxxx = Ind("ADMvars_dKxxx");
double *dK111 = box->v[index_ADMvars_dKxxx + 0];
double *dK112 = box->v[index_ADMvars_dKxxx + 1];
double *dK113 = box->v[index_ADMvars_dKxxx + 2];
double *dK121 = box->v[index_ADMvars_dKxxx + 3];
double *dK122 = box->v[index_ADMvars_dKxxx + 4];
double *dK123 = box->v[index_ADMvars_dKxxx + 5];
double *dK131 = box->v[index_ADMvars_dKxxx + 6];
double *dK132 = box->v[index_ADMvars_dKxxx + 7];
double *dK133 = box->v[index_ADMvars_dKxxx + 8];
double *dK221 = box->v[index_ADMvars_dKxxx + 9];
double *dK222 = box->v[index_ADMvars_dKxxx + 10];
double *dK223 = box->v[index_ADMvars_dKxxx + 11];
double *dK231 = box->v[index_ADMvars_dKxxx + 12];
double *dK232 = box->v[index_ADMvars_dKxxx + 13];
double *dK233 = box->v[index_ADMvars_dKxxx + 14];
double *dK331 = box->v[index_ADMvars_dKxxx + 15];
double *dK332 = box->v[index_ADMvars_dKxxx + 16];
double *dK333 = box->v[index_ADMvars_dKxxx + 17];

double denom;
double detginvf;
double f;
double hamrhs;
double K;
double KudKud;
double R;
double RA;
double RB;
double RC;
double RD;
double cdKdA1;
double cdKdA2;
double cdKdA3;
double cdKdB1;
double cdKdB2;
double cdKdB3;
double cdKdC1;
double cdKdC2;
double cdKdC3;
double cdKudd111;
double cdKudd112;
double cdKudd113;
double cdKudd122;
double cdKudd123;
double cdKudd133;
double cdKudd211;
double cdKudd212;
double cdKudd213;
double cdKudd222;
double cdKudd223;
double cdKudd233;
double cdKudd311;
double cdKudd312;
double cdKudd313;
double cdKudd322;
double cdKudd323;
double cdKudd333;
double codelK111;
double codelK112;
double codelK113;
double codelK122;
double codelK123;
double codelK133;
double codelK211;
double codelK212;
double codelK213;
double codelK222;
double codelK223;
double codelK233;
double codelK311;
double codelK312;
double codelK313;
double codelK322;
double codelK323;
double codelK333;
double codelKA111;
double codelKA112;
double codelKA113;
double codelKA122;
double codelKA123;
double codelKA133;
double codelKA211;
double codelKA212;
double codelKA213;
double codelKA222;
double codelKA223;
double codelKA233;
double codelKA311;
double codelKA312;
double codelKA313;
double codelKA322;
double codelKA323;
double codelKA333;
double codelKB111;
double codelKB112;
double codelKB113;
double codelKB122;
double codelKB123;
double codelKB133;
double codelKB211;
double codelKB212;
double codelKB213;
double codelKB222;
double codelKB223;
double codelKB233;
double codelKB311;
double codelKB312;
double codelKB313;
double codelKB322;
double codelKB323;
double codelKB333;
double codelKC111;
double codelKC112;
double codelKC113;
double codelKC122;
double codelKC123;
double codelKC133;
double codelKC211;
double codelKC212;
double codelKC213;
double codelKC222;
double codelKC223;
double codelKC233;
double codelKC311;
double codelKC312;
double codelKC313;
double codelKC322;
double codelKC323;
double codelKC333;
double codelTrKA1;
double codelTrKA2;
double codelTrKA3;
double codelTrKB1;
double codelTrKB2;
double codelTrKB3;
double codelTrKC1;
double codelTrKC2;
double codelTrKC3;
double deldelf11;
double deldelf12;
double deldelf13;
double deldelf22;
double deldelf23;
double deldelf33;
double deldelg1111;
double deldelg1112;
double deldelg1113;
double deldelg1122;
double deldelg1123;
double deldelg1133;
double deldelg1211;
double deldelg1212;
double deldelg1213;
double deldelg1222;
double deldelg1223;
double deldelg1233;
double deldelg1311;
double deldelg1312;
double deldelg1313;
double deldelg1322;
double deldelg1323;
double deldelg1333;
double deldelg2211;
double deldelg2212;
double deldelg2213;
double deldelg2222;
double deldelg2223;
double deldelg2233;
double deldelg2311;
double deldelg2312;
double deldelg2313;
double deldelg2322;
double deldelg2323;
double deldelg2333;
double deldelg3311;
double deldelg3312;
double deldelg3313;
double deldelg3322;
double deldelg3323;
double deldelg3333;
double delf1;
double delf2;
double delf3;
double delg111;
double delg112;
double delg113;
double delg122;
double delg123;
double delg133;
double delg211;
double delg212;
double delg213;
double delg222;
double delg223;
double delg233;
double delg311;
double delg312;
double delg313;
double delg322;
double delg323;
double delg333;
double delK111;
double delK112;
double delK113;
double delK122;
double delK123;
double delK133;
double delK211;
double delK212;
double delK213;
double delK222;
double delK223;
double delK233;
double delK311;
double delK312;
double delK313;
double delK322;
double delK323;
double delK333;
double denom1;
double denom2;
double denom3;
double gamma111;
double gamma112;
double gamma113;
double gamma121;
double gamma122;
double gamma123;
double gamma131;
double gamma132;
double gamma133;
double gamma211;
double gamma212;
double gamma213;
double gamma221;
double gamma222;
double gamma223;
double gamma231;
double gamma232;
double gamma233;
double gamma311;
double gamma312;
double gamma313;
double gamma321;
double gamma322;
double gamma323;
double gamma331;
double gamma332;
double gamma333;
double gammado111;
double gammado112;
double gammado113;
double gammado121;
double gammado122;
double gammado123;
double gammado131;
double gammado132;
double gammado133;
double gammado211;
double gammado212;
double gammado213;
double gammado221;
double gammado222;
double gammado223;
double gammado231;
double gammado232;
double gammado233;
double gammado311;
double gammado312;
double gammado313;
double gammado321;
double gammado322;
double gammado323;
double gammado331;
double gammado332;
double gammado333;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double Kud11;
double Kud12;
double Kud13;
double Kud21;
double Kud22;
double Kud23;
double Kud31;
double Kud32;
double Kud33;
double momrhs1;
double momrhs2;
double momrhs3;
double momu1;
double momu2;
double momu3;
double R11;
double R12;
double R13;
double R22;
double R23;
double R33;
double RA11;
double RA12;
double RA13;
double RA22;
double RA23;
double RA33;
double RB11;
double RB12;
double RB13;
double RB22;
double RB23;
double RB33;
double RC11;
double RC12;
double RC13;
double RC22;
double RC23;
double RC33;
double RD11;
double RD12;
double RD13;
double RD22;
double RD23;
double RD33;



/* Jetzt geht's los! */
delg111
=
dg111[ijk]
;

delg112
=
dg121[ijk]
;

delg113
=
dg131[ijk]
;

delg122
=
dg221[ijk]
;

delg123
=
dg231[ijk]
;

delg133
=
dg331[ijk]
;

delg211
=
dg112[ijk]
;

delg212
=
dg122[ijk]
;

delg213
=
dg132[ijk]
;

delg222
=
dg222[ijk]
;

delg223
=
dg232[ijk]
;

delg233
=
dg332[ijk]
;

delg311
=
dg113[ijk]
;

delg312
=
dg123[ijk]
;

delg313
=
dg133[ijk]
;

delg322
=
dg223[ijk]
;

delg323
=
dg233[ijk]
;

delg333
=
dg333[ijk]
;

deldelg1111
=
ddg1111[ijk]
;

deldelg1112
=
ddg1211[ijk]
;

deldelg1113
=
ddg1311[ijk]
;

deldelg1122
=
ddg2211[ijk]
;

deldelg1123
=
ddg2311[ijk]
;

deldelg1133
=
ddg3311[ijk]
;

deldelg1211
=
ddg1112[ijk]
;

deldelg1212
=
ddg1212[ijk]
;

deldelg1213
=
ddg1312[ijk]
;

deldelg1222
=
ddg2212[ijk]
;

deldelg1223
=
ddg2312[ijk]
;

deldelg1233
=
ddg3312[ijk]
;

deldelg1311
=
ddg1113[ijk]
;

deldelg1312
=
ddg1213[ijk]
;

deldelg1313
=
ddg1313[ijk]
;

deldelg1322
=
ddg2213[ijk]
;

deldelg1323
=
ddg2313[ijk]
;

deldelg1333
=
ddg3313[ijk]
;

deldelg2211
=
ddg1122[ijk]
;

deldelg2212
=
ddg1222[ijk]
;

deldelg2213
=
ddg1322[ijk]
;

deldelg2222
=
ddg2222[ijk]
;

deldelg2223
=
ddg2322[ijk]
;

deldelg2233
=
ddg3322[ijk]
;

deldelg2311
=
ddg1123[ijk]
;

deldelg2312
=
ddg1223[ijk]
;

deldelg2313
=
ddg1323[ijk]
;

deldelg2322
=
ddg2223[ijk]
;

deldelg2323
=
ddg2323[ijk]
;

deldelg2333
=
ddg3323[ijk]
;

deldelg3311
=
ddg1133[ijk]
;

deldelg3312
=
ddg1233[ijk]
;

deldelg3313
=
ddg1333[ijk]
;

deldelg3322
=
ddg2233[ijk]
;

deldelg3323
=
ddg2333[ijk]
;

deldelg3333
=
ddg3333[ijk]
;

delK111
=
dK111[ijk]
;

delK112
=
dK121[ijk]
;

delK113
=
dK131[ijk]
;

delK122
=
dK221[ijk]
;

delK123
=
dK231[ijk]
;

delK133
=
dK331[ijk]
;

delK211
=
dK112[ijk]
;

delK212
=
dK122[ijk]
;

delK213
=
dK132[ijk]
;

delK222
=
dK222[ijk]
;

delK223
=
dK232[ijk]
;

delK233
=
dK332[ijk]
;

delK311
=
dK113[ijk]
;

delK312
=
dK123[ijk]
;

delK313
=
dK133[ijk]
;

delK322
=
dK223[ijk]
;

delK323
=
dK233[ijk]
;

delK333
=
dK333[ijk]
;

f
=
Power(psi[ijk],4)
;

delf1
=
4.*f*dpop1[ijk]
;

delf2
=
4.*f*dpop2[ijk]
;

delf3
=
4.*f*dpop3[ijk]
;

deldelf11
=
f*(4.*ddpop11[ijk] + 12.*pow2(dpop1[ijk]))
;

deldelf12
=
f*(4.*ddpop12[ijk] + 12.*dpop1[ijk]*dpop2[ijk])
;

deldelf13
=
f*(4.*ddpop13[ijk] + 12.*dpop1[ijk]*dpop3[ijk])
;

deldelf22
=
f*(4.*ddpop22[ijk] + 12.*pow2(dpop2[ijk]))
;

deldelf23
=
f*(4.*ddpop23[ijk] + 12.*dpop2[ijk]*dpop3[ijk])
;

deldelf33
=
f*(4.*ddpop33[ijk] + 12.*pow2(dpop3[ijk]))
;

deldelg1111
=
2.*delf1*delg111 + deldelg1111*f + deldelf11*g11[ijk]
;

deldelg1112
=
2.*delf1*delg112 + deldelg1112*f + deldelf11*g12[ijk]
;

deldelg1113
=
2.*delf1*delg113 + deldelg1113*f + deldelf11*g13[ijk]
;

deldelg1122
=
2.*delf1*delg122 + deldelg1122*f + deldelf11*g22[ijk]
;

deldelg1123
=
2.*delf1*delg123 + deldelg1123*f + deldelf11*g23[ijk]
;

deldelg1133
=
2.*delf1*delg133 + deldelg1133*f + deldelf11*g33[ijk]
;

deldelg1211
=
delf2*delg111 + delf1*delg211 + deldelg1211*f + deldelf12*g11[ijk]
;

deldelg1212
=
delf2*delg112 + delf1*delg212 + deldelg1212*f + deldelf12*g12[ijk]
;

deldelg1213
=
delf2*delg113 + delf1*delg213 + deldelg1213*f + deldelf12*g13[ijk]
;

deldelg1222
=
delf2*delg122 + delf1*delg222 + deldelg1222*f + deldelf12*g22[ijk]
;

deldelg1223
=
delf2*delg123 + delf1*delg223 + deldelg1223*f + deldelf12*g23[ijk]
;

deldelg1233
=
delf2*delg133 + delf1*delg233 + deldelg1233*f + deldelf12*g33[ijk]
;

deldelg1311
=
delf3*delg111 + delf1*delg311 + deldelg1311*f + deldelf13*g11[ijk]
;

deldelg1312
=
delf3*delg112 + delf1*delg312 + deldelg1312*f + deldelf13*g12[ijk]
;

deldelg1313
=
delf3*delg113 + delf1*delg313 + deldelg1313*f + deldelf13*g13[ijk]
;

deldelg1322
=
delf3*delg122 + delf1*delg322 + deldelg1322*f + deldelf13*g22[ijk]
;

deldelg1323
=
delf3*delg123 + delf1*delg323 + deldelg1323*f + deldelf13*g23[ijk]
;

deldelg1333
=
delf3*delg133 + delf1*delg333 + deldelg1333*f + deldelf13*g33[ijk]
;

deldelg2211
=
2.*delf2*delg211 + deldelg2211*f + deldelf22*g11[ijk]
;

deldelg2212
=
2.*delf2*delg212 + deldelg2212*f + deldelf22*g12[ijk]
;

deldelg2213
=
2.*delf2*delg213 + deldelg2213*f + deldelf22*g13[ijk]
;

deldelg2222
=
2.*delf2*delg222 + deldelg2222*f + deldelf22*g22[ijk]
;

deldelg2223
=
2.*delf2*delg223 + deldelg2223*f + deldelf22*g23[ijk]
;

deldelg2233
=
2.*delf2*delg233 + deldelg2233*f + deldelf22*g33[ijk]
;

deldelg2311
=
delf3*delg211 + delf2*delg311 + deldelg2311*f + deldelf23*g11[ijk]
;

deldelg2312
=
delf3*delg212 + delf2*delg312 + deldelg2312*f + deldelf23*g12[ijk]
;

deldelg2313
=
delf3*delg213 + delf2*delg313 + deldelg2313*f + deldelf23*g13[ijk]
;

deldelg2322
=
delf3*delg222 + delf2*delg322 + deldelg2322*f + deldelf23*g22[ijk]
;

deldelg2323
=
delf3*delg223 + delf2*delg323 + deldelg2323*f + deldelf23*g23[ijk]
;

deldelg2333
=
delf3*delg233 + delf2*delg333 + deldelg2333*f + deldelf23*g33[ijk]
;

deldelg3311
=
2.*delf3*delg311 + deldelg3311*f + deldelf33*g11[ijk]
;

deldelg3312
=
2.*delf3*delg312 + deldelg3312*f + deldelf33*g12[ijk]
;

deldelg3313
=
2.*delf3*delg313 + deldelg3313*f + deldelf33*g13[ijk]
;

deldelg3322
=
2.*delf3*delg322 + deldelg3322*f + deldelf33*g22[ijk]
;

deldelg3323
=
2.*delf3*delg323 + deldelg3323*f + deldelf33*g23[ijk]
;

deldelg3333
=
2.*delf3*delg333 + deldelg3333*f + deldelf33*g33[ijk]
;

delg111
=
delg111*f + delf1*g11[ijk]
;

delg112
=
delg112*f + delf1*g12[ijk]
;

delg113
=
delg113*f + delf1*g13[ijk]
;

delg122
=
delg122*f + delf1*g22[ijk]
;

delg123
=
delg123*f + delf1*g23[ijk]
;

delg133
=
delg133*f + delf1*g33[ijk]
;

delg211
=
delg211*f + delf2*g11[ijk]
;

delg212
=
delg212*f + delf2*g12[ijk]
;

delg213
=
delg213*f + delf2*g13[ijk]
;

delg222
=
delg222*f + delf2*g22[ijk]
;

delg223
=
delg223*f + delf2*g23[ijk]
;

delg233
=
delg233*f + delf2*g33[ijk]
;

delg311
=
delg311*f + delf3*g11[ijk]
;

delg312
=
delg312*f + delf3*g12[ijk]
;

delg313
=
delg313*f + delf3*g13[ijk]
;

delg322
=
delg322*f + delf3*g22[ijk]
;

delg323
=
delg323*f + delf3*g23[ijk]
;

delg333
=
delg333*f + delf3*g33[ijk]
;

detginvf
=
1/(f*(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
      g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
      g11[ijk]*pow2(g23[ijk])))
;

ginv11
=
detginvf*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

ginv12
=
detginvf*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

ginv13
=
detginvf*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

ginv22
=
detginvf*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

ginv23
=
detginvf*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

ginv33
=
detginvf*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
;

gammado111
=
0.5*delg111
;

gammado112
=
0.5*delg211
;

gammado113
=
0.5*delg311
;

gammado121
=
0.5*delg211
;

gammado122
=
-0.5*delg122 + delg212
;

gammado123
=
0.5*(-delg123 + delg213 + delg312)
;

gammado131
=
0.5*delg311
;

gammado132
=
0.5*(-delg123 + delg213 + delg312)
;

gammado133
=
-0.5*delg133 + delg313
;

gammado211
=
delg112 - 0.5*delg211
;

gammado212
=
0.5*delg122
;

gammado213
=
0.5*(delg123 - delg213 + delg312)
;

gammado221
=
0.5*delg122
;

gammado222
=
0.5*delg222
;

gammado223
=
0.5*delg322
;

gammado231
=
0.5*(delg123 - delg213 + delg312)
;

gammado232
=
0.5*delg322
;

gammado233
=
-0.5*delg233 + delg323
;

gammado311
=
delg113 - 0.5*delg311
;

gammado312
=
0.5*(delg123 + delg213 - delg312)
;

gammado313
=
0.5*delg133
;

gammado321
=
0.5*(delg123 + delg213 - delg312)
;

gammado322
=
delg223 - 0.5*delg322
;

gammado323
=
0.5*delg233
;

gammado331
=
0.5*delg133
;

gammado332
=
0.5*delg233
;

gammado333
=
0.5*delg333
;

gamma111
=
gammado111*ginv11 + gammado211*ginv12 + gammado311*ginv13
;

gamma112
=
gammado112*ginv11 + gammado212*ginv12 + gammado312*ginv13
;

gamma113
=
gammado113*ginv11 + gammado213*ginv12 + gammado313*ginv13
;

gamma121
=
gammado121*ginv11 + gammado221*ginv12 + gammado321*ginv13
;

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
;

gamma131
=
gammado131*ginv11 + gammado231*ginv12 + gammado331*ginv13
;

gamma132
=
gammado132*ginv11 + gammado232*ginv12 + gammado332*ginv13
;

gamma133
=
gammado133*ginv11 + gammado233*ginv12 + gammado333*ginv13
;

gamma211
=
gammado111*ginv12 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv12 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv12 + gammado213*ginv22 + gammado313*ginv23
;

gamma221
=
gammado121*ginv12 + gammado221*ginv22 + gammado321*ginv23
;

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
;

gamma231
=
gammado131*ginv12 + gammado231*ginv22 + gammado331*ginv23
;

gamma232
=
gammado132*ginv12 + gammado232*ginv22 + gammado332*ginv23
;

gamma233
=
gammado133*ginv12 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv13 + gammado211*ginv23 + gammado311*ginv33
;

gamma312
=
gammado112*ginv13 + gammado212*ginv23 + gammado312*ginv33
;

gamma313
=
gammado113*ginv13 + gammado213*ginv23 + gammado313*ginv33
;

gamma321
=
gammado121*ginv13 + gammado221*ginv23 + gammado321*ginv33
;

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma331
=
gammado131*ginv13 + gammado231*ginv23 + gammado331*ginv33
;

gamma332
=
gammado132*ginv13 + gammado232*ginv23 + gammado332*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

R11
=
(gamma112*gammado111 - gamma111*gammado121 + gamma212*gammado211 - 
     gamma211*gammado221 + gamma312*gammado311 - gamma311*gammado321)*ginv12 \
+ (gamma113*gammado111 - gamma111*gammado131 + gamma213*gammado211 - 
     gamma211*gammado231 + gamma313*gammado311 - gamma311*gammado331)*ginv13 \
+ (deldelg1212 + 0.5*(-deldelg1122 - deldelg2211) + gamma112*gammado112 - 
     gamma111*gammado122 + gamma212*gammado212 - gamma211*gammado222 + 
     gamma312*gammado312 - gamma311*gammado322)*ginv22 + 
  (-deldelg1123 + deldelg1213 + deldelg1312 - deldelg2311 + 
     gamma113*gammado112 + gamma112*gammado113 - 
     gamma111*(gammado123 + gammado132) + gamma213*gammado212 + 
     gamma212*gammado213 - gamma211*(gammado223 + gammado232) + 
     gamma313*gammado312 + gamma312*gammado313 - 
     gamma311*(gammado323 + gammado332))*ginv23 + 
  (deldelg1313 + 0.5*(-deldelg1133 - deldelg3311) + gamma113*gammado113 - 
     gamma111*gammado133 + gamma213*gammado213 - gamma211*gammado233 + 
     gamma313*gammado313 - gamma311*gammado333)*ginv33
;

R12
=
(-(gamma112*gammado111) + gamma111*gammado121 - gamma212*gammado211 + 
     gamma211*gammado221 - gamma312*gammado311 + gamma311*gammado321)*ginv11 \
+ (-deldelg1212 + 0.5*(deldelg1122 + deldelg2211) - gamma112*gammado112 + 
     gamma111*gammado122 - gamma212*gammado212 + gamma211*gammado222 - 
     gamma312*gammado312 + gamma311*gammado322)*ginv12 + 
  (0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) + 
     gamma113*gammado121 + gamma111*gammado123 - 
     gamma112*(gammado113 + gammado131) + gamma213*gammado221 + 
     gamma211*gammado223 - gamma212*(gammado213 + gammado231) + 
     gamma313*gammado321 + gamma311*gammado323 - 
     gamma312*(gammado313 + gammado331))*ginv13 + 
  (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) + 
     gamma113*gammado122 - gamma112*gammado132 + gamma213*gammado222 - 
     gamma212*gammado232 + gamma313*gammado322 - gamma312*gammado332)*ginv23 \
+ (0.5*(-deldelg1233 + deldelg1323 + deldelg2313 - deldelg3312) + 
     gamma113*gammado123 - gamma112*gammado133 + gamma213*gammado223 - 
     gamma212*gammado233 + gamma313*gammado323 - gamma312*gammado333)*ginv33
;

R12
=
(-deldelg1212 + 0.5*(deldelg1122 + deldelg2211) + gamma122*gammado111 - 
     gamma121*gammado121 + gamma222*gammado211 - gamma221*gammado221 + 
     gamma322*gammado311 - gamma321*gammado321)*ginv12 + 
  (0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) + 
     gamma123*gammado111 - gamma121*gammado131 + gamma223*gammado211 - 
     gamma221*gammado231 + gamma323*gammado311 - gamma321*gammado331)*ginv13 \
+ (gamma122*gammado112 - gamma121*gammado122 + gamma222*gammado212 - 
     gamma221*gammado222 + gamma322*gammado312 - gamma321*gammado322)*ginv22 \
+ (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) + 
     gamma123*gammado112 + gamma122*gammado113 - 
     gamma121*(gammado123 + gammado132) + gamma223*gammado212 + 
     gamma222*gammado213 - gamma221*(gammado223 + gammado232) + 
     gamma323*gammado312 + gamma322*gammado313 - 
     gamma321*(gammado323 + gammado332))*ginv23 + 
  (0.5*(-deldelg1233 + deldelg1323 + deldelg2313 - deldelg3312) + 
     gamma123*gammado113 - gamma121*gammado133 + gamma223*gammado213 - 
     gamma221*gammado233 + gamma323*gammado313 - gamma321*gammado333)*ginv33
;

R13
=
(-(gamma113*gammado111) + gamma111*gammado131 - gamma213*gammado211 + 
     gamma211*gammado231 - gamma313*gammado311 + gamma311*gammado331)*ginv11 \
+ (0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) - 
     gamma113*(gammado112 + gammado121) + gamma112*gammado131 + 
     gamma111*gammado132 - gamma213*(gammado212 + gammado221) + 
     gamma212*gammado231 + gamma211*gammado232 - 
     gamma313*(gammado312 + gammado321) + gamma312*gammado331 + 
     gamma311*gammado332)*ginv12 + 
  (-deldelg1313 + 0.5*(deldelg1133 + deldelg3311) - gamma113*gammado113 + 
     gamma111*gammado133 - gamma213*gammado213 + gamma211*gammado233 - 
     gamma313*gammado313 + gamma311*gammado333)*ginv13 + 
  (0.5*(deldelg1223 - deldelg1322 - deldelg2213 + deldelg2312) - 
     gamma113*gammado122 + gamma112*gammado132 - gamma213*gammado222 + 
     gamma212*gammado232 - gamma313*gammado322 + gamma312*gammado332)*ginv22 \
+ (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) - 
     gamma113*gammado123 + gamma112*gammado133 - gamma213*gammado223 + 
     gamma212*gammado233 - gamma313*gammado323 + gamma312*gammado333)*ginv23
;

R13
=
(0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) + 
     gamma132*gammado111 - gamma131*gammado121 + gamma232*gammado211 - 
     gamma231*gammado221 + gamma332*gammado311 - gamma331*gammado321)*ginv12 \
+ (-deldelg1313 + 0.5*(deldelg1133 + deldelg3311) + gamma133*gammado111 - 
     gamma131*gammado131 + gamma233*gammado211 - gamma231*gammado231 + 
     gamma333*gammado311 - gamma331*gammado331)*ginv13 + 
  (0.5*(deldelg1223 - deldelg1322 - deldelg2213 + deldelg2312) + 
     gamma132*gammado112 - gamma131*gammado122 + gamma232*gammado212 - 
     gamma231*gammado222 + gamma332*gammado312 - gamma331*gammado322)*ginv22 \
+ (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) + 
     gamma133*gammado112 + gamma132*gammado113 - 
     gamma131*(gammado123 + gammado132) + gamma233*gammado212 + 
     gamma232*gammado213 - gamma231*(gammado223 + gammado232) + 
     gamma333*gammado312 + gamma332*gammado313 - 
     gamma331*(gammado323 + gammado332))*ginv23 + 
  (gamma133*gammado113 - gamma131*gammado133 + gamma233*gammado213 - 
     gamma231*gammado233 + gamma333*gammado313 - gamma331*gammado333)*ginv33
;

R22
=
(deldelg1212 + 0.5*(-deldelg1122 - deldelg2211) - gamma122*gammado111 + 
     gamma121*gammado121 - gamma222*gammado211 + gamma221*gammado221 - 
     gamma322*gammado311 + gamma321*gammado321)*ginv11 + 
  (-(gamma122*gammado112) + gamma121*gammado122 - gamma222*gammado212 + 
     gamma221*gammado222 - gamma322*gammado312 + gamma321*gammado322)*ginv12 \
+ (deldelg1223 - deldelg1322 - deldelg2213 + deldelg2312 + 
     gamma123*gammado121 + gamma121*gammado123 - 
     gamma122*(gammado113 + gammado131) + gamma223*gammado221 + 
     gamma221*gammado223 - gamma222*(gammado213 + gammado231) + 
     gamma323*gammado321 + gamma321*gammado323 - 
     gamma322*(gammado313 + gammado331))*ginv13 + 
  (gamma123*gammado122 - gamma122*gammado132 + gamma223*gammado222 - 
     gamma222*gammado232 + gamma323*gammado322 - gamma322*gammado332)*ginv23 \
+ (deldelg2323 + 0.5*(-deldelg2233 - deldelg3322) + gamma123*gammado123 - 
     gamma122*gammado133 + gamma223*gammado223 - gamma222*gammado233 + 
     gamma323*gammado323 - gamma322*gammado333)*ginv33
;

R23
=
(0.5*(-deldelg1123 + deldelg1213 + deldelg1312 - deldelg2311) - 
     gamma123*gammado111 + gamma121*gammado131 - gamma223*gammado211 + 
     gamma221*gammado231 - gamma323*gammado311 + gamma321*gammado331)*ginv11 \
+ (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) - 
     gamma123*(gammado112 + gammado121) + gamma122*gammado131 + 
     gamma121*gammado132 - gamma223*(gammado212 + gammado221) + 
     gamma222*gammado231 + gamma221*gammado232 - 
     gamma323*(gammado312 + gammado321) + gamma322*gammado331 + 
     gamma321*gammado332)*ginv12 + 
  (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) - 
     gamma123*gammado113 + gamma121*gammado133 - gamma223*gammado213 + 
     gamma221*gammado233 - gamma323*gammado313 + gamma321*gammado333)*ginv13 \
+ (-(gamma123*gammado122) + gamma122*gammado132 - gamma223*gammado222 + 
     gamma222*gammado232 - gamma323*gammado322 + gamma322*gammado332)*ginv22 \
+ (-deldelg2323 + 0.5*(deldelg2233 + deldelg3322) - gamma123*gammado123 + 
     gamma122*gammado133 - gamma223*gammado223 + gamma222*gammado233 - 
     gamma323*gammado323 + gamma322*gammado333)*ginv23
;

R23
=
(0.5*(-deldelg1123 + deldelg1213 + deldelg1312 - deldelg2311) - 
     gamma132*gammado111 + gamma131*gammado121 - gamma232*gammado211 + 
     gamma231*gammado221 - gamma332*gammado311 + gamma331*gammado321)*ginv11 \
+ (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) - 
     gamma132*gammado112 + gamma131*gammado122 - gamma232*gammado212 + 
     gamma231*gammado222 - gamma332*gammado312 + gamma331*gammado322)*ginv12 \
+ (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) + 
     gamma133*gammado121 + gamma131*gammado123 - 
     gamma132*(gammado113 + gammado131) + gamma233*gammado221 + 
     gamma231*gammado223 - gamma232*(gammado213 + gammado231) + 
     gamma333*gammado321 + gamma331*gammado323 - 
     gamma332*(gammado313 + gammado331))*ginv13 + 
  (-deldelg2323 + 0.5*(deldelg2233 + deldelg3322) + gamma133*gammado122 - 
     gamma132*gammado132 + gamma233*gammado222 - gamma232*gammado232 + 
     gamma333*gammado322 - gamma332*gammado332)*ginv23 + 
  (gamma133*gammado123 - gamma132*gammado133 + gamma233*gammado223 - 
     gamma232*gammado233 + gamma333*gammado323 - gamma332*gammado333)*ginv33
;

R33
=
(deldelg1313 + 0.5*(-deldelg1133 - deldelg3311) - gamma133*gammado111 + 
     gamma131*gammado131 - gamma233*gammado211 + gamma231*gammado231 - 
     gamma333*gammado311 + gamma331*gammado331)*ginv11 + 
  (-deldelg1233 + deldelg1323 + deldelg2313 - deldelg3312 - 
     gamma133*(gammado112 + gammado121) + gamma132*gammado131 + 
     gamma131*gammado132 - gamma233*(gammado212 + gammado221) + 
     gamma232*gammado231 + gamma231*gammado232 - 
     gamma333*(gammado312 + gammado321) + gamma332*gammado331 + 
     gamma331*gammado332)*ginv12 + 
  (-(gamma133*gammado113) + gamma131*gammado133 - gamma233*gammado213 + 
     gamma231*gammado233 - gamma333*gammado313 + gamma331*gammado333)*ginv13 \
+ (deldelg2323 + 0.5*(-deldelg2233 - deldelg3322) - gamma133*gammado122 + 
     gamma132*gammado132 - gamma233*gammado222 + gamma232*gammado232 - 
     gamma333*gammado322 + gamma332*gammado332)*ginv22 + 
  (-(gamma133*gammado123) + gamma132*gammado133 - gamma233*gammado223 + 
     gamma232*gammado233 - gamma333*gammado323 + gamma332*gammado333)*ginv23
;

R
=
ginv11*R11 + ginv22*R22 + 2.*(ginv12*R12 + ginv13*R13 + ginv23*R23) + 
  ginv33*R33
;

Kud11
=
ginv11*K11[ijk] + ginv12*K12[ijk] + ginv13*K13[ijk]
;

Kud12
=
ginv11*K12[ijk] + ginv12*K22[ijk] + ginv13*K23[ijk]
;

Kud13
=
ginv11*K13[ijk] + ginv12*K23[ijk] + ginv13*K33[ijk]
;

Kud21
=
ginv12*K11[ijk] + ginv22*K12[ijk] + ginv23*K13[ijk]
;

Kud22
=
ginv12*K12[ijk] + ginv22*K22[ijk] + ginv23*K23[ijk]
;

Kud23
=
ginv12*K13[ijk] + ginv22*K23[ijk] + ginv23*K33[ijk]
;

Kud31
=
ginv13*K11[ijk] + ginv23*K12[ijk] + ginv33*K13[ijk]
;

Kud32
=
ginv13*K12[ijk] + ginv23*K22[ijk] + ginv33*K23[ijk]
;

Kud33
=
ginv13*K13[ijk] + ginv23*K23[ijk] + ginv33*K33[ijk]
;

K
=
Kud11 + Kud22 + Kud33
;

KudKud
=
2.*(Kud12*Kud21 + Kud13*Kud31 + Kud23*Kud32) + pow2(Kud11) + pow2(Kud22) + 
  pow2(Kud33)
;

codelK111
=
delK111 - 2.*(gamma111*K11[ijk] + gamma211*K12[ijk] + gamma311*K13[ijk])
;

codelK112
=
delK112 - gamma112*K11[ijk] - (gamma111 + gamma212)*K12[ijk] - 
  gamma312*K13[ijk] - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelK113
=
delK113 - gamma113*K11[ijk] - gamma213*K12[ijk] - 
  (gamma111 + gamma313)*K13[ijk] - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelK122
=
delK122 - 2.*(gamma112*K12[ijk] + gamma212*K22[ijk] + gamma312*K23[ijk])
;

codelK123
=
delK123 - gamma113*K12[ijk] - gamma112*K13[ijk] - gamma213*K22[ijk] - 
  (gamma212 + gamma313)*K23[ijk] - gamma312*K33[ijk]
;

codelK133
=
delK133 - 2.*(gamma113*K13[ijk] + gamma213*K23[ijk] + gamma313*K33[ijk])
;

codelK211
=
delK211 - 2.*(gamma121*K11[ijk] + gamma221*K12[ijk] + gamma321*K13[ijk])
;

codelK212
=
delK212 - gamma122*K11[ijk] - (gamma121 + gamma222)*K12[ijk] - 
  gamma322*K13[ijk] - gamma221*K22[ijk] - gamma321*K23[ijk]
;

codelK213
=
delK213 - gamma123*K11[ijk] - gamma223*K12[ijk] - 
  (gamma121 + gamma323)*K13[ijk] - gamma221*K23[ijk] - gamma321*K33[ijk]
;

codelK222
=
delK222 - 2.*(gamma122*K12[ijk] + gamma222*K22[ijk] + gamma322*K23[ijk])
;

codelK223
=
delK223 - gamma123*K12[ijk] - gamma122*K13[ijk] - gamma223*K22[ijk] - 
  (gamma222 + gamma323)*K23[ijk] - gamma322*K33[ijk]
;

codelK233
=
delK233 - 2.*(gamma123*K13[ijk] + gamma223*K23[ijk] + gamma323*K33[ijk])
;

codelK311
=
delK311 - 2.*(gamma131*K11[ijk] + gamma231*K12[ijk] + gamma331*K13[ijk])
;

codelK312
=
delK312 - gamma132*K11[ijk] - (gamma131 + gamma232)*K12[ijk] - 
  gamma332*K13[ijk] - gamma231*K22[ijk] - gamma331*K23[ijk]
;

codelK313
=
delK313 - gamma133*K11[ijk] - gamma233*K12[ijk] - 
  (gamma131 + gamma333)*K13[ijk] - gamma231*K23[ijk] - gamma331*K33[ijk]
;

codelK322
=
delK322 - 2.*(gamma132*K12[ijk] + gamma232*K22[ijk] + gamma332*K23[ijk])
;

codelK323
=
delK323 - gamma133*K12[ijk] - gamma132*K13[ijk] - gamma233*K22[ijk] - 
  (gamma232 + gamma333)*K23[ijk] - gamma332*K33[ijk]
;

codelK333
=
delK333 - 2.*(gamma133*K13[ijk] + gamma233*K23[ijk] + gamma333*K33[ijk])
;

cdKudd111
=
codelK111*ginv11 + codelK211*ginv12 + codelK311*ginv13
;

cdKudd112
=
codelK112*ginv11 + codelK212*ginv12 + codelK312*ginv13
;

cdKudd113
=
codelK113*ginv11 + codelK213*ginv12 + codelK313*ginv13
;

cdKudd122
=
codelK122*ginv11 + codelK222*ginv12 + codelK322*ginv13
;

cdKudd123
=
codelK123*ginv11 + codelK223*ginv12 + codelK323*ginv13
;

cdKudd133
=
codelK133*ginv11 + codelK233*ginv12 + codelK333*ginv13
;

cdKudd211
=
codelK111*ginv12 + codelK211*ginv22 + codelK311*ginv23
;

cdKudd212
=
codelK112*ginv12 + codelK212*ginv22 + codelK312*ginv23
;

cdKudd213
=
codelK113*ginv12 + codelK213*ginv22 + codelK313*ginv23
;

cdKudd222
=
codelK122*ginv12 + codelK222*ginv22 + codelK322*ginv23
;

cdKudd223
=
codelK123*ginv12 + codelK223*ginv22 + codelK323*ginv23
;

cdKudd233
=
codelK133*ginv12 + codelK233*ginv22 + codelK333*ginv23
;

cdKudd311
=
codelK111*ginv13 + codelK211*ginv23 + codelK311*ginv33
;

cdKudd312
=
codelK112*ginv13 + codelK212*ginv23 + codelK312*ginv33
;

cdKudd313
=
codelK113*ginv13 + codelK213*ginv23 + codelK313*ginv33
;

cdKudd322
=
codelK122*ginv13 + codelK222*ginv23 + codelK322*ginv33
;

cdKudd323
=
codelK123*ginv13 + codelK223*ginv23 + codelK323*ginv33
;

cdKudd333
=
codelK133*ginv13 + codelK233*ginv23 + codelK333*ginv33
;


if(rho==NULL) { 

hamrhs
=
0
;


} else { 

hamrhs
=
50.26548245743669182*rho[ijk]
;


} 


if(j1==NULL) { 

momrhs1
=
0
;

momrhs2
=
0
;

momrhs3
=
0
;


} else { 

momrhs1
=
-25.132741228718345908*j1[ijk]
;

momrhs2
=
-25.132741228718345908*j2[ijk]
;

momrhs3
=
-25.132741228718345908*j3[ijk]
;


} 

ham[ijk]
=
-hamrhs - KudKud + R + pow2(K)
;

momu1
=
(cdKudd212 + cdKudd313)*ginv11 + 
  (-cdKudd112 + cdKudd222 + cdKudd323)*ginv12 + 
  (-cdKudd113 + cdKudd223 + cdKudd333)*ginv13 - cdKudd122*ginv22 - 
  2.*cdKudd123*ginv23 - cdKudd133*ginv33
;

momu2
=
-(cdKudd211*ginv11) + (cdKudd111 - cdKudd212 + cdKudd313)*ginv12 - 
  2.*cdKudd213*ginv13 + (cdKudd112 + cdKudd323)*ginv22 + 
  (cdKudd113 - cdKudd223 + cdKudd333)*ginv23 - cdKudd233*ginv33
;

momu3
=
-(cdKudd311*ginv11) - 2.*cdKudd312*ginv12 + 
  (cdKudd111 + cdKudd212 - cdKudd313)*ginv13 - cdKudd322*ginv22 + 
  (cdKudd112 + cdKudd222 - cdKudd323)*ginv23 + (cdKudd113 + cdKudd223)*ginv33
;

mom1[ijk]
=
-momrhs1 + f*(momu1*g11[ijk] + momu2*g12[ijk] + momu3*g13[ijk])
;

mom2[ijk]
=
-momrhs2 + f*(momu1*g12[ijk] + momu2*g22[ijk] + momu3*g23[ijk])
;

mom3[ijk]
=
-momrhs3 + f*(momu1*g13[ijk] + momu2*g23[ijk] + momu3*g33[ijk])
;

trK[ijk]
=
K
;



/* conditional */
if (normConstr) {



/* conditional */
if (TermByTerm) {

RA11
=
-(deldelg1111*ginv11) - deldelg2211*ginv22 - 
  2.*(deldelg1211*ginv12 + deldelg1311*ginv13 + deldelg2311*ginv23) - 
  deldelg3311*ginv33
;

RA12
=
-(deldelg1112*ginv11) - deldelg2212*ginv22 - 
  2.*(deldelg1212*ginv12 + deldelg1312*ginv13 + deldelg2312*ginv23) - 
  deldelg3312*ginv33
;

RA13
=
-(deldelg1113*ginv11) - deldelg2213*ginv22 - 
  2.*(deldelg1213*ginv12 + deldelg1313*ginv13 + deldelg2313*ginv23) - 
  deldelg3313*ginv33
;

RA22
=
-(deldelg1122*ginv11) - deldelg2222*ginv22 - 
  2.*(deldelg1222*ginv12 + deldelg1322*ginv13 + deldelg2322*ginv23) - 
  deldelg3322*ginv33
;

RA23
=
-(deldelg1123*ginv11) - deldelg2223*ginv22 - 
  2.*(deldelg1223*ginv12 + deldelg1323*ginv13 + deldelg2323*ginv23) - 
  deldelg3323*ginv33
;

RA33
=
-(deldelg1133*ginv11) - deldelg2233*ginv22 - 
  2.*(deldelg1233*ginv12 + deldelg1333*ginv13 + deldelg2333*ginv23) - 
  deldelg3333*ginv33
;

RB11
=
deldelg1111*ginv11 + (deldelg1112 + deldelg1211)*ginv12 + 
  (deldelg1113 + deldelg1311)*ginv13 + deldelg1212*ginv22 + 
  (deldelg1213 + deldelg1312)*ginv23 + deldelg1313*ginv33
;

RB12
=
deldelg1112*ginv11 + (deldelg1122 + deldelg1212)*ginv12 + 
  (deldelg1123 + deldelg1312)*ginv13 + deldelg1222*ginv22 + 
  (deldelg1223 + deldelg1322)*ginv23 + deldelg1323*ginv33
;

RB12
=
deldelg1211*ginv11 + (deldelg1212 + deldelg2211)*ginv12 + 
  (deldelg1213 + deldelg2311)*ginv13 + deldelg2212*ginv22 + 
  (deldelg2213 + deldelg2312)*ginv23 + deldelg2313*ginv33
;

RB13
=
deldelg1113*ginv11 + (deldelg1123 + deldelg1213)*ginv12 + 
  (deldelg1133 + deldelg1313)*ginv13 + deldelg1223*ginv22 + 
  (deldelg1233 + deldelg1323)*ginv23 + deldelg1333*ginv33
;

RB13
=
deldelg1311*ginv11 + (deldelg1312 + deldelg2311)*ginv12 + 
  (deldelg1313 + deldelg3311)*ginv13 + deldelg2312*ginv22 + 
  (deldelg2313 + deldelg3312)*ginv23 + deldelg3313*ginv33
;

RB22
=
deldelg1212*ginv11 + (deldelg1222 + deldelg2212)*ginv12 + 
  (deldelg1223 + deldelg2312)*ginv13 + deldelg2222*ginv22 + 
  (deldelg2223 + deldelg2322)*ginv23 + deldelg2323*ginv33
;

RB23
=
deldelg1213*ginv11 + (deldelg1223 + deldelg2213)*ginv12 + 
  (deldelg1233 + deldelg2313)*ginv13 + deldelg2223*ginv22 + 
  (deldelg2233 + deldelg2323)*ginv23 + deldelg2333*ginv33
;

RB23
=
deldelg1312*ginv11 + (deldelg1322 + deldelg2312)*ginv12 + 
  (deldelg1323 + deldelg3312)*ginv13 + deldelg2322*ginv22 + 
  (deldelg2323 + deldelg3322)*ginv23 + deldelg3323*ginv33
;

RB33
=
deldelg1313*ginv11 + (deldelg1323 + deldelg2313)*ginv12 + 
  (deldelg1333 + deldelg3313)*ginv13 + deldelg2323*ginv22 + 
  (deldelg2333 + deldelg3323)*ginv23 + deldelg3333*ginv33
;

RC11
=
(gamma111*gammado111 + gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (gamma112*gammado111 + gamma111*gammado112 + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (gamma113*gammado111 + gamma111*gammado113 + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (gamma112*gammado112 + gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (gamma113*gammado112 + gamma112*gammado113 + gamma213*gammado212 + 
     gamma212*gammado213 + gamma313*gammado312 + gamma312*gammado313)*ginv23 \
+ (gamma113*gammado113 + gamma213*gammado213 + gamma313*gammado313)*ginv33
;

RC12
=
(gamma121*gammado111 + gamma221*gammado211 + gamma321*gammado311)*ginv11 + 
  (gamma122*gammado111 + gamma121*gammado112 + gamma222*gammado211 + 
     gamma221*gammado212 + gamma322*gammado311 + gamma321*gammado312)*ginv12 \
+ (gamma123*gammado111 + gamma121*gammado113 + gamma223*gammado211 + 
     gamma221*gammado213 + gamma323*gammado311 + gamma321*gammado313)*ginv13 \
+ (gamma122*gammado112 + gamma222*gammado212 + gamma322*gammado312)*ginv22 + 
  (gamma123*gammado112 + gamma122*gammado113 + gamma223*gammado212 + 
     gamma222*gammado213 + gamma323*gammado312 + gamma322*gammado313)*ginv23 \
+ (gamma123*gammado113 + gamma223*gammado213 + gamma323*gammado313)*ginv33
;

RC12
=
(gamma111*gammado121 + gamma211*gammado221 + gamma311*gammado321)*ginv11 + 
  (gamma112*gammado121 + gamma111*gammado122 + gamma212*gammado221 + 
     gamma211*gammado222 + gamma312*gammado321 + gamma311*gammado322)*ginv12 \
+ (gamma113*gammado121 + gamma111*gammado123 + gamma213*gammado221 + 
     gamma211*gammado223 + gamma313*gammado321 + gamma311*gammado323)*ginv13 \
+ (gamma112*gammado122 + gamma212*gammado222 + gamma312*gammado322)*ginv22 + 
  (gamma113*gammado122 + gamma112*gammado123 + gamma213*gammado222 + 
     gamma212*gammado223 + gamma313*gammado322 + gamma312*gammado323)*ginv23 \
+ (gamma113*gammado123 + gamma213*gammado223 + gamma313*gammado323)*ginv33
;

RC13
=
(gamma131*gammado111 + gamma231*gammado211 + gamma331*gammado311)*ginv11 + 
  (gamma132*gammado111 + gamma131*gammado112 + gamma232*gammado211 + 
     gamma231*gammado212 + gamma332*gammado311 + gamma331*gammado312)*ginv12 \
+ (gamma133*gammado111 + gamma131*gammado113 + gamma233*gammado211 + 
     gamma231*gammado213 + gamma333*gammado311 + gamma331*gammado313)*ginv13 \
+ (gamma132*gammado112 + gamma232*gammado212 + gamma332*gammado312)*ginv22 + 
  (gamma133*gammado112 + gamma132*gammado113 + gamma233*gammado212 + 
     gamma232*gammado213 + gamma333*gammado312 + gamma332*gammado313)*ginv23 \
+ (gamma133*gammado113 + gamma233*gammado213 + gamma333*gammado313)*ginv33
;

RC13
=
(gamma111*gammado131 + gamma211*gammado231 + gamma311*gammado331)*ginv11 + 
  (gamma112*gammado131 + gamma111*gammado132 + gamma212*gammado231 + 
     gamma211*gammado232 + gamma312*gammado331 + gamma311*gammado332)*ginv12 \
+ (gamma113*gammado131 + gamma111*gammado133 + gamma213*gammado231 + 
     gamma211*gammado233 + gamma313*gammado331 + gamma311*gammado333)*ginv13 \
+ (gamma112*gammado132 + gamma212*gammado232 + gamma312*gammado332)*ginv22 + 
  (gamma113*gammado132 + gamma112*gammado133 + gamma213*gammado232 + 
     gamma212*gammado233 + gamma313*gammado332 + gamma312*gammado333)*ginv23 \
+ (gamma113*gammado133 + gamma213*gammado233 + gamma313*gammado333)*ginv33
;

RC22
=
(gamma121*gammado121 + gamma221*gammado221 + gamma321*gammado321)*ginv11 + 
  (gamma122*gammado121 + gamma121*gammado122 + gamma222*gammado221 + 
     gamma221*gammado222 + gamma322*gammado321 + gamma321*gammado322)*ginv12 \
+ (gamma123*gammado121 + gamma121*gammado123 + gamma223*gammado221 + 
     gamma221*gammado223 + gamma323*gammado321 + gamma321*gammado323)*ginv13 \
+ (gamma122*gammado122 + gamma222*gammado222 + gamma322*gammado322)*ginv22 + 
  (gamma123*gammado122 + gamma122*gammado123 + gamma223*gammado222 + 
     gamma222*gammado223 + gamma323*gammado322 + gamma322*gammado323)*ginv23 \
+ (gamma123*gammado123 + gamma223*gammado223 + gamma323*gammado323)*ginv33
;

RC23
=
(gamma131*gammado121 + gamma231*gammado221 + gamma331*gammado321)*ginv11 + 
  (gamma132*gammado121 + gamma131*gammado122 + gamma232*gammado221 + 
     gamma231*gammado222 + gamma332*gammado321 + gamma331*gammado322)*ginv12 \
+ (gamma133*gammado121 + gamma131*gammado123 + gamma233*gammado221 + 
     gamma231*gammado223 + gamma333*gammado321 + gamma331*gammado323)*ginv13 \
+ (gamma132*gammado122 + gamma232*gammado222 + gamma332*gammado322)*ginv22 + 
  (gamma133*gammado122 + gamma132*gammado123 + gamma233*gammado222 + 
     gamma232*gammado223 + gamma333*gammado322 + gamma332*gammado323)*ginv23 \
+ (gamma133*gammado123 + gamma233*gammado223 + gamma333*gammado323)*ginv33
;

RC23
=
(gamma121*gammado131 + gamma221*gammado231 + gamma321*gammado331)*ginv11 + 
  (gamma122*gammado131 + gamma121*gammado132 + gamma222*gammado231 + 
     gamma221*gammado232 + gamma322*gammado331 + gamma321*gammado332)*ginv12 \
+ (gamma123*gammado131 + gamma121*gammado133 + gamma223*gammado231 + 
     gamma221*gammado233 + gamma323*gammado331 + gamma321*gammado333)*ginv13 \
+ (gamma122*gammado132 + gamma222*gammado232 + gamma322*gammado332)*ginv22 + 
  (gamma123*gammado132 + gamma122*gammado133 + gamma223*gammado232 + 
     gamma222*gammado233 + gamma323*gammado332 + gamma322*gammado333)*ginv23 \
+ (gamma123*gammado133 + gamma223*gammado233 + gamma323*gammado333)*ginv33
;

RC33
=
(gamma131*gammado131 + gamma231*gammado231 + gamma331*gammado331)*ginv11 + 
  (gamma132*gammado131 + gamma131*gammado132 + gamma232*gammado231 + 
     gamma231*gammado232 + gamma332*gammado331 + gamma331*gammado332)*ginv12 \
+ (gamma133*gammado131 + gamma131*gammado133 + gamma233*gammado231 + 
     gamma231*gammado233 + gamma333*gammado331 + gamma331*gammado333)*ginv13 \
+ (gamma132*gammado132 + gamma232*gammado232 + gamma332*gammado332)*ginv22 + 
  (gamma133*gammado132 + gamma132*gammado133 + gamma233*gammado232 + 
     gamma232*gammado233 + gamma333*gammado332 + gamma332*gammado333)*ginv23 \
+ (gamma133*gammado133 + gamma233*gammado233 + gamma333*gammado333)*ginv33
;

RD11
=
-((gamma111*gammado111 + gamma211*gammado211 + gamma311*gammado311)*
     ginv11) - (gamma111*(gammado112 + gammado121) + 
     gamma211*(gammado212 + gammado221) + 
     gamma311*(gammado312 + gammado321))*ginv12 - 
  (gamma111*(gammado113 + gammado131) + 
     gamma211*(gammado213 + gammado231) + 
     gamma311*(gammado313 + gammado331))*ginv13 - 
  (gamma111*gammado122 + gamma211*gammado222 + gamma311*gammado322)*ginv22 - 
  (gamma111*(gammado123 + gammado132) + 
     gamma211*(gammado223 + gammado232) + 
     gamma311*(gammado323 + gammado332))*ginv23 - 
  (gamma111*gammado133 + gamma211*gammado233 + gamma311*gammado333)*ginv33
;

RD12
=
-((gamma112*gammado111 + gamma212*gammado211 + gamma312*gammado311)*
     ginv11) - (gamma112*(gammado112 + gammado121) + 
     gamma212*(gammado212 + gammado221) + 
     gamma312*(gammado312 + gammado321))*ginv12 - 
  (gamma112*(gammado113 + gammado131) + 
     gamma212*(gammado213 + gammado231) + 
     gamma312*(gammado313 + gammado331))*ginv13 - 
  (gamma112*gammado122 + gamma212*gammado222 + gamma312*gammado322)*ginv22 - 
  (gamma112*(gammado123 + gammado132) + 
     gamma212*(gammado223 + gammado232) + 
     gamma312*(gammado323 + gammado332))*ginv23 - 
  (gamma112*gammado133 + gamma212*gammado233 + gamma312*gammado333)*ginv33
;

RD12
=
-((gamma121*gammado111 + gamma221*gammado211 + gamma321*gammado311)*
     ginv11) - (gamma121*(gammado112 + gammado121) + 
     gamma221*(gammado212 + gammado221) + 
     gamma321*(gammado312 + gammado321))*ginv12 - 
  (gamma121*(gammado113 + gammado131) + 
     gamma221*(gammado213 + gammado231) + 
     gamma321*(gammado313 + gammado331))*ginv13 - 
  (gamma121*gammado122 + gamma221*gammado222 + gamma321*gammado322)*ginv22 - 
  (gamma121*(gammado123 + gammado132) + 
     gamma221*(gammado223 + gammado232) + 
     gamma321*(gammado323 + gammado332))*ginv23 - 
  (gamma121*gammado133 + gamma221*gammado233 + gamma321*gammado333)*ginv33
;

RD13
=
-((gamma113*gammado111 + gamma213*gammado211 + gamma313*gammado311)*
     ginv11) - (gamma113*(gammado112 + gammado121) + 
     gamma213*(gammado212 + gammado221) + 
     gamma313*(gammado312 + gammado321))*ginv12 - 
  (gamma113*(gammado113 + gammado131) + 
     gamma213*(gammado213 + gammado231) + 
     gamma313*(gammado313 + gammado331))*ginv13 - 
  (gamma113*gammado122 + gamma213*gammado222 + gamma313*gammado322)*ginv22 - 
  (gamma113*(gammado123 + gammado132) + 
     gamma213*(gammado223 + gammado232) + 
     gamma313*(gammado323 + gammado332))*ginv23 - 
  (gamma113*gammado133 + gamma213*gammado233 + gamma313*gammado333)*ginv33
;

RD13
=
-((gamma131*gammado111 + gamma231*gammado211 + gamma331*gammado311)*
     ginv11) - (gamma131*(gammado112 + gammado121) + 
     gamma231*(gammado212 + gammado221) + 
     gamma331*(gammado312 + gammado321))*ginv12 - 
  (gamma131*(gammado113 + gammado131) + 
     gamma231*(gammado213 + gammado231) + 
     gamma331*(gammado313 + gammado331))*ginv13 - 
  (gamma131*gammado122 + gamma231*gammado222 + gamma331*gammado322)*ginv22 - 
  (gamma131*(gammado123 + gammado132) + 
     gamma231*(gammado223 + gammado232) + 
     gamma331*(gammado323 + gammado332))*ginv23 - 
  (gamma131*gammado133 + gamma231*gammado233 + gamma331*gammado333)*ginv33
;

RD22
=
-((gamma122*gammado111 + gamma222*gammado211 + gamma322*gammado311)*
     ginv11) - (gamma122*(gammado112 + gammado121) + 
     gamma222*(gammado212 + gammado221) + 
     gamma322*(gammado312 + gammado321))*ginv12 - 
  (gamma122*(gammado113 + gammado131) + 
     gamma222*(gammado213 + gammado231) + 
     gamma322*(gammado313 + gammado331))*ginv13 - 
  (gamma122*gammado122 + gamma222*gammado222 + gamma322*gammado322)*ginv22 - 
  (gamma122*(gammado123 + gammado132) + 
     gamma222*(gammado223 + gammado232) + 
     gamma322*(gammado323 + gammado332))*ginv23 - 
  (gamma122*gammado133 + gamma222*gammado233 + gamma322*gammado333)*ginv33
;

RD23
=
-((gamma123*gammado111 + gamma223*gammado211 + gamma323*gammado311)*
     ginv11) - (gamma123*(gammado112 + gammado121) + 
     gamma223*(gammado212 + gammado221) + 
     gamma323*(gammado312 + gammado321))*ginv12 - 
  (gamma123*(gammado113 + gammado131) + 
     gamma223*(gammado213 + gammado231) + 
     gamma323*(gammado313 + gammado331))*ginv13 - 
  (gamma123*gammado122 + gamma223*gammado222 + gamma323*gammado322)*ginv22 - 
  (gamma123*(gammado123 + gammado132) + 
     gamma223*(gammado223 + gammado232) + 
     gamma323*(gammado323 + gammado332))*ginv23 - 
  (gamma123*gammado133 + gamma223*gammado233 + gamma323*gammado333)*ginv33
;

RD23
=
-((gamma132*gammado111 + gamma232*gammado211 + gamma332*gammado311)*
     ginv11) - (gamma132*(gammado112 + gammado121) + 
     gamma232*(gammado212 + gammado221) + 
     gamma332*(gammado312 + gammado321))*ginv12 - 
  (gamma132*(gammado113 + gammado131) + 
     gamma232*(gammado213 + gammado231) + 
     gamma332*(gammado313 + gammado331))*ginv13 - 
  (gamma132*gammado122 + gamma232*gammado222 + gamma332*gammado322)*ginv22 - 
  (gamma132*(gammado123 + gammado132) + 
     gamma232*(gammado223 + gammado232) + 
     gamma332*(gammado323 + gammado332))*ginv23 - 
  (gamma132*gammado133 + gamma232*gammado233 + gamma332*gammado333)*ginv33
;

RD33
=
-((gamma133*gammado111 + gamma233*gammado211 + gamma333*gammado311)*
     ginv11) - (gamma133*(gammado112 + gammado121) + 
     gamma233*(gammado212 + gammado221) + 
     gamma333*(gammado312 + gammado321))*ginv12 - 
  (gamma133*(gammado113 + gammado131) + 
     gamma233*(gammado213 + gammado231) + 
     gamma333*(gammado313 + gammado331))*ginv13 - 
  (gamma133*gammado122 + gamma233*gammado222 + gamma333*gammado322)*ginv22 - 
  (gamma133*(gammado123 + gammado132) + 
     gamma233*(gammado223 + gammado232) + 
     gamma333*(gammado323 + gammado332))*ginv23 - 
  (gamma133*gammado133 + gamma233*gammado233 + gamma333*gammado333)*ginv33
;

RA
=
ginv11*RA11 + ginv22*RA22 + 2.*(ginv12*RA12 + ginv13*RA13 + ginv23*RA23) + 
  ginv33*RA33
;

RB
=
ginv11*RB11 + ginv22*RB22 + 2.*(ginv12*RB12 + ginv13*RB13 + ginv23*RB23) + 
  ginv33*RB33
;

RC
=
ginv11*RC11 + ginv22*RC22 + 2.*(ginv12*RC12 + ginv13*RC13 + ginv23*RC23) + 
  ginv33*RC33
;

RD
=
ginv11*RD11 + ginv22*RD22 + 2.*(ginv12*RD12 + ginv13*RD13 + ginv23*RD23) + 
  ginv33*RD33
;

denom
=
fabs(hamrhs) + fabs(KudKud) + fabs(RA) + fabs(RB) + fabs(RC) + fabs(RD) + 
  pow2(K)
;


} else { /* if (!TermByTerm) */

denom
=
fabs(hamrhs) + fabs(KudKud) + fabs(R) + pow2(K)
;

}
/* if (TermByTerm) */


normham[ijk]
=
Cal(denom <= 0.,0.,ham[ijk]/denom)
;



/* conditional */
if (TermByTerm) {

codelKA111
=
delK111
;

codelKA112
=
delK112
;

codelKA113
=
delK113
;

codelKA122
=
delK122
;

codelKA123
=
delK123
;

codelKA133
=
delK133
;

codelKA211
=
delK211
;

codelKA212
=
delK212
;

codelKA213
=
delK213
;

codelKA222
=
delK222
;

codelKA223
=
delK223
;

codelKA233
=
delK233
;

codelKA311
=
delK311
;

codelKA312
=
delK312
;

codelKA313
=
delK313
;

codelKA322
=
delK322
;

codelKA323
=
delK323
;

codelKA333
=
delK333
;

codelKB111
=
-(gamma111*K11[ijk]) - gamma211*K12[ijk] - gamma311*K13[ijk]
;

codelKB112
=
-(gamma112*K11[ijk]) - gamma212*K12[ijk] - gamma312*K13[ijk]
;

codelKB112
=
-(gamma111*K12[ijk]) - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelKB113
=
-(gamma113*K11[ijk]) - gamma213*K12[ijk] - gamma313*K13[ijk]
;

codelKB113
=
-(gamma111*K13[ijk]) - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelKB122
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKB123
=
-(gamma113*K12[ijk]) - gamma213*K22[ijk] - gamma313*K23[ijk]
;

codelKB123
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKB133
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKB211
=
-(gamma121*K11[ijk]) - gamma221*K12[ijk] - gamma321*K13[ijk]
;

codelKB212
=
-(gamma122*K11[ijk]) - gamma222*K12[ijk] - gamma322*K13[ijk]
;

codelKB212
=
-(gamma121*K12[ijk]) - gamma221*K22[ijk] - gamma321*K23[ijk]
;

codelKB213
=
-(gamma123*K11[ijk]) - gamma223*K12[ijk] - gamma323*K13[ijk]
;

codelKB213
=
-(gamma121*K13[ijk]) - gamma221*K23[ijk] - gamma321*K33[ijk]
;

codelKB222
=
-(gamma122*K12[ijk]) - gamma222*K22[ijk] - gamma322*K23[ijk]
;

codelKB223
=
-(gamma123*K12[ijk]) - gamma223*K22[ijk] - gamma323*K23[ijk]
;

codelKB223
=
-(gamma122*K13[ijk]) - gamma222*K23[ijk] - gamma322*K33[ijk]
;

codelKB233
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKB311
=
-(gamma131*K11[ijk]) - gamma231*K12[ijk] - gamma331*K13[ijk]
;

codelKB312
=
-(gamma132*K11[ijk]) - gamma232*K12[ijk] - gamma332*K13[ijk]
;

codelKB312
=
-(gamma131*K12[ijk]) - gamma231*K22[ijk] - gamma331*K23[ijk]
;

codelKB313
=
-(gamma133*K11[ijk]) - gamma233*K12[ijk] - gamma333*K13[ijk]
;

codelKB313
=
-(gamma131*K13[ijk]) - gamma231*K23[ijk] - gamma331*K33[ijk]
;

codelKB322
=
-(gamma132*K12[ijk]) - gamma232*K22[ijk] - gamma332*K23[ijk]
;

codelKB323
=
-(gamma133*K12[ijk]) - gamma233*K22[ijk] - gamma333*K23[ijk]
;

codelKB323
=
-(gamma132*K13[ijk]) - gamma232*K23[ijk] - gamma332*K33[ijk]
;

codelKB333
=
-(gamma133*K13[ijk]) - gamma233*K23[ijk] - gamma333*K33[ijk]
;

codelKC111
=
-(gamma111*K11[ijk]) - gamma211*K12[ijk] - gamma311*K13[ijk]
;

codelKC112
=
-(gamma112*K11[ijk]) - gamma212*K12[ijk] - gamma312*K13[ijk]
;

codelKC112
=
-(gamma111*K12[ijk]) - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelKC113
=
-(gamma113*K11[ijk]) - gamma213*K12[ijk] - gamma313*K13[ijk]
;

codelKC113
=
-(gamma111*K13[ijk]) - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelKC122
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKC123
=
-(gamma113*K12[ijk]) - gamma213*K22[ijk] - gamma313*K23[ijk]
;

codelKC123
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKC133
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKC211
=
-(gamma121*K11[ijk]) - gamma221*K12[ijk] - gamma321*K13[ijk]
;

codelKC212
=
-(gamma122*K11[ijk]) - gamma222*K12[ijk] - gamma322*K13[ijk]
;

codelKC212
=
-(gamma121*K12[ijk]) - gamma221*K22[ijk] - gamma321*K23[ijk]
;

codelKC213
=
-(gamma123*K11[ijk]) - gamma223*K12[ijk] - gamma323*K13[ijk]
;

codelKC213
=
-(gamma121*K13[ijk]) - gamma221*K23[ijk] - gamma321*K33[ijk]
;

codelKC222
=
-(gamma122*K12[ijk]) - gamma222*K22[ijk] - gamma322*K23[ijk]
;

codelKC223
=
-(gamma123*K12[ijk]) - gamma223*K22[ijk] - gamma323*K23[ijk]
;

codelKC223
=
-(gamma122*K13[ijk]) - gamma222*K23[ijk] - gamma322*K33[ijk]
;

codelKC233
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKC311
=
-(gamma131*K11[ijk]) - gamma231*K12[ijk] - gamma331*K13[ijk]
;

codelKC312
=
-(gamma132*K11[ijk]) - gamma232*K12[ijk] - gamma332*K13[ijk]
;

codelKC312
=
-(gamma131*K12[ijk]) - gamma231*K22[ijk] - gamma331*K23[ijk]
;

codelKC313
=
-(gamma133*K11[ijk]) - gamma233*K12[ijk] - gamma333*K13[ijk]
;

codelKC313
=
-(gamma131*K13[ijk]) - gamma231*K23[ijk] - gamma331*K33[ijk]
;

codelKC322
=
-(gamma132*K12[ijk]) - gamma232*K22[ijk] - gamma332*K23[ijk]
;

codelKC323
=
-(gamma133*K12[ijk]) - gamma233*K22[ijk] - gamma333*K23[ijk]
;

codelKC323
=
-(gamma132*K13[ijk]) - gamma232*K23[ijk] - gamma332*K33[ijk]
;

codelKC333
=
-(gamma133*K13[ijk]) - gamma233*K23[ijk] - gamma333*K33[ijk]
;

cdKdA1
=
codelKA111*ginv11 + (codelKA112 + codelKA211)*ginv12 + 
  (codelKA113 + codelKA311)*ginv13 + codelKA212*ginv22 + 
  (codelKA213 + codelKA312)*ginv23 + codelKA313*ginv33
;

cdKdA2
=
codelKA112*ginv11 + (codelKA122 + codelKA212)*ginv12 + 
  (codelKA123 + codelKA312)*ginv13 + codelKA222*ginv22 + 
  (codelKA223 + codelKA322)*ginv23 + codelKA323*ginv33
;

cdKdA3
=
codelKA113*ginv11 + (codelKA123 + codelKA213)*ginv12 + 
  (codelKA133 + codelKA313)*ginv13 + codelKA223*ginv22 + 
  (codelKA233 + codelKA323)*ginv23 + codelKA333*ginv33
;

cdKdB1
=
codelKB111*ginv11 + (codelKB112 + codelKB211)*ginv12 + 
  (codelKB113 + codelKB311)*ginv13 + codelKB212*ginv22 + 
  (codelKB213 + codelKB312)*ginv23 + codelKB313*ginv33
;

cdKdB2
=
codelKB112*ginv11 + (codelKB122 + codelKB212)*ginv12 + 
  (codelKB123 + codelKB312)*ginv13 + codelKB222*ginv22 + 
  (codelKB223 + codelKB322)*ginv23 + codelKB323*ginv33
;

cdKdB3
=
codelKB113*ginv11 + (codelKB123 + codelKB213)*ginv12 + 
  (codelKB133 + codelKB313)*ginv13 + codelKB223*ginv22 + 
  (codelKB233 + codelKB323)*ginv23 + codelKB333*ginv33
;

cdKdC1
=
codelKC111*ginv11 + (codelKC112 + codelKC211)*ginv12 + 
  (codelKC113 + codelKC311)*ginv13 + codelKC212*ginv22 + 
  (codelKC213 + codelKC312)*ginv23 + codelKC313*ginv33
;

cdKdC2
=
codelKC112*ginv11 + (codelKC122 + codelKC212)*ginv12 + 
  (codelKC123 + codelKC312)*ginv13 + codelKC222*ginv22 + 
  (codelKC223 + codelKC322)*ginv23 + codelKC323*ginv33
;

cdKdC3
=
codelKC113*ginv11 + (codelKC123 + codelKC213)*ginv12 + 
  (codelKC133 + codelKC313)*ginv13 + codelKC223*ginv22 + 
  (codelKC233 + codelKC323)*ginv23 + codelKC333*ginv33
;

codelTrKA1
=
codelKA111*ginv11 + codelKA122*ginv22 + 
  2.*(codelKA112*ginv12 + codelKA113*ginv13 + codelKA123*ginv23) + 
  codelKA133*ginv33
;

codelTrKA2
=
codelKA211*ginv11 + codelKA222*ginv22 + 
  2.*(codelKA212*ginv12 + codelKA213*ginv13 + codelKA223*ginv23) + 
  codelKA233*ginv33
;

codelTrKA3
=
codelKA311*ginv11 + codelKA322*ginv22 + 
  2.*(codelKA312*ginv12 + codelKA313*ginv13 + codelKA323*ginv23) + 
  codelKA333*ginv33
;

codelTrKB1
=
codelKB111*ginv11 + codelKB122*ginv22 + 
  2.*(codelKB112*ginv12 + codelKB113*ginv13 + codelKB123*ginv23) + 
  codelKB133*ginv33
;

codelTrKB2
=
codelKB211*ginv11 + codelKB222*ginv22 + 
  2.*(codelKB212*ginv12 + codelKB213*ginv13 + codelKB223*ginv23) + 
  codelKB233*ginv33
;

codelTrKB3
=
codelKB311*ginv11 + codelKB322*ginv22 + 
  2.*(codelKB312*ginv12 + codelKB313*ginv13 + codelKB323*ginv23) + 
  codelKB333*ginv33
;

codelTrKC1
=
codelKC111*ginv11 + codelKC122*ginv22 + 
  2.*(codelKC112*ginv12 + codelKC113*ginv13 + codelKC123*ginv23) + 
  codelKC133*ginv33
;

codelTrKC2
=
codelKC211*ginv11 + codelKC222*ginv22 + 
  2.*(codelKC212*ginv12 + codelKC213*ginv13 + codelKC223*ginv23) + 
  codelKC233*ginv33
;

codelTrKC3
=
codelKC311*ginv11 + codelKC322*ginv22 + 
  2.*(codelKC312*ginv12 + codelKC313*ginv13 + codelKC323*ginv23) + 
  codelKC333*ginv33
;

denom1
=
fabs(cdKdA1) + fabs(cdKdB1) + fabs(cdKdC1) + fabs(codelTrKA1) + 
  fabs(codelTrKB1) + fabs(codelTrKC1) + fabs(momrhs1)
;

denom2
=
fabs(cdKdA2) + fabs(cdKdB2) + fabs(cdKdC2) + fabs(codelTrKA2) + 
  fabs(codelTrKB2) + fabs(codelTrKC2) + fabs(momrhs2)
;

denom3
=
fabs(cdKdA3) + fabs(cdKdB3) + fabs(cdKdC3) + fabs(codelTrKA3) + 
  fabs(codelTrKB3) + fabs(codelTrKC3) + fabs(momrhs3)
;


} else { /* if (!TermByTerm) */

denom1
=
fabs(cdKudd111 + cdKudd212 + cdKudd313) + 
  fabs(codelK111*ginv11 + codelK122*ginv22 + 
    2.*(codelK112*ginv12 + codelK113*ginv13 + codelK123*ginv23) + 
    codelK133*ginv33) + fabs(momrhs1)
;

denom2
=
fabs(cdKudd112 + cdKudd222 + cdKudd323) + 
  fabs(codelK211*ginv11 + codelK222*ginv22 + 
    2.*(codelK212*ginv12 + codelK213*ginv13 + codelK223*ginv23) + 
    codelK233*ginv33) + fabs(momrhs2)
;

denom3
=
fabs(cdKudd113 + cdKudd223 + cdKudd333) + 
  fabs(codelK311*ginv11 + codelK322*ginv22 + 
    2.*(codelK312*ginv12 + codelK313*ginv13 + codelK323*ginv23) + 
    codelK333*ginv33) + fabs(momrhs3)
;

}
/* if (TermByTerm) */


denom
=
denom1 + denom2 + denom3
;

normmom1[ijk]
=
Cal(denom <= 0.,0.,mom1[ijk]/denom)
;

normmom2[ijk]
=
Cal(denom <= 0.,0.,mom2[ijk]/denom)
;

normmom3[ijk]
=
Cal(denom <= 0.,0.,mom3[ijk]/denom)
;

}
/* if (TermByTerm) */


} /* end of points */
} /* end of boxes */


}  /* end of function */

/* ADMconstraints.c */
/* nvars = 107, n* = 2000,  n/ = 36,  n+ = 1954, n = 3990, O = 1 */
