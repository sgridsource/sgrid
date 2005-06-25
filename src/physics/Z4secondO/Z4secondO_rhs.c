/* Z4secondO_rhs.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 25.6.2005 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "Z4secondO.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void Z4secondO_rhs(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)
{
tGrid *grid = ucur->grid;
int bi;

int addlinear = (dt != 0.0l);
/* int usepsi = 1; */
int useDD               = Getv("Z4secondO_useDD", "yes");
double RtoRminusHfactor = Getd("Z4secondO_RtoRminusHfactor");
int nonconstantlapse    =!Getv("Z4secondO_lapse", "constant");
int oploglapse          = Getv("Z4secondO_lapse", "1+log");
int withshift           = Getv("Z4secondO_lapse", "withshift");
int harmoniclapse       = Getv("Z4secondO_lapse", "harmonic");
int subtractK0          = Getv("Z4secondO_subtractK0", "yes");
double gamma0factor    = Getv("Z4secondO_shift", "gamma0");
double gamma2factor    = Getv("Z4secondO_shift", "gamma2");
double shift_st        = Getd("Z4secondO_shift_stop_time");
double evolveshift     = Cal( (grid->time>=shift_st)*(shift_st>=0), 0.0, 1.0);
double lapsepsipower   = Getd("Z4secondO_lapsepsipower");
double lapseharmonicf  = Getd("Z4secondO_lapseharmonicf");
double lapseharmonicm  = Getd("Z4secondO_lapseharmonicm");
double shiftpsipower   = Getd("Z4secondO_shiftpsipower");
double shiftalphapower = Getd("Z4secondO_shiftalphapower");
double shiftgammacoeff = Getd("Z4secondO_shiftgammacoeff");
double shiftdriver     = Getd("Z4secondO_shiftdriver");

for(bi = 0; bi < grid->nboxes; bi++)
{
tBox *box = grid->box[bi];
int ijk;



double *K0 = vlldataptr(K_initial, box, 0);
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
double *g11 = vlldataptr(ucur, box, 0);
double *g12 = vlldataptr(ucur, box, 1);
double *g13 = vlldataptr(ucur, box, 2);
double *g22 = vlldataptr(ucur, box, 3);
double *g23 = vlldataptr(ucur, box, 4);
double *g33 = vlldataptr(ucur, box, 5);
double *K11 = vlldataptr(ucur, box, 6);
double *K12 = vlldataptr(ucur, box, 7);
double *K13 = vlldataptr(ucur, box, 8);
double *K22 = vlldataptr(ucur, box, 9);
double *K23 = vlldataptr(ucur, box, 10);
double *K33 = vlldataptr(ucur, box, 11);
double *Theta = vlldataptr(ucur, box, 12);
double *Z1 = vlldataptr(ucur, box, 13);
double *Z2 = vlldataptr(ucur, box, 14);
double *Z3 = vlldataptr(ucur, box, 15);
double *alpha = vlldataptr(ucur, box, 16);
double *beta1 = vlldataptr(ucur, box, 17);
double *beta2 = vlldataptr(ucur, box, 18);
double *beta3 = vlldataptr(ucur, box, 19);
double *B1 = vlldataptr(ucur, box, 20);
double *B2 = vlldataptr(ucur, box, 21);
double *B3 = vlldataptr(ucur, box, 22);
double *ng11 = vlldataptr(unew, box, 0);
double *ng12 = vlldataptr(unew, box, 1);
double *ng13 = vlldataptr(unew, box, 2);
double *ng22 = vlldataptr(unew, box, 3);
double *ng23 = vlldataptr(unew, box, 4);
double *ng33 = vlldataptr(unew, box, 5);
double *nK11 = vlldataptr(unew, box, 6);
double *nK12 = vlldataptr(unew, box, 7);
double *nK13 = vlldataptr(unew, box, 8);
double *nK22 = vlldataptr(unew, box, 9);
double *nK23 = vlldataptr(unew, box, 10);
double *nK33 = vlldataptr(unew, box, 11);
double *nTheta = vlldataptr(unew, box, 12);
double *nZ1 = vlldataptr(unew, box, 13);
double *nZ2 = vlldataptr(unew, box, 14);
double *nZ3 = vlldataptr(unew, box, 15);
double *nalpha = vlldataptr(unew, box, 16);
double *nbeta1 = vlldataptr(unew, box, 17);
double *nbeta2 = vlldataptr(unew, box, 18);
double *nbeta3 = vlldataptr(unew, box, 19);
double *nB1 = vlldataptr(unew, box, 20);
double *nB2 = vlldataptr(unew, box, 21);
double *nB3 = vlldataptr(unew, box, 22);
double *pg11 = vlldataptr(upre, box, 0);
double *pg12 = vlldataptr(upre, box, 1);
double *pg13 = vlldataptr(upre, box, 2);
double *pg22 = vlldataptr(upre, box, 3);
double *pg23 = vlldataptr(upre, box, 4);
double *pg33 = vlldataptr(upre, box, 5);
double *pK11 = vlldataptr(upre, box, 6);
double *pK12 = vlldataptr(upre, box, 7);
double *pK13 = vlldataptr(upre, box, 8);
double *pK22 = vlldataptr(upre, box, 9);
double *pK23 = vlldataptr(upre, box, 10);
double *pK33 = vlldataptr(upre, box, 11);
double *pTheta = vlldataptr(upre, box, 12);
double *pZ1 = vlldataptr(upre, box, 13);
double *pZ2 = vlldataptr(upre, box, 14);
double *pZ3 = vlldataptr(upre, box, 15);
double *palpha = vlldataptr(upre, box, 16);
double *pbeta1 = vlldataptr(upre, box, 17);
double *pbeta2 = vlldataptr(upre, box, 18);
double *pbeta3 = vlldataptr(upre, box, 19);
double *pB1 = vlldataptr(upre, box, 20);
double *pB2 = vlldataptr(upre, box, 21);
double *pB3 = vlldataptr(upre, box, 22);
int index_g11 = (ucur)->index[0];
int index_g12 = (ucur)->index[1];
int index_g13 = (ucur)->index[2];
int index_g22 = (ucur)->index[3];
int index_g23 = (ucur)->index[4];
int index_g33 = (ucur)->index[5];
int index_K11 = (ucur)->index[6];
int index_K12 = (ucur)->index[7];
int index_K13 = (ucur)->index[8];
int index_K22 = (ucur)->index[9];
int index_K23 = (ucur)->index[10];
int index_K33 = (ucur)->index[11];
int index_Theta = (ucur)->index[12];
int index_Z1 = (ucur)->index[13];
int index_Z2 = (ucur)->index[14];
int index_Z3 = (ucur)->index[15];
int index_alpha = (ucur)->index[16];
int index_beta1 = (ucur)->index[17];
int index_beta2 = (ucur)->index[18];
int index_beta3 = (ucur)->index[19];
int index_B1 = (ucur)->index[20];
int index_B2 = (ucur)->index[21];
int index_B3 = (ucur)->index[22];
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
int index_Z4secondO_dThetax = Ind("Z4secondO_dThetax");
double *dTheta1 = box->v[index_Z4secondO_dThetax + 0];
double *dTheta2 = box->v[index_Z4secondO_dThetax + 1];
double *dTheta3 = box->v[index_Z4secondO_dThetax + 2];
int index_Z4secondO_dZxx = Ind("Z4secondO_dZxx");
double *dZ11 = box->v[index_Z4secondO_dZxx + 0];
double *dZ12 = box->v[index_Z4secondO_dZxx + 1];
double *dZ13 = box->v[index_Z4secondO_dZxx + 2];
double *dZ21 = box->v[index_Z4secondO_dZxx + 3];
double *dZ22 = box->v[index_Z4secondO_dZxx + 4];
double *dZ23 = box->v[index_Z4secondO_dZxx + 5];
double *dZ31 = box->v[index_Z4secondO_dZxx + 6];
double *dZ32 = box->v[index_Z4secondO_dZxx + 7];
double *dZ33 = box->v[index_Z4secondO_dZxx + 8];
int index_Z4secondO_dalpx = Ind("Z4secondO_dalpx");
double *dalp1 = box->v[index_Z4secondO_dalpx + 0];
double *dalp2 = box->v[index_Z4secondO_dalpx + 1];
double *dalp3 = box->v[index_Z4secondO_dalpx + 2];
int index_Z4secondO_ddalpxx = Ind("Z4secondO_ddalpxx");
double *ddalp11 = box->v[index_Z4secondO_ddalpxx + 0];
double *ddalp12 = box->v[index_Z4secondO_ddalpxx + 1];
double *ddalp13 = box->v[index_Z4secondO_ddalpxx + 2];
double *ddalp22 = box->v[index_Z4secondO_ddalpxx + 3];
double *ddalp23 = box->v[index_Z4secondO_ddalpxx + 4];
double *ddalp33 = box->v[index_Z4secondO_ddalpxx + 5];
int index_Z4secondO_dbetaxx = Ind("Z4secondO_dbetaxx");
double *dbeta11 = box->v[index_Z4secondO_dbetaxx + 0];
double *dbeta12 = box->v[index_Z4secondO_dbetaxx + 1];
double *dbeta13 = box->v[index_Z4secondO_dbetaxx + 2];
double *dbeta21 = box->v[index_Z4secondO_dbetaxx + 3];
double *dbeta22 = box->v[index_Z4secondO_dbetaxx + 4];
double *dbeta23 = box->v[index_Z4secondO_dbetaxx + 5];
double *dbeta31 = box->v[index_Z4secondO_dbetaxx + 6];
double *dbeta32 = box->v[index_Z4secondO_dbetaxx + 7];
double *dbeta33 = box->v[index_Z4secondO_dbetaxx + 8];
int index_Z4secondO_ddbetaxxx = Ind("Z4secondO_ddbetaxxx");
double *ddbeta111 = box->v[index_Z4secondO_ddbetaxxx + 0];
double *ddbeta112 = box->v[index_Z4secondO_ddbetaxxx + 1];
double *ddbeta113 = box->v[index_Z4secondO_ddbetaxxx + 2];
double *ddbeta122 = box->v[index_Z4secondO_ddbetaxxx + 3];
double *ddbeta123 = box->v[index_Z4secondO_ddbetaxxx + 4];
double *ddbeta133 = box->v[index_Z4secondO_ddbetaxxx + 5];
double *ddbeta211 = box->v[index_Z4secondO_ddbetaxxx + 6];
double *ddbeta212 = box->v[index_Z4secondO_ddbetaxxx + 7];
double *ddbeta213 = box->v[index_Z4secondO_ddbetaxxx + 8];
double *ddbeta222 = box->v[index_Z4secondO_ddbetaxxx + 9];
double *ddbeta223 = box->v[index_Z4secondO_ddbetaxxx + 10];
double *ddbeta233 = box->v[index_Z4secondO_ddbetaxxx + 11];
double *ddbeta311 = box->v[index_Z4secondO_ddbetaxxx + 12];
double *ddbeta312 = box->v[index_Z4secondO_ddbetaxxx + 13];
double *ddbeta313 = box->v[index_Z4secondO_ddbetaxxx + 14];
double *ddbeta322 = box->v[index_Z4secondO_ddbetaxxx + 15];
double *ddbeta323 = box->v[index_Z4secondO_ddbetaxxx + 16];
double *ddbeta333 = box->v[index_Z4secondO_ddbetaxxx + 17];


double betadalp;
double betadg11;
double betadg12;
double betadg13;
double betadg21;
double betadg22;
double betadg23;
double betadg31;
double betadg32;
double betadg33;
double betadK11;
double betadK12;
double betadK13;
double betadK21;
double betadK22;
double betadK23;
double betadK31;
double betadK32;
double betadK33;
double betadTheta;
double betadZ1;
double betadZ2;
double betadZ3;
double cddalp11;
double cddalp12;
double cddalp13;
double cddalp22;
double cddalp23;
double cddalp33;
double cdK111;
double cdK112;
double cdK113;
double cdK121;
double cdK122;
double cdK123;
double cdK131;
double cdK132;
double cdK133;
double cdK211;
double cdK212;
double cdK213;
double cdK221;
double cdK222;
double cdK223;
double cdK231;
double cdK232;
double cdK233;
double cdK311;
double cdK312;
double cdK313;
double cdK321;
double cdK322;
double cdK323;
double cdK331;
double cdK332;
double cdK333;
double cdZ11;
double cdZ12;
double cdZ13;
double cdZ21;
double cdZ22;
double cdZ23;
double cdZ31;
double cdZ32;
double cdZ33;
double detginv;
double dginv111;
double dginv112;
double dginv113;
double dginv121;
double dginv122;
double dginv123;
double dginv131;
double dginv132;
double dginv133;
double dginv221;
double dginv222;
double dginv223;
double dginv231;
double dginv232;
double dginv233;
double dginv331;
double dginv332;
double dginv333;
double dK1;
double dK2;
double dK3;
double gamma111;
double gamma112;
double gamma113;
double gamma122;
double gamma123;
double gamma133;
double gamma211;
double gamma212;
double gamma213;
double gamma222;
double gamma223;
double gamma233;
double gamma311;
double gamma312;
double gamma313;
double gamma322;
double gamma323;
double gamma333;
double gammado111;
double gammado112;
double gammado113;
double gammado122;
double gammado123;
double gammado133;
double gammado211;
double gammado212;
double gammado213;
double gammado222;
double gammado223;
double gammado233;
double gammado311;
double gammado312;
double gammado313;
double gammado322;
double gammado323;
double gammado333;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double hamil;
double KK11;
double KK12;
double KK13;
double KK22;
double KK23;
double KK33;
double liealpha;
double lieg11;
double lieg12;
double lieg13;
double lieg22;
double lieg23;
double lieg33;
double lieK11;
double lieK12;
double lieK13;
double lieK22;
double lieK23;
double lieK33;
double lieTheta;
double lieZ1;
double lieZ2;
double lieZ3;
double R;
double R11;
double R12;
double R13;
double R22;
double R23;
double R33;
double ralpha;
double ralpha0;
double rB1;
double rB2;
double rB3;
double rbeta1;
double rbeta2;
double rbeta3;
double rg11;
double rg12;
double rg13;
double rg22;
double rg23;
double rg33;
double rK11;
double rK12;
double rK13;
double rK22;
double rK23;
double rK33;
double rTheta;
double rZ1;
double rZ2;
double rZ3;
double trK;
double trKK;



/* Jetzt geht's los! */


/* conditional */
if (useDD) {


allDerivsOf_Sab(box, index_g11,                     Ind("ADMvars_dgxxx"), Ind("ADMvars_ddgxxxx")); 


FirstDerivsOf_Sab(box, index_K11, Ind("ADMvars_dKxxx")); 


FirstDerivsOf_Sa(box, index_Z1, Ind("Z4secondO_dZxx")); 


FirstDerivsOf_S(box, index_Theta, Ind("Z4secondO_dThetax")); 


allDerivsOf_S(box, index_alpha,                     Ind("Z4secondO_dalpx"), Ind("Z4secondO_ddalpxx")); 


allDerivsOf_Sa(box, index_beta1,                     Ind("Z4secondO_dbetaxx"), Ind("Z4secondO_ddbetaxxx")); 


} else { /* if (!useDD) */


FirstAndSecondDerivsOf_Sab(box, index_g11,                     Ind("ADMvars_dgxxx"), Ind("ADMvars_ddgxxxx")); 


FirstDerivsOf_Sab(box, index_K11, Ind("ADMvars_dKxxx")); 


FirstDerivsOf_Sa(box, index_Z1, Ind("Z4secondO_dZxx")); 


FirstDerivsOf_S(box, index_Theta, Ind("Z4secondO_dThetax")); 


FirstAndSecondDerivsOf_S(box, index_alpha,                     Ind("Z4secondO_dalpx"), Ind("Z4secondO_ddalpxx")); 


FirstAndSecondDerivsOf_Sa(box, index_beta1,                     Ind("Z4secondO_dbetaxx"), Ind("Z4secondO_ddbetaxxx")); 

}
/* if (useDD) */



forallpoints(box, ijk) { 

betadg11
=
beta1[ijk]*dg111[ijk] + beta2[ijk]*dg112[ijk] + beta3[ijk]*dg113[ijk]
;

betadg12
=
beta1[ijk]*dg121[ijk] + beta2[ijk]*dg122[ijk] + beta3[ijk]*dg123[ijk]
;

betadg13
=
beta1[ijk]*dg131[ijk] + beta2[ijk]*dg132[ijk] + beta3[ijk]*dg133[ijk]
;

betadg21
=
beta1[ijk]*dg121[ijk] + beta2[ijk]*dg122[ijk] + beta3[ijk]*dg123[ijk]
;

betadg22
=
beta1[ijk]*dg221[ijk] + beta2[ijk]*dg222[ijk] + beta3[ijk]*dg223[ijk]
;

betadg23
=
beta1[ijk]*dg231[ijk] + beta2[ijk]*dg232[ijk] + beta3[ijk]*dg233[ijk]
;

betadg31
=
beta1[ijk]*dg131[ijk] + beta2[ijk]*dg132[ijk] + beta3[ijk]*dg133[ijk]
;

betadg32
=
beta1[ijk]*dg231[ijk] + beta2[ijk]*dg232[ijk] + beta3[ijk]*dg233[ijk]
;

betadg33
=
beta1[ijk]*dg331[ijk] + beta2[ijk]*dg332[ijk] + beta3[ijk]*dg333[ijk]
;

betadK11
=
beta1[ijk]*dK111[ijk] + beta2[ijk]*dK112[ijk] + beta3[ijk]*dK113[ijk]
;

betadK12
=
beta1[ijk]*dK121[ijk] + beta2[ijk]*dK122[ijk] + beta3[ijk]*dK123[ijk]
;

betadK13
=
beta1[ijk]*dK131[ijk] + beta2[ijk]*dK132[ijk] + beta3[ijk]*dK133[ijk]
;

betadK21
=
beta1[ijk]*dK121[ijk] + beta2[ijk]*dK122[ijk] + beta3[ijk]*dK123[ijk]
;

betadK22
=
beta1[ijk]*dK221[ijk] + beta2[ijk]*dK222[ijk] + beta3[ijk]*dK223[ijk]
;

betadK23
=
beta1[ijk]*dK231[ijk] + beta2[ijk]*dK232[ijk] + beta3[ijk]*dK233[ijk]
;

betadK31
=
beta1[ijk]*dK131[ijk] + beta2[ijk]*dK132[ijk] + beta3[ijk]*dK133[ijk]
;

betadK32
=
beta1[ijk]*dK231[ijk] + beta2[ijk]*dK232[ijk] + beta3[ijk]*dK233[ijk]
;

betadK33
=
beta1[ijk]*dK331[ijk] + beta2[ijk]*dK332[ijk] + beta3[ijk]*dK333[ijk]
;

betadTheta
=
beta1[ijk]*dTheta1[ijk] + beta2[ijk]*dTheta2[ijk] + beta3[ijk]*dTheta3[ijk]
;

betadZ1
=
beta1[ijk]*dZ11[ijk] + beta2[ijk]*dZ12[ijk] + beta3[ijk]*dZ13[ijk]
;

betadZ2
=
beta1[ijk]*dZ21[ijk] + beta2[ijk]*dZ22[ijk] + beta3[ijk]*dZ23[ijk]
;

betadZ3
=
beta1[ijk]*dZ31[ijk] + beta2[ijk]*dZ32[ijk] + beta3[ijk]*dZ33[ijk]
;

betadalp
=
beta1[ijk]*dalp1[ijk] + beta2[ijk]*dalp2[ijk] + beta3[ijk]*dalp3[ijk]
;

detginv
=
1/(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
    g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
    g11[ijk]*pow2(g23[ijk]))
;

ginv11
=
detginv*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

ginv12
=
detginv*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

ginv13
=
detginv*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

ginv22
=
detginv*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

ginv23
=
detginv*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

ginv33
=
detginv*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
;

gammado111
=
0.5*dg111[ijk]
;

gammado112
=
0.5*dg112[ijk]
;

gammado113
=
0.5*dg113[ijk]
;

gammado122
=
dg122[ijk] - 0.5*dg221[ijk]
;

gammado123
=
0.5*(dg123[ijk] + dg132[ijk] - dg231[ijk])
;

gammado133
=
dg133[ijk] - 0.5*dg331[ijk]
;

gammado211
=
-0.5*dg112[ijk] + dg121[ijk]
;

gammado212
=
0.5*dg221[ijk]
;

gammado213
=
0.5*(dg123[ijk] - dg132[ijk] + dg231[ijk])
;

gammado222
=
0.5*dg222[ijk]
;

gammado223
=
0.5*dg223[ijk]
;

gammado233
=
dg233[ijk] - 0.5*dg332[ijk]
;

gammado311
=
-0.5*dg113[ijk] + dg131[ijk]
;

gammado312
=
0.5*(-dg123[ijk] + dg132[ijk] + dg231[ijk])
;

gammado313
=
0.5*dg331[ijk]
;

gammado322
=
-0.5*dg223[ijk] + dg232[ijk]
;

gammado323
=
0.5*dg332[ijk]
;

gammado333
=
0.5*dg333[ijk]
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

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
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

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
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

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

R11
=
(gamma112*gammado111 - gamma111*gammado112 + gamma212*gammado211 - 
     gamma211*gammado212 + gamma312*gammado311 - gamma311*gammado312)*ginv12 \
+ (gamma113*gammado111 - gamma111*gammado113 + gamma213*gammado211 - 
     gamma211*gammado213 + gamma313*gammado311 - gamma311*gammado313)*ginv13 \
+ ginv22*(gamma112*gammado112 - gamma111*gammado122 + gamma212*gammado212 - 
     gamma211*gammado222 + gamma312*gammado312 - gamma311*gammado322 + 
     ddg1212[ijk] + 0.5*(-ddg1122[ijk] - ddg2211[ijk])) + 
  ginv23*(gamma113*gammado112 + gamma112*gammado113 + gamma213*gammado212 + 
     gamma212*gammado213 + gamma313*gammado312 + gamma312*gammado313 - 
     2.*(gamma111*gammado123 + gamma211*gammado223 + gamma311*gammado323) - 
     ddg1123[ijk] + ddg1213[ijk] + ddg1312[ijk] - ddg2311[ijk]) + 
  ginv33*(gamma113*gammado113 - gamma111*gammado133 + gamma213*gammado213 - 
     gamma211*gammado233 + gamma313*gammado313 - gamma311*gammado333 + 
     ddg1313[ijk] + 0.5*(-ddg1133[ijk] - ddg3311[ijk]))
;

R12
=
(gamma122*gammado112 - gamma112*gammado122 + gamma222*gammado212 - 
     gamma212*gammado222 + gamma322*gammado312 - gamma312*gammado322)*ginv22 \
+ ginv12*(gamma122*gammado111 - gamma112*gammado112 + gamma222*gammado211 - 
     gamma212*gammado212 + gamma322*gammado311 - gamma312*gammado312 - 
     ddg1212[ijk] + 0.5*(ddg1122[ijk] + ddg2211[ijk])) + 
  ginv13*(gamma123*gammado111 - gamma112*gammado113 + gamma223*gammado211 - 
     gamma212*gammado213 + gamma323*gammado311 - gamma312*gammado313 - 
     0.5*ddg1213[ijk] + 0.5*(ddg1123[ijk] - ddg1312[ijk] + ddg2311[ijk])) + 
  ginv23*(gamma123*gammado112 + gamma122*gammado113 + gamma223*gammado212 + 
     gamma222*gammado213 + gamma323*gammado312 + gamma322*gammado313 - 
     2.*(gamma112*gammado123 + gamma212*gammado223 + gamma312*gammado323) + 
     0.5*(-ddg1223[ijk] + ddg1322[ijk] + ddg2213[ijk] - ddg2312[ijk])) + 
  ginv33*(gamma123*gammado113 - gamma112*gammado133 + gamma223*gammado213 - 
     gamma212*gammado233 + gamma323*gammado313 - gamma312*gammado333 + 
     0.5*(ddg1323[ijk] + ddg2313[ijk]) + 0.5*(-ddg1233[ijk] - ddg3312[ijk]))
;

R12
=
(-(gamma112*gammado111) + gamma111*gammado112 - gamma212*gammado211 + 
     gamma211*gammado212 - gamma312*gammado311 + gamma311*gammado312)*ginv11 \
+ ginv12*(-(gamma112*gammado112) + gamma111*gammado122 - 
     gamma212*gammado212 + gamma211*gammado222 - gamma312*gammado312 + 
     gamma311*gammado322 - ddg1212[ijk] + 0.5*(ddg1122[ijk] + ddg2211[ijk])) \
+ ginv13*(gamma113*gammado112 + gamma111*gammado123 + gamma213*gammado212 + 
     gamma211*gammado223 + gamma313*gammado312 - 
     2.*(gamma112*gammado113 + gamma212*gammado213 + gamma312*gammado313) + 
     gamma311*gammado323 + 0.5*
      (ddg1123[ijk] - ddg1213[ijk] - ddg1312[ijk] + ddg2311[ijk])) + 
  ginv23*(gamma113*gammado122 - gamma112*gammado123 + gamma213*gammado222 - 
     gamma212*gammado223 + gamma313*gammado322 - gamma312*gammado323 + 
     0.5*(ddg1322[ijk] + ddg2213[ijk]) + 0.5*(-ddg1223[ijk] - ddg2312[ijk])) \
+ ginv33*(gamma113*gammado123 - gamma112*gammado133 + gamma213*gammado223 - 
     gamma212*gammado233 + gamma313*gammado323 - gamma312*gammado333 + 
     0.5*(ddg1323[ijk] + ddg2313[ijk]) + 0.5*(-ddg1233[ijk] - ddg3312[ijk]))
;

R13
=
(-(gamma113*gammado111) + gamma111*gammado113 - gamma213*gammado211 + 
     gamma211*gammado213 - gamma313*gammado311 + gamma311*gammado313)*ginv11 \
+ ginv12*(gamma112*gammado113 + gamma111*gammado123 + gamma212*gammado213 + 
     gamma211*gammado223 - 2.*(gamma113*gammado112 + gamma213*gammado212 + 
        gamma313*gammado312) + gamma312*gammado313 + gamma311*gammado323 + 
     0.5*(ddg1123[ijk] - ddg1213[ijk] - ddg1312[ijk] + ddg2311[ijk])) + 
  ginv22*(-(gamma113*gammado122) + gamma112*gammado123 - 
     gamma213*gammado222 + gamma212*gammado223 - gamma313*gammado322 + 
     gamma312*gammado323 - 0.5*ddg1322[ijk] + 
     0.5*(ddg1223[ijk] - ddg2213[ijk] + ddg2312[ijk])) + 
  ginv13*(-(gamma113*gammado113) + gamma111*gammado133 - 
     gamma213*gammado213 + gamma211*gammado233 - gamma313*gammado313 + 
     gamma311*gammado333 - ddg1313[ijk] + 0.5*(ddg1133[ijk] + ddg3311[ijk])) \
+ ginv23*(-(gamma113*gammado123) + gamma112*gammado133 - 
     gamma213*gammado223 + gamma212*gammado233 - gamma313*gammado323 + 
     gamma312*gammado333 - 0.5*ddg1323[ijk] + 
     0.5*(ddg1233[ijk] - ddg2313[ijk] + ddg3312[ijk]))
;

R13
=
(gamma133*gammado113 - gamma113*gammado133 + gamma233*gammado213 - 
     gamma213*gammado233 + gamma333*gammado313 - gamma313*gammado333)*ginv33 \
+ ginv12*(gamma123*gammado111 - gamma113*gammado112 + gamma223*gammado211 - 
     gamma213*gammado212 + gamma323*gammado311 - gamma313*gammado312 - 
     0.5*ddg1213[ijk] + 0.5*(ddg1123[ijk] - ddg1312[ijk] + ddg2311[ijk])) + 
  ginv22*(gamma123*gammado112 - gamma113*gammado122 + gamma223*gammado212 - 
     gamma213*gammado222 + gamma323*gammado312 - gamma313*gammado322 - 
     0.5*ddg1322[ijk] + 0.5*(ddg1223[ijk] - ddg2213[ijk] + ddg2312[ijk])) + 
  ginv13*(gamma133*gammado111 - gamma113*gammado113 + gamma233*gammado211 - 
     gamma213*gammado213 + gamma333*gammado311 - gamma313*gammado313 - 
     ddg1313[ijk] + 0.5*(ddg1133[ijk] + ddg3311[ijk])) + 
  ginv23*(gamma133*gammado112 + gamma123*gammado113 + gamma233*gammado212 + 
     gamma223*gammado213 + gamma333*gammado312 + gamma323*gammado313 - 
     2.*(gamma113*gammado123 + gamma213*gammado223 + gamma313*gammado323) + 
     0.5*(ddg1233[ijk] - ddg1323[ijk] - ddg2313[ijk] + ddg3312[ijk]))
;

R22
=
(-(gamma122*gammado112) + gamma112*gammado122 - gamma222*gammado212 + 
     gamma212*gammado222 - gamma322*gammado312 + gamma312*gammado322)*ginv12 \
+ (gamma123*gammado122 - gamma122*gammado123 + gamma223*gammado222 - 
     gamma222*gammado223 + gamma323*gammado322 - gamma322*gammado323)*ginv23 \
+ ginv11*(-(gamma122*gammado111) + gamma112*gammado112 - 
     gamma222*gammado211 + gamma212*gammado212 - gamma322*gammado311 + 
     gamma312*gammado312 + ddg1212[ijk] + 0.5*(-ddg1122[ijk] - ddg2211[ijk])\
) + ginv13*(gamma123*gammado112 + gamma112*gammado123 + 
     gamma223*gammado212 + gamma212*gammado223 + gamma323*gammado312 - 
     2.*(gamma122*gammado113 + gamma222*gammado213 + gamma322*gammado313) + 
     gamma312*gammado323 + ddg1223[ijk] - ddg1322[ijk] - ddg2213[ijk] + 
     ddg2312[ijk]) + ginv33*(gamma123*gammado123 - gamma122*gammado133 + 
     gamma223*gammado223 - gamma222*gammado233 + gamma323*gammado323 - 
     gamma322*gammado333 + ddg2323[ijk] + 0.5*(-ddg2233[ijk] - ddg3322[ijk]))
;

R23
=
(-(gamma123*gammado122) + gamma122*gammado123 - gamma223*gammado222 + 
     gamma222*gammado223 - gamma323*gammado322 + gamma322*gammado323)*ginv22 \
+ ginv11*(-(gamma123*gammado111) + gamma112*gammado113 - 
     gamma223*gammado211 + gamma212*gammado213 - gamma323*gammado311 + 
     gamma312*gammado313 + 0.5*(ddg1213[ijk] + ddg1312[ijk]) + 
     0.5*(-ddg1123[ijk] - ddg2311[ijk])) + 
  ginv12*(gamma122*gammado113 + gamma112*gammado123 + gamma222*gammado213 + 
     gamma212*gammado223 - 2.*(gamma123*gammado112 + gamma223*gammado212 + 
        gamma323*gammado312) + gamma322*gammado313 + gamma312*gammado323 + 
     0.5*(-ddg1223[ijk] + ddg1322[ijk] + ddg2213[ijk] - ddg2312[ijk])) + 
  ginv13*(-(gamma123*gammado113) + gamma112*gammado133 - 
     gamma223*gammado213 + gamma212*gammado233 - gamma323*gammado313 + 
     gamma312*gammado333 - 0.5*ddg1323[ijk] + 
     0.5*(ddg1233[ijk] - ddg2313[ijk] + ddg3312[ijk])) + 
  ginv23*(-(gamma123*gammado123) + gamma122*gammado133 - 
     gamma223*gammado223 + gamma222*gammado233 - gamma323*gammado323 + 
     gamma322*gammado333 - ddg2323[ijk] + 0.5*(ddg2233[ijk] + ddg3322[ijk]))
;

R23
=
(gamma133*gammado123 - gamma123*gammado133 + gamma233*gammado223 - 
     gamma223*gammado233 + gamma333*gammado323 - gamma323*gammado333)*ginv33 \
+ ginv11*(-(gamma123*gammado111) + gamma113*gammado112 - 
     gamma223*gammado211 + gamma213*gammado212 - gamma323*gammado311 + 
     gamma313*gammado312 + 0.5*(ddg1213[ijk] + ddg1312[ijk]) + 
     0.5*(-ddg1123[ijk] - ddg2311[ijk])) + 
  ginv12*(-(gamma123*gammado112) + gamma113*gammado122 - 
     gamma223*gammado212 + gamma213*gammado222 - gamma323*gammado312 + 
     gamma313*gammado322 + 0.5*(ddg1322[ijk] + ddg2213[ijk]) + 
     0.5*(-ddg1223[ijk] - ddg2312[ijk])) + 
  ginv13*(gamma133*gammado112 + gamma113*gammado123 + gamma233*gammado212 + 
     gamma213*gammado223 + gamma333*gammado312 - 
     2.*(gamma123*gammado113 + gamma223*gammado213 + gamma323*gammado313) + 
     gamma313*gammado323 + 0.5*
      (ddg1233[ijk] - ddg1323[ijk] - ddg2313[ijk] + ddg3312[ijk])) + 
  ginv23*(gamma133*gammado122 - gamma123*gammado123 + gamma233*gammado222 - 
     gamma223*gammado223 + gamma333*gammado322 - gamma323*gammado323 - 
     ddg2323[ijk] + 0.5*(ddg2233[ijk] + ddg3322[ijk]))
;

R33
=
(-(gamma133*gammado113) + gamma113*gammado133 - gamma233*gammado213 + 
     gamma213*gammado233 - gamma333*gammado313 + gamma313*gammado333)*ginv13 \
+ (-(gamma133*gammado123) + gamma123*gammado133 - gamma233*gammado223 + 
     gamma223*gammado233 - gamma333*gammado323 + gamma323*gammado333)*ginv23 \
+ ginv11*(-(gamma133*gammado111) + gamma113*gammado113 - 
     gamma233*gammado211 + gamma213*gammado213 - gamma333*gammado311 + 
     gamma313*gammado313 + ddg1313[ijk] + 0.5*(-ddg1133[ijk] - ddg3311[ijk])\
) + ginv12*(gamma123*gammado113 + gamma113*gammado123 + 
     gamma223*gammado213 + gamma213*gammado223 - 
     2.*(gamma133*gammado112 + gamma233*gammado212 + gamma333*gammado312) + 
     gamma323*gammado313 + gamma313*gammado323 - ddg1233[ijk] + 
     ddg1323[ijk] + ddg2313[ijk] - ddg3312[ijk]) + 
  ginv22*(-(gamma133*gammado122) + gamma123*gammado123 - 
     gamma233*gammado222 + gamma223*gammado223 - gamma333*gammado322 + 
     gamma323*gammado323 + ddg2323[ijk] + 0.5*(-ddg2233[ijk] - ddg3322[ijk]))
;

R
=
ginv11*R11 + ginv22*R22 + 2.*(ginv12*R12 + ginv13*R13 + ginv23*R23) + 
  ginv33*R33
;

KK11
=
2.*(ginv23*K12[ijk]*K13[ijk] + 
     K11[ijk]*(ginv12*K12[ijk] + ginv13*K13[ijk])) + ginv11*pow2(K11[ijk]) + 
  ginv22*pow2(K12[ijk]) + ginv33*pow2(K13[ijk])
;

KK12
=
K12[ijk]*(ginv11*K11[ijk] + ginv22*K22[ijk]) + ginv33*K13[ijk]*K23[ijk] + 
  ginv13*(K12[ijk]*K13[ijk] + K11[ijk]*K23[ijk]) + 
  ginv23*(K13[ijk]*K22[ijk] + K12[ijk]*K23[ijk]) + 
  ginv12*(K11[ijk]*K22[ijk] + pow2(K12[ijk]))
;

KK13
=
ginv22*K12[ijk]*K23[ijk] + ginv12*(K12[ijk]*K13[ijk] + K11[ijk]*K23[ijk]) + 
  K13[ijk]*(ginv11*K11[ijk] + ginv33*K33[ijk]) + 
  ginv23*(K13[ijk]*K23[ijk] + K12[ijk]*K33[ijk]) + 
  ginv13*(K11[ijk]*K33[ijk] + pow2(K13[ijk]))
;

KK22
=
2.*(ginv23*K22[ijk]*K23[ijk] + 
     K12[ijk]*(ginv12*K22[ijk] + ginv13*K23[ijk])) + ginv11*pow2(K12[ijk]) + 
  ginv22*pow2(K22[ijk]) + ginv33*pow2(K23[ijk])
;

KK23
=
ginv11*K12[ijk]*K13[ijk] + ginv12*(K13[ijk]*K22[ijk] + K12[ijk]*K23[ijk]) + 
  K23[ijk]*(ginv22*K22[ijk] + ginv33*K33[ijk]) + 
  ginv13*(K13[ijk]*K23[ijk] + K12[ijk]*K33[ijk]) + 
  ginv23*(K22[ijk]*K33[ijk] + pow2(K23[ijk]))
;

KK33
=
2.*(ginv23*K23[ijk]*K33[ijk] + 
     K13[ijk]*(ginv12*K23[ijk] + ginv13*K33[ijk])) + ginv11*pow2(K13[ijk]) + 
  ginv22*pow2(K23[ijk]) + ginv33*pow2(K33[ijk])
;

trK
=
ginv11*K11[ijk] + ginv22*K22[ijk] + 
  2.*(ginv12*K12[ijk] + ginv13*K13[ijk] + ginv23*K23[ijk]) + ginv33*K33[ijk]
;

trKK
=
ginv11*KK11 + ginv22*KK22 + 2.*(ginv12*KK12 + ginv13*KK13 + ginv23*KK23) + 
  ginv33*KK33
;

hamil
=
R - trKK + pow2(trK)
;

R
=
R - hamil*RtoRminusHfactor
;

dginv111
=
-2.*(ginv11*(ginv12*dg121[ijk] + ginv13*dg131[ijk]) + 
     ginv12*ginv13*dg231[ijk]) - dg111[ijk]*pow2(ginv11) - 
  dg221[ijk]*pow2(ginv12) - dg331[ijk]*pow2(ginv13)
;

dginv112
=
-2.*(ginv11*(ginv12*dg122[ijk] + ginv13*dg132[ijk]) + 
     ginv12*ginv13*dg232[ijk]) - dg112[ijk]*pow2(ginv11) - 
  dg222[ijk]*pow2(ginv12) - dg332[ijk]*pow2(ginv13)
;

dginv113
=
-2.*(ginv11*(ginv12*dg123[ijk] + ginv13*dg133[ijk]) + 
     ginv12*ginv13*dg233[ijk]) - dg113[ijk]*pow2(ginv11) - 
  dg223[ijk]*pow2(ginv12) - dg333[ijk]*pow2(ginv13)
;

dginv121
=
-(ginv11*(ginv12*dg111[ijk] + ginv22*dg121[ijk] + ginv23*dg131[ijk])) - 
  ginv12*(ginv13*dg131[ijk] + ginv22*dg221[ijk] + ginv23*dg231[ijk]) - 
  ginv13*(ginv22*dg231[ijk] + ginv23*dg331[ijk]) - dg121[ijk]*pow2(ginv12)
;

dginv122
=
-(ginv11*(ginv12*dg112[ijk] + ginv22*dg122[ijk] + ginv23*dg132[ijk])) - 
  ginv12*(ginv13*dg132[ijk] + ginv22*dg222[ijk] + ginv23*dg232[ijk]) - 
  ginv13*(ginv22*dg232[ijk] + ginv23*dg332[ijk]) - dg122[ijk]*pow2(ginv12)
;

dginv123
=
-(ginv11*(ginv12*dg113[ijk] + ginv22*dg123[ijk] + ginv23*dg133[ijk])) - 
  ginv12*(ginv13*dg133[ijk] + ginv22*dg223[ijk] + ginv23*dg233[ijk]) - 
  ginv13*(ginv22*dg233[ijk] + ginv23*dg333[ijk]) - dg123[ijk]*pow2(ginv12)
;

dginv131
=
-(ginv11*(ginv13*dg111[ijk] + ginv23*dg121[ijk] + ginv33*dg131[ijk])) - 
  ginv12*(ginv13*dg121[ijk] + ginv23*dg221[ijk] + ginv33*dg231[ijk]) - 
  ginv13*(ginv23*dg231[ijk] + ginv33*dg331[ijk]) - dg131[ijk]*pow2(ginv13)
;

dginv132
=
-(ginv11*(ginv13*dg112[ijk] + ginv23*dg122[ijk] + ginv33*dg132[ijk])) - 
  ginv12*(ginv13*dg122[ijk] + ginv23*dg222[ijk] + ginv33*dg232[ijk]) - 
  ginv13*(ginv23*dg232[ijk] + ginv33*dg332[ijk]) - dg132[ijk]*pow2(ginv13)
;

dginv133
=
-(ginv11*(ginv13*dg113[ijk] + ginv23*dg123[ijk] + ginv33*dg133[ijk])) - 
  ginv12*(ginv13*dg123[ijk] + ginv23*dg223[ijk] + ginv33*dg233[ijk]) - 
  ginv13*(ginv23*dg233[ijk] + ginv33*dg333[ijk]) - dg133[ijk]*pow2(ginv13)
;

dginv221
=
-2.*(ginv12*(ginv22*dg121[ijk] + ginv23*dg131[ijk]) + 
     ginv22*ginv23*dg231[ijk]) - dg111[ijk]*pow2(ginv12) - 
  dg221[ijk]*pow2(ginv22) - dg331[ijk]*pow2(ginv23)
;

dginv222
=
-2.*(ginv12*(ginv22*dg122[ijk] + ginv23*dg132[ijk]) + 
     ginv22*ginv23*dg232[ijk]) - dg112[ijk]*pow2(ginv12) - 
  dg222[ijk]*pow2(ginv22) - dg332[ijk]*pow2(ginv23)
;

dginv223
=
-2.*(ginv12*(ginv22*dg123[ijk] + ginv23*dg133[ijk]) + 
     ginv22*ginv23*dg233[ijk]) - dg113[ijk]*pow2(ginv12) - 
  dg223[ijk]*pow2(ginv22) - dg333[ijk]*pow2(ginv23)
;

dginv231
=
-(ginv13*(ginv22*dg121[ijk] + ginv23*dg131[ijk])) - 
  ginv12*(ginv13*dg111[ijk] + ginv23*dg121[ijk] + ginv33*dg131[ijk]) - 
  ginv22*(ginv23*dg221[ijk] + ginv33*dg231[ijk]) - 
  ginv23*ginv33*dg331[ijk] - dg231[ijk]*pow2(ginv23)
;

dginv232
=
-(ginv13*(ginv22*dg122[ijk] + ginv23*dg132[ijk])) - 
  ginv12*(ginv13*dg112[ijk] + ginv23*dg122[ijk] + ginv33*dg132[ijk]) - 
  ginv22*(ginv23*dg222[ijk] + ginv33*dg232[ijk]) - 
  ginv23*ginv33*dg332[ijk] - dg232[ijk]*pow2(ginv23)
;

dginv233
=
-(ginv13*(ginv22*dg123[ijk] + ginv23*dg133[ijk])) - 
  ginv12*(ginv13*dg113[ijk] + ginv23*dg123[ijk] + ginv33*dg133[ijk]) - 
  ginv22*(ginv23*dg223[ijk] + ginv33*dg233[ijk]) - 
  ginv23*ginv33*dg333[ijk] - dg233[ijk]*pow2(ginv23)
;

dginv331
=
-2.*(ginv13*(ginv23*dg121[ijk] + ginv33*dg131[ijk]) + 
     ginv23*ginv33*dg231[ijk]) - dg111[ijk]*pow2(ginv13) - 
  dg221[ijk]*pow2(ginv23) - dg331[ijk]*pow2(ginv33)
;

dginv332
=
-2.*(ginv13*(ginv23*dg122[ijk] + ginv33*dg132[ijk]) + 
     ginv23*ginv33*dg232[ijk]) - dg112[ijk]*pow2(ginv13) - 
  dg222[ijk]*pow2(ginv23) - dg332[ijk]*pow2(ginv33)
;

dginv333
=
-2.*(ginv13*(ginv23*dg123[ijk] + ginv33*dg133[ijk]) + 
     ginv23*ginv33*dg233[ijk]) - dg113[ijk]*pow2(ginv13) - 
  dg223[ijk]*pow2(ginv23) - dg333[ijk]*pow2(ginv33)
;

dK1
=
ginv11*dK111[ijk] + ginv22*dK221[ijk] + ginv33*dK331[ijk] + 
  dginv111*K11[ijk] + dginv221*K22[ijk] + 
  2.*(ginv12*dK121[ijk] + ginv13*dK131[ijk] + ginv23*dK231[ijk] + 
     dginv121*K12[ijk] + dginv131*K13[ijk] + dginv231*K23[ijk]) + 
  dginv331*K33[ijk]
;

dK2
=
ginv11*dK112[ijk] + ginv22*dK222[ijk] + ginv33*dK332[ijk] + 
  dginv112*K11[ijk] + dginv222*K22[ijk] + 
  2.*(ginv12*dK122[ijk] + ginv13*dK132[ijk] + ginv23*dK232[ijk] + 
     dginv122*K12[ijk] + dginv132*K13[ijk] + dginv232*K23[ijk]) + 
  dginv332*K33[ijk]
;

dK3
=
ginv11*dK113[ijk] + ginv22*dK223[ijk] + ginv33*dK333[ijk] + 
  dginv113*K11[ijk] + dginv223*K22[ijk] + 
  2.*(ginv12*dK123[ijk] + ginv13*dK133[ijk] + ginv23*dK233[ijk] + 
     dginv123*K12[ijk] + dginv133*K13[ijk] + dginv233*K23[ijk]) + 
  dginv333*K33[ijk]
;

cddalp11
=
-(gamma111*dalp1[ijk]) - gamma211*dalp2[ijk] - gamma311*dalp3[ijk] + 
  ddalp11[ijk]
;

cddalp12
=
-(gamma112*dalp1[ijk]) - gamma212*dalp2[ijk] - gamma312*dalp3[ijk] + 
  ddalp12[ijk]
;

cddalp13
=
-(gamma113*dalp1[ijk]) - gamma213*dalp2[ijk] - gamma313*dalp3[ijk] + 
  ddalp13[ijk]
;

cddalp22
=
-(gamma122*dalp1[ijk]) - gamma222*dalp2[ijk] - gamma322*dalp3[ijk] + 
  ddalp22[ijk]
;

cddalp23
=
-(gamma123*dalp1[ijk]) - gamma223*dalp2[ijk] - gamma323*dalp3[ijk] + 
  ddalp23[ijk]
;

cddalp33
=
-(gamma133*dalp1[ijk]) - gamma233*dalp2[ijk] - gamma333*dalp3[ijk] + 
  ddalp33[ijk]
;

cdZ11
=
dZ11[ijk] - gamma111*Z1[ijk] - gamma211*Z2[ijk] - gamma311*Z3[ijk]
;

cdZ12
=
dZ12[ijk] - gamma112*Z1[ijk] - gamma212*Z2[ijk] - gamma312*Z3[ijk]
;

cdZ13
=
dZ13[ijk] - gamma113*Z1[ijk] - gamma213*Z2[ijk] - gamma313*Z3[ijk]
;

cdZ21
=
dZ21[ijk] - gamma112*Z1[ijk] - gamma212*Z2[ijk] - gamma312*Z3[ijk]
;

cdZ22
=
dZ22[ijk] - gamma122*Z1[ijk] - gamma222*Z2[ijk] - gamma322*Z3[ijk]
;

cdZ23
=
dZ23[ijk] - gamma123*Z1[ijk] - gamma223*Z2[ijk] - gamma323*Z3[ijk]
;

cdZ31
=
dZ31[ijk] - gamma113*Z1[ijk] - gamma213*Z2[ijk] - gamma313*Z3[ijk]
;

cdZ32
=
dZ32[ijk] - gamma123*Z1[ijk] - gamma223*Z2[ijk] - gamma323*Z3[ijk]
;

cdZ33
=
dZ33[ijk] - gamma133*Z1[ijk] - gamma233*Z2[ijk] - gamma333*Z3[ijk]
;

cdK111
=
dK111[ijk] - 2.*(gamma111*K11[ijk] + gamma211*K12[ijk] + gamma311*K13[ijk])
;

cdK112
=
dK112[ijk] - gamma112*K11[ijk] - (gamma111 + gamma212)*K12[ijk] - 
  gamma312*K13[ijk] - gamma211*K22[ijk] - gamma311*K23[ijk]
;

cdK113
=
dK113[ijk] - gamma113*K11[ijk] - gamma213*K12[ijk] - 
  (gamma111 + gamma313)*K13[ijk] - gamma211*K23[ijk] - gamma311*K33[ijk]
;

cdK121
=
dK121[ijk] - gamma112*K11[ijk] - (gamma111 + gamma212)*K12[ijk] - 
  gamma312*K13[ijk] - gamma211*K22[ijk] - gamma311*K23[ijk]
;

cdK122
=
dK122[ijk] - 2.*(gamma112*K12[ijk] + gamma212*K22[ijk] + gamma312*K23[ijk])
;

cdK123
=
dK123[ijk] - gamma113*K12[ijk] - gamma112*K13[ijk] - gamma213*K22[ijk] - 
  (gamma212 + gamma313)*K23[ijk] - gamma312*K33[ijk]
;

cdK131
=
dK131[ijk] - gamma113*K11[ijk] - gamma213*K12[ijk] - 
  (gamma111 + gamma313)*K13[ijk] - gamma211*K23[ijk] - gamma311*K33[ijk]
;

cdK132
=
dK132[ijk] - gamma113*K12[ijk] - gamma112*K13[ijk] - gamma213*K22[ijk] - 
  (gamma212 + gamma313)*K23[ijk] - gamma312*K33[ijk]
;

cdK133
=
dK133[ijk] - 2.*(gamma113*K13[ijk] + gamma213*K23[ijk] + gamma313*K33[ijk])
;

cdK211
=
dK121[ijk] - 2.*(gamma112*K11[ijk] + gamma212*K12[ijk] + gamma312*K13[ijk])
;

cdK212
=
dK122[ijk] - gamma122*K11[ijk] - (gamma112 + gamma222)*K12[ijk] - 
  gamma322*K13[ijk] - gamma212*K22[ijk] - gamma312*K23[ijk]
;

cdK213
=
dK123[ijk] - gamma123*K11[ijk] - gamma223*K12[ijk] - 
  (gamma112 + gamma323)*K13[ijk] - gamma212*K23[ijk] - gamma312*K33[ijk]
;

cdK221
=
dK221[ijk] - gamma122*K11[ijk] - (gamma112 + gamma222)*K12[ijk] - 
  gamma322*K13[ijk] - gamma212*K22[ijk] - gamma312*K23[ijk]
;

cdK222
=
dK222[ijk] - 2.*(gamma122*K12[ijk] + gamma222*K22[ijk] + gamma322*K23[ijk])
;

cdK223
=
dK223[ijk] - gamma123*K12[ijk] - gamma122*K13[ijk] - gamma223*K22[ijk] - 
  (gamma222 + gamma323)*K23[ijk] - gamma322*K33[ijk]
;

cdK231
=
dK231[ijk] - gamma123*K11[ijk] - gamma223*K12[ijk] - 
  (gamma112 + gamma323)*K13[ijk] - gamma212*K23[ijk] - gamma312*K33[ijk]
;

cdK232
=
dK232[ijk] - gamma123*K12[ijk] - gamma122*K13[ijk] - gamma223*K22[ijk] - 
  (gamma222 + gamma323)*K23[ijk] - gamma322*K33[ijk]
;

cdK233
=
dK233[ijk] - 2.*(gamma123*K13[ijk] + gamma223*K23[ijk] + gamma323*K33[ijk])
;

cdK311
=
dK131[ijk] - 2.*(gamma113*K11[ijk] + gamma213*K12[ijk] + gamma313*K13[ijk])
;

cdK312
=
dK132[ijk] - gamma123*K11[ijk] - (gamma113 + gamma223)*K12[ijk] - 
  gamma323*K13[ijk] - gamma213*K22[ijk] - gamma313*K23[ijk]
;

cdK313
=
dK133[ijk] - gamma133*K11[ijk] - gamma233*K12[ijk] - 
  (gamma113 + gamma333)*K13[ijk] - gamma213*K23[ijk] - gamma313*K33[ijk]
;

cdK321
=
dK231[ijk] - gamma123*K11[ijk] - (gamma113 + gamma223)*K12[ijk] - 
  gamma323*K13[ijk] - gamma213*K22[ijk] - gamma313*K23[ijk]
;

cdK322
=
dK232[ijk] - 2.*(gamma123*K12[ijk] + gamma223*K22[ijk] + gamma323*K23[ijk])
;

cdK323
=
dK233[ijk] - gamma133*K12[ijk] - gamma123*K13[ijk] - gamma233*K22[ijk] - 
  (gamma223 + gamma333)*K23[ijk] - gamma323*K33[ijk]
;

cdK331
=
dK331[ijk] - gamma133*K11[ijk] - gamma233*K12[ijk] - 
  (gamma113 + gamma333)*K13[ijk] - gamma213*K23[ijk] - gamma313*K33[ijk]
;

cdK332
=
dK332[ijk] - gamma133*K12[ijk] - gamma123*K13[ijk] - gamma233*K22[ijk] - 
  (gamma223 + gamma333)*K23[ijk] - gamma323*K33[ijk]
;

cdK333
=
dK333[ijk] - 2.*(gamma133*K13[ijk] + gamma233*K23[ijk] + gamma333*K33[ijk])
;

lieg11
=
betadg11 + 2.*(dbeta11[ijk]*g11[ijk] + dbeta21[ijk]*g12[ijk] + 
     dbeta31[ijk]*g13[ijk])
;

lieg12
=
betadg12 + dbeta12[ijk]*g11[ijk] + (dbeta11[ijk] + dbeta22[ijk])*g12[ijk] + 
  dbeta32[ijk]*g13[ijk] + dbeta21[ijk]*g22[ijk] + dbeta31[ijk]*g23[ijk]
;

lieg12
=
betadg21 + dbeta12[ijk]*g11[ijk] + (dbeta11[ijk] + dbeta22[ijk])*g12[ijk] + 
  dbeta32[ijk]*g13[ijk] + dbeta21[ijk]*g22[ijk] + dbeta31[ijk]*g23[ijk]
;

lieg13
=
betadg13 + dbeta13[ijk]*g11[ijk] + dbeta23[ijk]*g12[ijk] + 
  (dbeta11[ijk] + dbeta33[ijk])*g13[ijk] + dbeta21[ijk]*g23[ijk] + 
  dbeta31[ijk]*g33[ijk]
;

lieg13
=
betadg31 + dbeta13[ijk]*g11[ijk] + dbeta23[ijk]*g12[ijk] + 
  (dbeta11[ijk] + dbeta33[ijk])*g13[ijk] + dbeta21[ijk]*g23[ijk] + 
  dbeta31[ijk]*g33[ijk]
;

lieg22
=
betadg22 + 2.*(dbeta12[ijk]*g12[ijk] + dbeta22[ijk]*g22[ijk] + 
     dbeta32[ijk]*g23[ijk])
;

lieg23
=
betadg23 + dbeta13[ijk]*g12[ijk] + dbeta12[ijk]*g13[ijk] + 
  dbeta23[ijk]*g22[ijk] + (dbeta22[ijk] + dbeta33[ijk])*g23[ijk] + 
  dbeta32[ijk]*g33[ijk]
;

lieg23
=
betadg32 + dbeta13[ijk]*g12[ijk] + dbeta12[ijk]*g13[ijk] + 
  dbeta23[ijk]*g22[ijk] + (dbeta22[ijk] + dbeta33[ijk])*g23[ijk] + 
  dbeta32[ijk]*g33[ijk]
;

lieg33
=
betadg33 + 2.*(dbeta13[ijk]*g13[ijk] + dbeta23[ijk]*g23[ijk] + 
     dbeta33[ijk]*g33[ijk])
;

lieK11
=
betadK11 + 2.*(dbeta11[ijk]*K11[ijk] + dbeta21[ijk]*K12[ijk] + 
     dbeta31[ijk]*K13[ijk])
;

lieK12
=
betadK12 + dbeta12[ijk]*K11[ijk] + (dbeta11[ijk] + dbeta22[ijk])*K12[ijk] + 
  dbeta32[ijk]*K13[ijk] + dbeta21[ijk]*K22[ijk] + dbeta31[ijk]*K23[ijk]
;

lieK12
=
betadK21 + dbeta12[ijk]*K11[ijk] + (dbeta11[ijk] + dbeta22[ijk])*K12[ijk] + 
  dbeta32[ijk]*K13[ijk] + dbeta21[ijk]*K22[ijk] + dbeta31[ijk]*K23[ijk]
;

lieK13
=
betadK13 + dbeta13[ijk]*K11[ijk] + dbeta23[ijk]*K12[ijk] + 
  (dbeta11[ijk] + dbeta33[ijk])*K13[ijk] + dbeta21[ijk]*K23[ijk] + 
  dbeta31[ijk]*K33[ijk]
;

lieK13
=
betadK31 + dbeta13[ijk]*K11[ijk] + dbeta23[ijk]*K12[ijk] + 
  (dbeta11[ijk] + dbeta33[ijk])*K13[ijk] + dbeta21[ijk]*K23[ijk] + 
  dbeta31[ijk]*K33[ijk]
;

lieK22
=
betadK22 + 2.*(dbeta12[ijk]*K12[ijk] + dbeta22[ijk]*K22[ijk] + 
     dbeta32[ijk]*K23[ijk])
;

lieK23
=
betadK23 + dbeta13[ijk]*K12[ijk] + dbeta12[ijk]*K13[ijk] + 
  dbeta23[ijk]*K22[ijk] + (dbeta22[ijk] + dbeta33[ijk])*K23[ijk] + 
  dbeta32[ijk]*K33[ijk]
;

lieK23
=
betadK32 + dbeta13[ijk]*K12[ijk] + dbeta12[ijk]*K13[ijk] + 
  dbeta23[ijk]*K22[ijk] + (dbeta22[ijk] + dbeta33[ijk])*K23[ijk] + 
  dbeta32[ijk]*K33[ijk]
;

lieK33
=
betadK33 + 2.*(dbeta13[ijk]*K13[ijk] + dbeta23[ijk]*K23[ijk] + 
     dbeta33[ijk]*K33[ijk])
;

lieTheta
=
betadTheta
;

lieZ1
=
betadZ1 + dbeta11[ijk]*Z1[ijk] + dbeta21[ijk]*Z2[ijk] + dbeta31[ijk]*Z3[ijk]
;

lieZ2
=
betadZ2 + dbeta12[ijk]*Z1[ijk] + dbeta22[ijk]*Z2[ijk] + dbeta32[ijk]*Z3[ijk]
;

lieZ3
=
betadZ3 + dbeta13[ijk]*Z1[ijk] + dbeta23[ijk]*Z2[ijk] + dbeta33[ijk]*Z3[ijk]
;

rg11
=
lieg11 - 2.*alpha[ijk]*K11[ijk]
;

rg12
=
lieg12 - 2.*alpha[ijk]*K12[ijk]
;

rg13
=
lieg13 - 2.*alpha[ijk]*K13[ijk]
;

rg22
=
lieg22 - 2.*alpha[ijk]*K22[ijk]
;

rg23
=
lieg23 - 2.*alpha[ijk]*K23[ijk]
;

rg33
=
lieg33 - 2.*alpha[ijk]*K33[ijk]
;

rK11
=
-cddalp11 + lieK11 + alpha[ijk]*
   (R11 + trK*K11[ijk] + 2.*(cdZ11 - KK11 - K11[ijk]*Theta[ijk]))
;

rK12
=
-cddalp12 + lieK12 + alpha[ijk]*
   (cdZ12 + cdZ21 + R12 + trK*K12[ijk] - 2.*(KK12 + K12[ijk]*Theta[ijk]))
;

rK13
=
-cddalp13 + lieK13 + alpha[ijk]*
   (cdZ13 + cdZ31 + R13 + trK*K13[ijk] - 2.*(KK13 + K13[ijk]*Theta[ijk]))
;

rK22
=
-cddalp22 + lieK22 + alpha[ijk]*
   (R22 + trK*K22[ijk] + 2.*(cdZ22 - KK22 - K22[ijk]*Theta[ijk]))
;

rK23
=
-cddalp23 + lieK23 + alpha[ijk]*
   (cdZ23 + cdZ32 + R23 + trK*K23[ijk] - 2.*(KK23 + K23[ijk]*Theta[ijk]))
;

rK33
=
-cddalp33 + lieK33 + alpha[ijk]*
   (R33 + trK*K33[ijk] + 2.*(cdZ33 - KK33 - K33[ijk]*Theta[ijk]))
;

rTheta
=
lieTheta + ginv11*(cdZ11*alpha[ijk] - dalp1[ijk]*Z1[ijk]) + 
  ginv12*((cdZ12 + cdZ21)*alpha[ijk] - dalp2[ijk]*Z1[ijk] - 
     dalp1[ijk]*Z2[ijk]) + ginv22*(cdZ22*alpha[ijk] - dalp2[ijk]*Z2[ijk]) + 
  ginv13*((cdZ13 + cdZ31)*alpha[ijk] - dalp3[ijk]*Z1[ijk] - 
     dalp1[ijk]*Z3[ijk]) + ginv23*
   ((cdZ23 + cdZ32)*alpha[ijk] - dalp3[ijk]*Z2[ijk] - dalp2[ijk]*Z3[ijk]) + 
  ginv33*(cdZ33*alpha[ijk] - dalp3[ijk]*Z3[ijk]) + 
  alpha[ijk]*(-(trK*Theta[ijk]) + 0.5*(R - trKK + pow2(trK)))
;

rZ1
=
lieZ1 - dalp1[ijk]*Theta[ijk] + 
  alpha[ijk]*(-dK1 + dTheta1[ijk] + ginv11*(cdK111 - 2.*K11[ijk]*Z1[ijk]) + 
     ginv22*(cdK212 - 2.*K12[ijk]*Z2[ijk]) + 
     ginv12*(cdK112 + cdK211 - 2.*(K12[ijk]*Z1[ijk] + K11[ijk]*Z2[ijk])) + 
     ginv33*(cdK313 - 2.*K13[ijk]*Z3[ijk]) + 
     ginv13*(cdK113 + cdK311 - 2.*(K13[ijk]*Z1[ijk] + K11[ijk]*Z3[ijk])) + 
     ginv23*(cdK213 + cdK312 - 2.*(K13[ijk]*Z2[ijk] + K12[ijk]*Z3[ijk])))
;

rZ2
=
lieZ2 - dalp2[ijk]*Theta[ijk] + 
  alpha[ijk]*(-dK2 + dTheta2[ijk] + ginv11*(cdK121 - 2.*K12[ijk]*Z1[ijk]) + 
     ginv22*(cdK222 - 2.*K22[ijk]*Z2[ijk]) + 
     ginv12*(cdK122 + cdK221 - 2.*(K22[ijk]*Z1[ijk] + K12[ijk]*Z2[ijk])) + 
     ginv33*(cdK323 - 2.*K23[ijk]*Z3[ijk]) + 
     ginv13*(cdK123 + cdK321 - 2.*(K23[ijk]*Z1[ijk] + K12[ijk]*Z3[ijk])) + 
     ginv23*(cdK223 + cdK322 - 2.*(K23[ijk]*Z2[ijk] + K22[ijk]*Z3[ijk])))
;

rZ3
=
lieZ3 - dalp3[ijk]*Theta[ijk] + 
  alpha[ijk]*(-dK3 + dTheta3[ijk] + ginv11*(cdK131 - 2.*K13[ijk]*Z1[ijk]) + 
     ginv22*(cdK232 - 2.*K23[ijk]*Z2[ijk]) + 
     ginv12*(cdK132 + cdK231 - 2.*(K23[ijk]*Z1[ijk] + K13[ijk]*Z2[ijk])) + 
     ginv33*(cdK333 - 2.*K33[ijk]*Z3[ijk]) + 
     ginv13*(cdK133 + cdK331 - 2.*(K33[ijk]*Z1[ijk] + K13[ijk]*Z3[ijk])) + 
     ginv23*(cdK233 + cdK332 - 2.*(K33[ijk]*Z2[ijk] + K23[ijk]*Z3[ijk])))
;

liealpha
=
betadalp
;

ralpha0
=
alpha[ijk]*Power(psi[ijk],lapsepsipower)*
  (-trK + subtractK0*K0[ijk] + lapseharmonicm*Theta[ijk])
;

ralpha
=
nonconstantlapse*(liealpha*withshift + 
    ralpha0*(lapseharmonicf*oploglapse + harmoniclapse*alpha[ijk]))
;

rbeta1
=
0
;

rbeta2
=
0
;

rbeta3
=
0
;

rB1
=
0
;

rB2
=
0
;

rB3
=
0
;



/* conditional */
if (addlinear) {

ng11[ijk]
=
dt*rg11 + pg11[ijk]
;

ng12[ijk]
=
dt*rg12 + pg12[ijk]
;

ng13[ijk]
=
dt*rg13 + pg13[ijk]
;

ng22[ijk]
=
dt*rg22 + pg22[ijk]
;

ng23[ijk]
=
dt*rg23 + pg23[ijk]
;

ng33[ijk]
=
dt*rg33 + pg33[ijk]
;

nK11[ijk]
=
dt*rK11 + pK11[ijk]
;

nK12[ijk]
=
dt*rK12 + pK12[ijk]
;

nK13[ijk]
=
dt*rK13 + pK13[ijk]
;

nK22[ijk]
=
dt*rK22 + pK22[ijk]
;

nK23[ijk]
=
dt*rK23 + pK23[ijk]
;

nK33[ijk]
=
dt*rK33 + pK33[ijk]
;

nZ1[ijk]
=
dt*rZ1 + pZ1[ijk]
;

nZ2[ijk]
=
dt*rZ2 + pZ2[ijk]
;

nZ3[ijk]
=
dt*rZ3 + pZ3[ijk]
;

nTheta[ijk]
=
dt*rTheta + pTheta[ijk]
;

nalpha[ijk]
=
dt*ralpha + palpha[ijk]
;

nbeta1[ijk]
=
dt*rbeta1 + pbeta1[ijk]
;

nbeta2[ijk]
=
dt*rbeta2 + pbeta2[ijk]
;

nbeta3[ijk]
=
dt*rbeta3 + pbeta3[ijk]
;

nB1[ijk]
=
dt*rB1 + pB1[ijk]
;

nB2[ijk]
=
dt*rB2 + pB2[ijk]
;

nB3[ijk]
=
dt*rB3 + pB3[ijk]
;


} else { /* if (!addlinear) */

ng11[ijk]
=
rg11
;

ng12[ijk]
=
rg12
;

ng13[ijk]
=
rg13
;

ng22[ijk]
=
rg22
;

ng23[ijk]
=
rg23
;

ng33[ijk]
=
rg33
;

nK11[ijk]
=
rK11
;

nK12[ijk]
=
rK12
;

nK13[ijk]
=
rK13
;

nK22[ijk]
=
rK22
;

nK23[ijk]
=
rK23
;

nK33[ijk]
=
rK33
;

nZ1[ijk]
=
rZ1
;

nZ2[ijk]
=
rZ2
;

nZ3[ijk]
=
rZ3
;

nTheta[ijk]
=
rTheta
;

nalpha[ijk]
=
ralpha
;

nbeta1[ijk]
=
rbeta1
;

nbeta2[ijk]
=
rbeta2
;

nbeta3[ijk]
=
rbeta3
;

nB1[ijk]
=
rB1
;

nB2[ijk]
=
rB2
;

nB3[ijk]
=
rB3
;

}
/* if (addlinear) */


} /* end of points */
} /* end of boxes */


}  /* end of function */

/* Z4secondO_rhs.c */
/* nvars = 200, n* = 1537,  n/ = 30,  n+ = 1497, n = 3064, O = 1 */
