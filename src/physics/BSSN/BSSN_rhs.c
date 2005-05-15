/* BSSN_rhs.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 14.5.2005 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BSSN.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BSSN_rhs(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)
{
tGrid *grid = ucur->grid;
int bi;

int addlinear = (dt != 0.0l);
int usepsi = 1;
double forceKzerofactor = Getv("BSSN_forceKzero", "no");
int subtractA           = Getv("BSSN_subtractA", "yes");
int normalizedetg       = Getv("BSSN_normalizedetg", "yes");
int nonconstantlapse    =!Getv("BSSN_lapse", "constant");
int oploglapse          = Getv("BSSN_lapse", "1+log");
int oploglapse2         = Getv("BSSN_lapse", "1+log2");
int oplogwithshift      = Getv("BSSN_lapse", "withshift");
int harmoniclapse       = Getv("BSSN_lapse", "harmonic");
int subtractK0          = Getv("BSSN_subtractK0", "yes");
int densitizedLapse = !Getv("BSSN_densitizedLapse", "no");
int densitizedoplogWithoutShift = Getv("BSSN_densitizedLapse", "1+log_withoutShift");
double alphaDensityWeight = Getd("BSSN_alphaDensityWeight");
double gamma0factor    = Getv("BSSN_shift", "gamma0");
double gamma2factor    = Getv("BSSN_shift", "gamma2");
double shift_st        = Getd("BSSN_shift_stop_time");
double evolveshift     = Cal( (grid->time>=shift_st)*(shift_st>=0), 0.0, 1.0);
double lapsepsipower   = Getd("BSSN_lapsepsipower");
double lapseharmonicf  = Getd("BSSN_lapseharmonicf");
double shiftpsipower   = Getd("BSSN_shiftpsipower");
double shiftalphapower = Getd("BSSN_shiftalphapower");
double shiftgammacoeff = Getd("BSSN_shiftgammacoeff");
double shiftdriver     = Getd("BSSN_shiftdriver");

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
double *phi = vlldataptr(ucur, box, 6);
double *A11 = vlldataptr(ucur, box, 7);
double *A12 = vlldataptr(ucur, box, 8);
double *A13 = vlldataptr(ucur, box, 9);
double *A22 = vlldataptr(ucur, box, 10);
double *A23 = vlldataptr(ucur, box, 11);
double *A33 = vlldataptr(ucur, box, 12);
double *K = vlldataptr(ucur, box, 13);
double *G1 = vlldataptr(ucur, box, 14);
double *G2 = vlldataptr(ucur, box, 15);
double *G3 = vlldataptr(ucur, box, 16);
double *alpha = vlldataptr(ucur, box, 17);
double *beta1 = vlldataptr(ucur, box, 18);
double *beta2 = vlldataptr(ucur, box, 19);
double *beta3 = vlldataptr(ucur, box, 20);
double *B1 = vlldataptr(ucur, box, 21);
double *B2 = vlldataptr(ucur, box, 22);
double *B3 = vlldataptr(ucur, box, 23);
double *alphaDensity = vlldataptr(ucur, box, 24);
double *ng11 = vlldataptr(unew, box, 0);
double *ng12 = vlldataptr(unew, box, 1);
double *ng13 = vlldataptr(unew, box, 2);
double *ng22 = vlldataptr(unew, box, 3);
double *ng23 = vlldataptr(unew, box, 4);
double *ng33 = vlldataptr(unew, box, 5);
double *nphi = vlldataptr(unew, box, 6);
double *nA11 = vlldataptr(unew, box, 7);
double *nA12 = vlldataptr(unew, box, 8);
double *nA13 = vlldataptr(unew, box, 9);
double *nA22 = vlldataptr(unew, box, 10);
double *nA23 = vlldataptr(unew, box, 11);
double *nA33 = vlldataptr(unew, box, 12);
double *nK = vlldataptr(unew, box, 13);
double *nG1 = vlldataptr(unew, box, 14);
double *nG2 = vlldataptr(unew, box, 15);
double *nG3 = vlldataptr(unew, box, 16);
double *nalpha = vlldataptr(unew, box, 17);
double *nbeta1 = vlldataptr(unew, box, 18);
double *nbeta2 = vlldataptr(unew, box, 19);
double *nbeta3 = vlldataptr(unew, box, 20);
double *nB1 = vlldataptr(unew, box, 21);
double *nB2 = vlldataptr(unew, box, 22);
double *nB3 = vlldataptr(unew, box, 23);
double *nalphaDensity = vlldataptr(unew, box, 24);
double *pg11 = vlldataptr(upre, box, 0);
double *pg12 = vlldataptr(upre, box, 1);
double *pg13 = vlldataptr(upre, box, 2);
double *pg22 = vlldataptr(upre, box, 3);
double *pg23 = vlldataptr(upre, box, 4);
double *pg33 = vlldataptr(upre, box, 5);
double *pphi = vlldataptr(upre, box, 6);
double *pA11 = vlldataptr(upre, box, 7);
double *pA12 = vlldataptr(upre, box, 8);
double *pA13 = vlldataptr(upre, box, 9);
double *pA22 = vlldataptr(upre, box, 10);
double *pA23 = vlldataptr(upre, box, 11);
double *pA33 = vlldataptr(upre, box, 12);
double *pK = vlldataptr(upre, box, 13);
double *pG1 = vlldataptr(upre, box, 14);
double *pG2 = vlldataptr(upre, box, 15);
double *pG3 = vlldataptr(upre, box, 16);
double *palpha = vlldataptr(upre, box, 17);
double *pbeta1 = vlldataptr(upre, box, 18);
double *pbeta2 = vlldataptr(upre, box, 19);
double *pbeta3 = vlldataptr(upre, box, 20);
double *pB1 = vlldataptr(upre, box, 21);
double *pB2 = vlldataptr(upre, box, 22);
double *pB3 = vlldataptr(upre, box, 23);
double *palphaDensity = vlldataptr(upre, box, 24);
int index_g11 = (ucur)->index[0];
int index_g12 = (ucur)->index[1];
int index_g13 = (ucur)->index[2];
int index_g22 = (ucur)->index[3];
int index_g23 = (ucur)->index[4];
int index_g33 = (ucur)->index[5];
int index_phi = (ucur)->index[6];
int index_A11 = (ucur)->index[7];
int index_A12 = (ucur)->index[8];
int index_A13 = (ucur)->index[9];
int index_A22 = (ucur)->index[10];
int index_A23 = (ucur)->index[11];
int index_A33 = (ucur)->index[12];
int index_K = (ucur)->index[13];
int index_G1 = (ucur)->index[14];
int index_G2 = (ucur)->index[15];
int index_G3 = (ucur)->index[16];
int index_alpha = (ucur)->index[17];
int index_beta1 = (ucur)->index[18];
int index_beta2 = (ucur)->index[19];
int index_beta3 = (ucur)->index[20];
int index_B1 = (ucur)->index[21];
int index_B2 = (ucur)->index[22];
int index_B3 = (ucur)->index[23];
int index_alphaDensity = (ucur)->index[24];
int index_ADMvars_dgxxx = Ind("ADMvars_dgxxx");
double *dgt111 = box->v[index_ADMvars_dgxxx + 0];
double *dgt112 = box->v[index_ADMvars_dgxxx + 1];
double *dgt113 = box->v[index_ADMvars_dgxxx + 2];
double *dgt121 = box->v[index_ADMvars_dgxxx + 3];
double *dgt122 = box->v[index_ADMvars_dgxxx + 4];
double *dgt123 = box->v[index_ADMvars_dgxxx + 5];
double *dgt131 = box->v[index_ADMvars_dgxxx + 6];
double *dgt132 = box->v[index_ADMvars_dgxxx + 7];
double *dgt133 = box->v[index_ADMvars_dgxxx + 8];
double *dgt221 = box->v[index_ADMvars_dgxxx + 9];
double *dgt222 = box->v[index_ADMvars_dgxxx + 10];
double *dgt223 = box->v[index_ADMvars_dgxxx + 11];
double *dgt231 = box->v[index_ADMvars_dgxxx + 12];
double *dgt232 = box->v[index_ADMvars_dgxxx + 13];
double *dgt233 = box->v[index_ADMvars_dgxxx + 14];
double *dgt331 = box->v[index_ADMvars_dgxxx + 15];
double *dgt332 = box->v[index_ADMvars_dgxxx + 16];
double *dgt333 = box->v[index_ADMvars_dgxxx + 17];
int index_ADMvars_ddgxxxx = Ind("ADMvars_ddgxxxx");
double *ddgt1111 = box->v[index_ADMvars_ddgxxxx + 0];
double *ddgt1112 = box->v[index_ADMvars_ddgxxxx + 1];
double *ddgt1113 = box->v[index_ADMvars_ddgxxxx + 2];
double *ddgt1122 = box->v[index_ADMvars_ddgxxxx + 3];
double *ddgt1123 = box->v[index_ADMvars_ddgxxxx + 4];
double *ddgt1133 = box->v[index_ADMvars_ddgxxxx + 5];
double *ddgt1211 = box->v[index_ADMvars_ddgxxxx + 6];
double *ddgt1212 = box->v[index_ADMvars_ddgxxxx + 7];
double *ddgt1213 = box->v[index_ADMvars_ddgxxxx + 8];
double *ddgt1222 = box->v[index_ADMvars_ddgxxxx + 9];
double *ddgt1223 = box->v[index_ADMvars_ddgxxxx + 10];
double *ddgt1233 = box->v[index_ADMvars_ddgxxxx + 11];
double *ddgt1311 = box->v[index_ADMvars_ddgxxxx + 12];
double *ddgt1312 = box->v[index_ADMvars_ddgxxxx + 13];
double *ddgt1313 = box->v[index_ADMvars_ddgxxxx + 14];
double *ddgt1322 = box->v[index_ADMvars_ddgxxxx + 15];
double *ddgt1323 = box->v[index_ADMvars_ddgxxxx + 16];
double *ddgt1333 = box->v[index_ADMvars_ddgxxxx + 17];
double *ddgt2211 = box->v[index_ADMvars_ddgxxxx + 18];
double *ddgt2212 = box->v[index_ADMvars_ddgxxxx + 19];
double *ddgt2213 = box->v[index_ADMvars_ddgxxxx + 20];
double *ddgt2222 = box->v[index_ADMvars_ddgxxxx + 21];
double *ddgt2223 = box->v[index_ADMvars_ddgxxxx + 22];
double *ddgt2233 = box->v[index_ADMvars_ddgxxxx + 23];
double *ddgt2311 = box->v[index_ADMvars_ddgxxxx + 24];
double *ddgt2312 = box->v[index_ADMvars_ddgxxxx + 25];
double *ddgt2313 = box->v[index_ADMvars_ddgxxxx + 26];
double *ddgt2322 = box->v[index_ADMvars_ddgxxxx + 27];
double *ddgt2323 = box->v[index_ADMvars_ddgxxxx + 28];
double *ddgt2333 = box->v[index_ADMvars_ddgxxxx + 29];
double *ddgt3311 = box->v[index_ADMvars_ddgxxxx + 30];
double *ddgt3312 = box->v[index_ADMvars_ddgxxxx + 31];
double *ddgt3313 = box->v[index_ADMvars_ddgxxxx + 32];
double *ddgt3322 = box->v[index_ADMvars_ddgxxxx + 33];
double *ddgt3323 = box->v[index_ADMvars_ddgxxxx + 34];
double *ddgt3333 = box->v[index_ADMvars_ddgxxxx + 35];
int index_ADMvars_dKxxx = Ind("ADMvars_dKxxx");
double *dAt111 = box->v[index_ADMvars_dKxxx + 0];
double *dAt112 = box->v[index_ADMvars_dKxxx + 1];
double *dAt113 = box->v[index_ADMvars_dKxxx + 2];
double *dAt121 = box->v[index_ADMvars_dKxxx + 3];
double *dAt122 = box->v[index_ADMvars_dKxxx + 4];
double *dAt123 = box->v[index_ADMvars_dKxxx + 5];
double *dAt131 = box->v[index_ADMvars_dKxxx + 6];
double *dAt132 = box->v[index_ADMvars_dKxxx + 7];
double *dAt133 = box->v[index_ADMvars_dKxxx + 8];
double *dAt221 = box->v[index_ADMvars_dKxxx + 9];
double *dAt222 = box->v[index_ADMvars_dKxxx + 10];
double *dAt223 = box->v[index_ADMvars_dKxxx + 11];
double *dAt231 = box->v[index_ADMvars_dKxxx + 12];
double *dAt232 = box->v[index_ADMvars_dKxxx + 13];
double *dAt233 = box->v[index_ADMvars_dKxxx + 14];
double *dAt331 = box->v[index_ADMvars_dKxxx + 15];
double *dAt332 = box->v[index_ADMvars_dKxxx + 16];
double *dAt333 = box->v[index_ADMvars_dKxxx + 17];
int index_BSSN_dphix = Ind("BSSN_dphix");
double *dphi1 = box->v[index_BSSN_dphix + 0];
double *dphi2 = box->v[index_BSSN_dphix + 1];
double *dphi3 = box->v[index_BSSN_dphix + 2];
int index_BSSN_ddphixx = Ind("BSSN_ddphixx");
double *ddphi11 = box->v[index_BSSN_ddphixx + 0];
double *ddphi12 = box->v[index_BSSN_ddphixx + 1];
double *ddphi13 = box->v[index_BSSN_ddphixx + 2];
double *ddphi22 = box->v[index_BSSN_ddphixx + 3];
double *ddphi23 = box->v[index_BSSN_ddphixx + 4];
double *ddphi33 = box->v[index_BSSN_ddphixx + 5];
int index_BSSN_dGxx = Ind("BSSN_dGxx");
double *dGt11 = box->v[index_BSSN_dGxx + 0];
double *dGt12 = box->v[index_BSSN_dGxx + 1];
double *dGt13 = box->v[index_BSSN_dGxx + 2];
double *dGt21 = box->v[index_BSSN_dGxx + 3];
double *dGt22 = box->v[index_BSSN_dGxx + 4];
double *dGt23 = box->v[index_BSSN_dGxx + 5];
double *dGt31 = box->v[index_BSSN_dGxx + 6];
double *dGt32 = box->v[index_BSSN_dGxx + 7];
double *dGt33 = box->v[index_BSSN_dGxx + 8];
int index_BSSN_dKx = Ind("BSSN_dKx");
double *dK1 = box->v[index_BSSN_dKx + 0];
double *dK2 = box->v[index_BSSN_dKx + 1];
double *dK3 = box->v[index_BSSN_dKx + 2];
int index_BSSN_dalpx = Ind("BSSN_dalpx");
double *dalp1 = box->v[index_BSSN_dalpx + 0];
double *dalp2 = box->v[index_BSSN_dalpx + 1];
double *dalp3 = box->v[index_BSSN_dalpx + 2];
int index_BSSN_ddalpxx = Ind("BSSN_ddalpxx");
double *ddalp11 = box->v[index_BSSN_ddalpxx + 0];
double *ddalp12 = box->v[index_BSSN_ddalpxx + 1];
double *ddalp13 = box->v[index_BSSN_ddalpxx + 2];
double *ddalp22 = box->v[index_BSSN_ddalpxx + 3];
double *ddalp23 = box->v[index_BSSN_ddalpxx + 4];
double *ddalp33 = box->v[index_BSSN_ddalpxx + 5];
int index_BSSN_dbetaxx = Ind("BSSN_dbetaxx");
double *dbeta11 = box->v[index_BSSN_dbetaxx + 0];
double *dbeta12 = box->v[index_BSSN_dbetaxx + 1];
double *dbeta13 = box->v[index_BSSN_dbetaxx + 2];
double *dbeta21 = box->v[index_BSSN_dbetaxx + 3];
double *dbeta22 = box->v[index_BSSN_dbetaxx + 4];
double *dbeta23 = box->v[index_BSSN_dbetaxx + 5];
double *dbeta31 = box->v[index_BSSN_dbetaxx + 6];
double *dbeta32 = box->v[index_BSSN_dbetaxx + 7];
double *dbeta33 = box->v[index_BSSN_dbetaxx + 8];
int index_BSSN_ddbetaxxx = Ind("BSSN_ddbetaxxx");
double *ddbeta111 = box->v[index_BSSN_ddbetaxxx + 0];
double *ddbeta112 = box->v[index_BSSN_ddbetaxxx + 1];
double *ddbeta113 = box->v[index_BSSN_ddbetaxxx + 2];
double *ddbeta122 = box->v[index_BSSN_ddbetaxxx + 3];
double *ddbeta123 = box->v[index_BSSN_ddbetaxxx + 4];
double *ddbeta133 = box->v[index_BSSN_ddbetaxxx + 5];
double *ddbeta211 = box->v[index_BSSN_ddbetaxxx + 6];
double *ddbeta212 = box->v[index_BSSN_ddbetaxxx + 7];
double *ddbeta213 = box->v[index_BSSN_ddbetaxxx + 8];
double *ddbeta222 = box->v[index_BSSN_ddbetaxxx + 9];
double *ddbeta223 = box->v[index_BSSN_ddbetaxxx + 10];
double *ddbeta233 = box->v[index_BSSN_ddbetaxxx + 11];
double *ddbeta311 = box->v[index_BSSN_ddbetaxxx + 12];
double *ddbeta312 = box->v[index_BSSN_ddbetaxxx + 13];
double *ddbeta313 = box->v[index_BSSN_ddbetaxxx + 14];
double *ddbeta322 = box->v[index_BSSN_ddbetaxxx + 15];
double *ddbeta323 = box->v[index_BSSN_ddbetaxxx + 16];
double *ddbeta333 = box->v[index_BSSN_ddbetaxxx + 17];


double AA;
double aux;
double betadf;
double betadK;
double betaF;
double detginv;
double detnginv;
double divbeta;
double E6alphaDensityWeightf;
double f;
double lieK;
double liephi;
double logpsi;
double nf;
double ooE6alphaDensityWeightf;
double psim4;
double R;
double ralpha;
double ralpha0;
double ralphaDensity;
double rK;
double rphi;
double totdivbeta;
double traceA;
double trcdda;
double trcddf;
double AA11;
double AA12;
double AA13;
double AA22;
double AA23;
double AA33;
double Ainv11;
double Ainv12;
double Ainv13;
double Ainv22;
double Ainv23;
double Ainv33;
double betadA11;
double betadA12;
double betadA13;
double betadA21;
double betadA22;
double betadA23;
double betadA31;
double betadA32;
double betadA33;
double betadg11;
double betadg12;
double betadg13;
double betadg21;
double betadg22;
double betadg23;
double betadg31;
double betadg32;
double betadg33;
double betadG1;
double betadG2;
double betadG3;
double cdda11;
double cdda12;
double cdda13;
double cdda22;
double cdda23;
double cdda33;
double cddf11;
double cddf12;
double cddf13;
double cddf22;
double cddf23;
double cddf33;
double da1;
double da2;
double da3;
double dA111;
double dA112;
double dA113;
double dA122;
double dA123;
double dA133;
double dA211;
double dA212;
double dA213;
double dA222;
double dA223;
double dA233;
double dA311;
double dA312;
double dA313;
double dA322;
double dA323;
double dA333;
double db11;
double db12;
double db13;
double db21;
double db22;
double db23;
double db31;
double db32;
double db33;
double dda11;
double dda12;
double dda13;
double dda22;
double dda23;
double dda33;
double ddb111;
double ddb112;
double ddb113;
double ddb121;
double ddb122;
double ddb123;
double ddb131;
double ddb132;
double ddb133;
double ddb221;
double ddb222;
double ddb223;
double ddb231;
double ddb232;
double ddb233;
double ddb331;
double ddb332;
double ddb333;
double ddf11;
double ddf12;
double ddf13;
double ddf22;
double ddf23;
double ddf33;
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
double delG11;
double delG12;
double delG13;
double delG21;
double delG22;
double delG23;
double delG31;
double delG32;
double delG33;
double df1;
double df2;
double df3;
double dginv111;
double dginv112;
double dginv113;
double dginv122;
double dginv123;
double dginv133;
double dginv211;
double dginv212;
double dginv213;
double dginv222;
double dginv223;
double dginv233;
double dginv311;
double dginv312;
double dginv313;
double dginv322;
double dginv323;
double dginv333;
double divAinv1;
double divAinv2;
double divAinv3;
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
double Gfromg1;
double Gfromg2;
double Gfromg3;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double lieA11;
double lieA12;
double lieA13;
double lieA22;
double lieA23;
double lieA33;
double lieg11;
double lieg12;
double lieg13;
double lieg22;
double lieg23;
double lieg33;
double ootddivbeta1;
double ootddivbeta2;
double ootddivbeta3;
double pseudolieG1;
double pseudolieG2;
double pseudolieG3;
double R11;
double R12;
double R13;
double R22;
double R23;
double R33;
double rA11;
double rA12;
double rA13;
double rA22;
double rA23;
double rA33;
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
double rG1;
double rG2;
double rG3;
double Rphi11;
double Rphi12;
double Rphi13;
double Rphi22;
double Rphi23;
double Rphi33;



/* Jetzt geht's los! */

FirstAndSecondDerivsOf_Sab(box, index_g11,                     Ind("ADMvars_dgxxx"), Ind("ADMvars_ddgxxxx")); 


FirstDerivsOf_Sab(box, index_A11, Ind("ADMvars_dKxxx")); 


FirstAndSecondDerivsOf_S(box, index_phi,                     Ind("BSSN_dphix"), Ind("BSSN_ddphixx")); 


FirstDerivsOf_Sa(box, index_G1, Ind("BSSN_dGxx")); 


FirstDerivsOf_S(box, index_K, Ind("BSSN_dKx")); 


FirstAndSecondDerivsOf_S(box, index_alpha,                     Ind("BSSN_dalpx"), Ind("BSSN_ddalpxx")); 


FirstAndSecondDerivsOf_Sa(box, index_beta1,                     Ind("BSSN_dbetaxx"), Ind("BSSN_ddbetaxxx")); 


forallpoints(box, ijk) { 

df1
=
dphi1[ijk]
;

df2
=
dphi2[ijk]
;

df3
=
dphi3[ijk]
;

ddf11
=
ddphi11[ijk]
;

ddf12
=
ddphi12[ijk]
;

ddf13
=
ddphi13[ijk]
;

ddf22
=
ddphi22[ijk]
;

ddf23
=
ddphi23[ijk]
;

ddf33
=
ddphi33[ijk]
;

da1
=
dalp1[ijk]
;

da2
=
dalp2[ijk]
;

da3
=
dalp3[ijk]
;

dda11
=
ddalp11[ijk]
;

dda12
=
ddalp12[ijk]
;

dda13
=
ddalp13[ijk]
;

dda22
=
ddalp22[ijk]
;

dda23
=
ddalp23[ijk]
;

dda33
=
ddalp33[ijk]
;

db11
=
dbeta11[ijk]
;

db12
=
dbeta21[ijk]
;

db13
=
dbeta31[ijk]
;

db21
=
dbeta12[ijk]
;

db22
=
dbeta22[ijk]
;

db23
=
dbeta32[ijk]
;

db31
=
dbeta13[ijk]
;

db32
=
dbeta23[ijk]
;

db33
=
dbeta33[ijk]
;

ddb111
=
ddbeta111[ijk]
;

ddb112
=
ddbeta211[ijk]
;

ddb113
=
ddbeta311[ijk]
;

ddb121
=
ddbeta112[ijk]
;

ddb122
=
ddbeta212[ijk]
;

ddb123
=
ddbeta312[ijk]
;

ddb131
=
ddbeta113[ijk]
;

ddb132
=
ddbeta213[ijk]
;

ddb133
=
ddbeta313[ijk]
;

ddb221
=
ddbeta122[ijk]
;

ddb222
=
ddbeta222[ijk]
;

ddb223
=
ddbeta322[ijk]
;

ddb231
=
ddbeta123[ijk]
;

ddb232
=
ddbeta223[ijk]
;

ddb233
=
ddbeta323[ijk]
;

ddb331
=
ddbeta133[ijk]
;

ddb332
=
ddbeta233[ijk]
;

ddb333
=
ddbeta333[ijk]
;

delg111
=
dgt111[ijk]
;

delg112
=
dgt121[ijk]
;

delg113
=
dgt131[ijk]
;

delg122
=
dgt221[ijk]
;

delg123
=
dgt231[ijk]
;

delg133
=
dgt331[ijk]
;

delg211
=
dgt112[ijk]
;

delg212
=
dgt122[ijk]
;

delg213
=
dgt132[ijk]
;

delg222
=
dgt222[ijk]
;

delg223
=
dgt232[ijk]
;

delg233
=
dgt332[ijk]
;

delg311
=
dgt113[ijk]
;

delg312
=
dgt123[ijk]
;

delg313
=
dgt133[ijk]
;

delg322
=
dgt223[ijk]
;

delg323
=
dgt233[ijk]
;

delg333
=
dgt333[ijk]
;

deldelg1111
=
ddgt1111[ijk]
;

deldelg1112
=
ddgt1211[ijk]
;

deldelg1113
=
ddgt1311[ijk]
;

deldelg1122
=
ddgt2211[ijk]
;

deldelg1123
=
ddgt2311[ijk]
;

deldelg1133
=
ddgt3311[ijk]
;

deldelg1211
=
ddgt1112[ijk]
;

deldelg1212
=
ddgt1212[ijk]
;

deldelg1213
=
ddgt1312[ijk]
;

deldelg1222
=
ddgt2212[ijk]
;

deldelg1223
=
ddgt2312[ijk]
;

deldelg1233
=
ddgt3312[ijk]
;

deldelg1311
=
ddgt1113[ijk]
;

deldelg1312
=
ddgt1213[ijk]
;

deldelg1313
=
ddgt1313[ijk]
;

deldelg1322
=
ddgt2213[ijk]
;

deldelg1323
=
ddgt2313[ijk]
;

deldelg1333
=
ddgt3313[ijk]
;

deldelg2211
=
ddgt1122[ijk]
;

deldelg2212
=
ddgt1222[ijk]
;

deldelg2213
=
ddgt1322[ijk]
;

deldelg2222
=
ddgt2222[ijk]
;

deldelg2223
=
ddgt2322[ijk]
;

deldelg2233
=
ddgt3322[ijk]
;

deldelg2311
=
ddgt1123[ijk]
;

deldelg2312
=
ddgt1223[ijk]
;

deldelg2313
=
ddgt1323[ijk]
;

deldelg2322
=
ddgt2223[ijk]
;

deldelg2323
=
ddgt2323[ijk]
;

deldelg2333
=
ddgt3323[ijk]
;

deldelg3311
=
ddgt1133[ijk]
;

deldelg3312
=
ddgt1233[ijk]
;

deldelg3313
=
ddgt1333[ijk]
;

deldelg3322
=
ddgt2233[ijk]
;

deldelg3323
=
ddgt2333[ijk]
;

deldelg3333
=
ddgt3333[ijk]
;

delG11
=
dGt11[ijk]
;

delG12
=
dGt21[ijk]
;

delG13
=
dGt31[ijk]
;

delG21
=
dGt12[ijk]
;

delG22
=
dGt22[ijk]
;

delG23
=
dGt32[ijk]
;

delG31
=
dGt13[ijk]
;

delG32
=
dGt23[ijk]
;

delG33
=
dGt33[ijk]
;

dA111
=
dAt111[ijk]
;

dA112
=
dAt121[ijk]
;

dA113
=
dAt131[ijk]
;

dA122
=
dAt221[ijk]
;

dA123
=
dAt231[ijk]
;

dA133
=
dAt331[ijk]
;

dA211
=
dAt112[ijk]
;

dA212
=
dAt122[ijk]
;

dA213
=
dAt132[ijk]
;

dA222
=
dAt222[ijk]
;

dA223
=
dAt232[ijk]
;

dA233
=
dAt332[ijk]
;

dA311
=
dAt113[ijk]
;

dA312
=
dAt123[ijk]
;

dA313
=
dAt133[ijk]
;

dA322
=
dAt223[ijk]
;

dA323
=
dAt233[ijk]
;

dA333
=
dAt333[ijk]
;

betadf
=
beta1[ijk]*dphi1[ijk] + beta2[ijk]*dphi2[ijk] + beta3[ijk]*dphi3[ijk]
;

betadg11
=
beta1[ijk]*dgt111[ijk] + beta2[ijk]*dgt112[ijk] + beta3[ijk]*dgt113[ijk]
;

betadg12
=
beta1[ijk]*dgt121[ijk] + beta2[ijk]*dgt122[ijk] + beta3[ijk]*dgt123[ijk]
;

betadg13
=
beta1[ijk]*dgt131[ijk] + beta2[ijk]*dgt132[ijk] + beta3[ijk]*dgt133[ijk]
;

betadg21
=
beta1[ijk]*dgt121[ijk] + beta2[ijk]*dgt122[ijk] + beta3[ijk]*dgt123[ijk]
;

betadg22
=
beta1[ijk]*dgt221[ijk] + beta2[ijk]*dgt222[ijk] + beta3[ijk]*dgt223[ijk]
;

betadg23
=
beta1[ijk]*dgt231[ijk] + beta2[ijk]*dgt232[ijk] + beta3[ijk]*dgt233[ijk]
;

betadg31
=
beta1[ijk]*dgt131[ijk] + beta2[ijk]*dgt132[ijk] + beta3[ijk]*dgt133[ijk]
;

betadg32
=
beta1[ijk]*dgt231[ijk] + beta2[ijk]*dgt232[ijk] + beta3[ijk]*dgt233[ijk]
;

betadg33
=
beta1[ijk]*dgt331[ijk] + beta2[ijk]*dgt332[ijk] + beta3[ijk]*dgt333[ijk]
;

betadA11
=
beta1[ijk]*dAt111[ijk] + beta2[ijk]*dAt112[ijk] + beta3[ijk]*dAt113[ijk]
;

betadA12
=
beta1[ijk]*dAt121[ijk] + beta2[ijk]*dAt122[ijk] + beta3[ijk]*dAt123[ijk]
;

betadA13
=
beta1[ijk]*dAt131[ijk] + beta2[ijk]*dAt132[ijk] + beta3[ijk]*dAt133[ijk]
;

betadA21
=
beta1[ijk]*dAt121[ijk] + beta2[ijk]*dAt122[ijk] + beta3[ijk]*dAt123[ijk]
;

betadA22
=
beta1[ijk]*dAt221[ijk] + beta2[ijk]*dAt222[ijk] + beta3[ijk]*dAt223[ijk]
;

betadA23
=
beta1[ijk]*dAt231[ijk] + beta2[ijk]*dAt232[ijk] + beta3[ijk]*dAt233[ijk]
;

betadA31
=
beta1[ijk]*dAt131[ijk] + beta2[ijk]*dAt132[ijk] + beta3[ijk]*dAt133[ijk]
;

betadA32
=
beta1[ijk]*dAt231[ijk] + beta2[ijk]*dAt232[ijk] + beta3[ijk]*dAt233[ijk]
;

betadA33
=
beta1[ijk]*dAt331[ijk] + beta2[ijk]*dAt332[ijk] + beta3[ijk]*dAt333[ijk]
;

betadK
=
beta1[ijk]*dK1[ijk] + beta2[ijk]*dK2[ijk] + beta3[ijk]*dK3[ijk]
;

betadG1
=
beta1[ijk]*dGt11[ijk] + beta2[ijk]*dGt12[ijk] + beta3[ijk]*dGt13[ijk]
;

betadG2
=
beta1[ijk]*dGt21[ijk] + beta2[ijk]*dGt22[ijk] + beta3[ijk]*dGt23[ijk]
;

betadG3
=
beta1[ijk]*dGt31[ijk] + beta2[ijk]*dGt32[ijk] + beta3[ijk]*dGt33[ijk]
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

dginv111
=
-2.*(delg123*ginv12*ginv13 + ginv11*(delg112*ginv12 + delg113*ginv13)) - 
  delg111*pow2(ginv11) - delg122*pow2(ginv12) - delg133*pow2(ginv13)
;

dginv112
=
-(ginv11*(delg111*ginv12 + delg112*ginv22 + delg113*ginv23)) - 
  ginv12*(delg113*ginv13 + delg122*ginv22 + delg123*ginv23) - 
  ginv13*(delg123*ginv22 + delg133*ginv23) - delg112*pow2(ginv12)
;

dginv113
=
-(ginv11*(delg111*ginv13 + delg112*ginv23 + delg113*ginv33)) - 
  ginv12*(delg112*ginv13 + delg122*ginv23 + delg123*ginv33) - 
  ginv13*(delg123*ginv23 + delg133*ginv33) - delg113*pow2(ginv13)
;

dginv122
=
-2.*(delg123*ginv22*ginv23 + ginv12*(delg112*ginv22 + delg113*ginv23)) - 
  delg111*pow2(ginv12) - delg122*pow2(ginv22) - delg133*pow2(ginv23)
;

dginv123
=
-(ginv13*(delg112*ginv22 + delg113*ginv23)) - delg133*ginv23*ginv33 - 
  ginv12*(delg111*ginv13 + delg112*ginv23 + delg113*ginv33) - 
  ginv22*(delg122*ginv23 + delg123*ginv33) - delg123*pow2(ginv23)
;

dginv133
=
-2.*(delg123*ginv23*ginv33 + ginv13*(delg112*ginv23 + delg113*ginv33)) - 
  delg111*pow2(ginv13) - delg122*pow2(ginv23) - delg133*pow2(ginv33)
;

dginv211
=
-2.*(delg223*ginv12*ginv13 + ginv11*(delg212*ginv12 + delg213*ginv13)) - 
  delg211*pow2(ginv11) - delg222*pow2(ginv12) - delg233*pow2(ginv13)
;

dginv212
=
-(ginv11*(delg211*ginv12 + delg212*ginv22 + delg213*ginv23)) - 
  ginv12*(delg213*ginv13 + delg222*ginv22 + delg223*ginv23) - 
  ginv13*(delg223*ginv22 + delg233*ginv23) - delg212*pow2(ginv12)
;

dginv213
=
-(ginv11*(delg211*ginv13 + delg212*ginv23 + delg213*ginv33)) - 
  ginv12*(delg212*ginv13 + delg222*ginv23 + delg223*ginv33) - 
  ginv13*(delg223*ginv23 + delg233*ginv33) - delg213*pow2(ginv13)
;

dginv222
=
-2.*(delg223*ginv22*ginv23 + ginv12*(delg212*ginv22 + delg213*ginv23)) - 
  delg211*pow2(ginv12) - delg222*pow2(ginv22) - delg233*pow2(ginv23)
;

dginv223
=
-(ginv13*(delg212*ginv22 + delg213*ginv23)) - delg233*ginv23*ginv33 - 
  ginv12*(delg211*ginv13 + delg212*ginv23 + delg213*ginv33) - 
  ginv22*(delg222*ginv23 + delg223*ginv33) - delg223*pow2(ginv23)
;

dginv233
=
-2.*(delg223*ginv23*ginv33 + ginv13*(delg212*ginv23 + delg213*ginv33)) - 
  delg211*pow2(ginv13) - delg222*pow2(ginv23) - delg233*pow2(ginv33)
;

dginv311
=
-2.*(delg323*ginv12*ginv13 + ginv11*(delg312*ginv12 + delg313*ginv13)) - 
  delg311*pow2(ginv11) - delg322*pow2(ginv12) - delg333*pow2(ginv13)
;

dginv312
=
-(ginv11*(delg311*ginv12 + delg312*ginv22 + delg313*ginv23)) - 
  ginv12*(delg313*ginv13 + delg322*ginv22 + delg323*ginv23) - 
  ginv13*(delg323*ginv22 + delg333*ginv23) - delg312*pow2(ginv12)
;

dginv313
=
-(ginv11*(delg311*ginv13 + delg312*ginv23 + delg313*ginv33)) - 
  ginv12*(delg312*ginv13 + delg322*ginv23 + delg323*ginv33) - 
  ginv13*(delg323*ginv23 + delg333*ginv33) - delg313*pow2(ginv13)
;

dginv322
=
-2.*(delg323*ginv22*ginv23 + ginv12*(delg312*ginv22 + delg313*ginv23)) - 
  delg311*pow2(ginv12) - delg322*pow2(ginv22) - delg333*pow2(ginv23)
;

dginv323
=
-(ginv13*(delg312*ginv22 + delg313*ginv23)) - delg333*ginv23*ginv33 - 
  ginv12*(delg311*ginv13 + delg312*ginv23 + delg313*ginv33) - 
  ginv22*(delg322*ginv23 + delg323*ginv33) - delg323*pow2(ginv23)
;

dginv333
=
-2.*(delg323*ginv23*ginv33 + ginv13*(delg312*ginv23 + delg313*ginv33)) - 
  delg311*pow2(ginv13) - delg322*pow2(ginv23) - delg333*pow2(ginv33)
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

Gfromg1
=
gamma111*ginv11 + (gamma112 + gamma121)*ginv12 + 
  (gamma113 + gamma131)*ginv13 + gamma122*ginv22 + 
  (gamma123 + gamma132)*ginv23 + gamma133*ginv33
;

Gfromg2
=
gamma211*ginv11 + (gamma212 + gamma221)*ginv12 + 
  (gamma213 + gamma231)*ginv13 + gamma222*ginv22 + 
  (gamma223 + gamma232)*ginv23 + gamma233*ginv33
;

Gfromg3
=
gamma311*ginv11 + (gamma312 + gamma321)*ginv12 + 
  (gamma313 + gamma331)*ginv13 + gamma322*ginv22 + 
  (gamma323 + gamma332)*ginv23 + gamma333*ginv33
;

R11
=
gammado111*Gfromg1 + gammado112*Gfromg2 + gammado113*Gfromg3 + 
  (-0.5*deldelg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado121 + gamma311*gammado131) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-deldelg1211 + (gamma112 + 2.*gamma121)*gammado111 + 
     gamma111*(2.*gammado112 + gammado121) + 
     2.*(gamma221*gammado121 + gamma211*gammado122 + 
        gamma321*gammado131 + gamma311*gammado132) + gamma212*gammado211 + 
     gamma211*gammado221 + gamma312*gammado311 + gamma311*gammado321)*ginv12 \
+ (-deldelg1311 + (gamma113 + 2.*gamma131)*gammado111 + 
     gamma111*gammado131 + 2.*
      (gamma111*gammado113 + gamma231*gammado121 + gamma211*gammado123 + 
        gamma331*gammado131 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado231 + gamma313*gammado311 + gamma311*gammado331)*ginv13 \
+ (-0.5*deldelg2211 + gamma112*gammado121 + 
     2.*(gamma121*gammado112 + gamma221*gammado122 + 
        gamma321*gammado132) + gamma212*gammado221 + gamma312*gammado321)*
   ginv22 + (-deldelg2311 + gamma113*gammado121 + gamma112*gammado131 + 
     2.*(gamma131*gammado112 + gamma121*gammado113 + 
        gamma231*gammado122 + gamma221*gammado123 + gamma331*gammado132 + 
        gamma321*gammado133) + gamma213*gammado221 + gamma212*gammado231 + 
     gamma313*gammado321 + gamma312*gammado331)*ginv23 + 
  (-0.5*deldelg3311 + gamma113*gammado131 + 
     2.*(gamma131*gammado113 + gamma231*gammado123 + 
        gamma331*gammado133) + gamma213*gammado231 + gamma313*gammado331)*
   ginv33 + delG11*g11[ijk] + delG12*g12[ijk] + delG13*g13[ijk]
;

R12
=
(-0.5*deldelg1112 + (gamma112 + gamma121)*gammado111 + 
     gamma212*gammado121 + gamma312*gammado131 + 
     (gamma111 + gamma221)*gammado211 + gamma211*gammado221 + 
     gamma311*gammado231 + gamma321*gammado311)*ginv11 + 
  (-deldelg1212 + gamma112*gammado112 + gamma212*gammado122 + 
     gamma312*gammado132 + (gamma121 + gamma222)*
      (gammado121 + gammado211) + gamma111*gammado212 + 
     2.*(gamma122*gammado111 + gamma221*gammado221) + 
     gamma211*gammado222 + gamma311*gammado232 + 
     gamma322*(gammado131 + gammado311) + 
     gamma321*(gammado231 + gammado321))*ginv12 + 
  (-deldelg1312 + (gamma123 + gamma132)*gammado111 + gamma112*gammado113 + 
     gamma232*gammado121 + gamma212*gammado123 + 
     (gamma121 + gamma332)*gammado131 + gamma312*gammado133 + 
     (gamma131 + gamma223)*gammado211 + gamma111*gammado213 + 
     gamma231*gammado221 + gamma211*gammado223 + 
     (gamma221 + gamma331)*gammado231 + gamma311*gammado233 + 
     gamma323*gammado311 + gamma321*gammado331)*ginv13 + 
  (-0.5*deldelg2212 + gamma122*(gammado112 + gammado121) + 
     gamma121*gammado212 + gamma222*(gammado122 + gammado221) + 
     gamma221*gammado222 + gamma321*gammado232 + 
     gamma322*(gammado132 + gammado321))*ginv22 + 
  (-deldelg2312 + gamma132*gammado112 + gamma123*gammado121 + 
     gamma232*gammado122 + gamma122*(gammado113 + gammado131) + 
     gamma332*gammado132 + gamma131*gammado212 + gamma121*gammado213 + 
     gamma223*gammado221 + gamma231*gammado222 + gamma221*gammado223 + 
     gamma222*(gammado123 + gammado231) + gamma331*gammado232 + 
     gamma321*gammado233 + gamma323*gammado321 + 
     gamma322*(gammado133 + gammado331))*ginv23 + 
  (-0.5*deldelg3312 + gamma132*gammado113 + gamma232*gammado123 + 
     gamma123*gammado131 + gamma332*gammado133 + gamma131*gammado213 + 
     gamma231*gammado223 + gamma223*gammado231 + gamma331*gammado233 + 
     gamma323*gammado331)*ginv33 + 
  0.5*((gammado121 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3 + 
     delG21*g11[ijk] + (delG11 + delG22)*g12[ijk] + delG23*g13[ijk] + 
     delG12*g22[ijk] + delG13*g23[ijk])
;

R12
=
(-0.5*deldelg1112 + gamma112*gammado111 + gamma212*gammado121 + 
     gamma312*gammado131 + gamma111*(gammado112 + gammado211) + 
     gamma211*(gammado212 + gammado221) + 
     gamma311*(gammado231 + gammado312))*ginv11 + 
  (-deldelg1212 + gamma122*gammado111 + gamma222*gammado121 + 
     gamma322*gammado131 + gamma121*gammado211 + 
     (gamma111 + gamma212)*(gammado122 + gammado212) + 
     gamma221*gammado221 + 2.*
      (gamma112*gammado112 + gamma211*gammado222) + gamma321*gammado231 + 
     gamma312*(gammado132 + gammado312) + 
     gamma311*(gammado232 + gammado322))*ginv12 + 
  (-deldelg1312 + gamma132*gammado111 + gamma113*gammado112 + 
     gamma112*gammado113 + gamma232*gammado121 + gamma212*gammado123 + 
     gamma332*gammado131 + gamma312*gammado133 + gamma131*gammado211 + 
     gamma213*gammado212 + gamma111*(gammado132 + gammado213) + 
     gamma231*gammado221 + gamma331*gammado231 + 
     gamma211*(gammado223 + gammado232) + gamma313*gammado312 + 
     gamma311*(gammado233 + gammado332))*ginv13 + 
  (-0.5*deldelg2212 + gamma122*gammado112 + 
     (gamma112 + gamma222)*gammado122 + gamma322*gammado132 + 
     gamma121*gammado212 + (gamma212 + gamma221)*gammado222 + 
     gamma321*gammado232 + gamma312*gammado322)*ginv22 + 
  (-deldelg2312 + gamma132*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma232)*gammado122 + gamma222*gammado123 + 
     (gamma112 + gamma332)*gammado132 + gamma322*gammado133 + 
     gamma131*gammado212 + gamma121*gammado213 + 
     (gamma213 + gamma231)*gammado222 + gamma221*gammado223 + 
     (gamma212 + gamma331)*gammado232 + gamma321*gammado233 + 
     gamma313*gammado322 + gamma312*gammado332)*ginv23 + 
  (-0.5*deldelg3312 + gamma132*gammado113 + gamma232*gammado123 + 
     gamma113*gammado132 + gamma332*gammado133 + gamma131*gammado213 + 
     gamma231*gammado223 + gamma213*gammado232 + gamma331*gammado233 + 
     gamma313*gammado332)*ginv33 + 
  0.5*((gammado121 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3 + 
     delG21*g11[ijk] + (delG11 + delG22)*g12[ijk] + delG23*g13[ijk] + 
     delG12*g22[ijk] + delG13*g23[ijk])
;

R13
=
(-0.5*deldelg1113 + (gamma113 + gamma131)*gammado111 + 
     gamma213*gammado121 + gamma313*gammado131 + gamma231*gammado211 + 
     (gamma111 + gamma331)*gammado311 + gamma211*gammado321 + 
     gamma311*gammado331)*ginv11 + 
  (-deldelg1213 + (gamma123 + gamma132)*gammado111 + gamma113*gammado112 + 
     (gamma131 + gamma223)*gammado121 + gamma213*gammado122 + 
     gamma323*gammado131 + gamma313*gammado132 + gamma232*gammado211 + 
     gamma231*gammado221 + (gamma121 + gamma332)*gammado311 + 
     gamma111*gammado312 + (gamma221 + gamma331)*gammado321 + 
     gamma211*gammado322 + gamma321*gammado331 + gamma311*gammado332)*ginv12 \
+ (-deldelg1313 + gamma113*gammado113 + gamma213*gammado123 + 
     gamma313*gammado133 + gamma233*(gammado121 + gammado211) + 
     (gamma131 + gamma333)*(gammado131 + gammado311) + 
     gamma111*gammado313 + gamma231*(gammado231 + gammado321) + 
     gamma211*gammado323 + 2.*
      (gamma133*gammado111 + gamma331*gammado331) + gamma311*gammado333)*
   ginv13 + (-0.5*deldelg2213 + gamma123*gammado112 + 
     gamma132*gammado121 + gamma223*gammado122 + gamma323*gammado132 + 
     gamma232*gammado221 + gamma121*gammado312 + gamma332*gammado321 + 
     gamma221*gammado322 + gamma321*gammado332)*ginv22 + 
  (-deldelg2313 + gamma123*gammado113 + 
     gamma133*(gammado112 + gammado121) + gamma223*gammado123 + 
     gamma132*gammado131 + gamma323*gammado133 + 
     gamma233*(gammado122 + gammado221) + gamma232*gammado231 + 
     gamma131*gammado312 + gamma121*gammado313 + 
     gamma333*(gammado132 + gammado321) + gamma231*gammado322 + 
     gamma221*gammado323 + gamma332*gammado331 + gamma331*gammado332 + 
     gamma321*gammado333)*ginv23 + 
  (-0.5*deldelg3313 + gamma133*(gammado113 + gammado131) + 
     gamma233*(gammado123 + gammado231) + gamma131*gammado313 + 
     gamma231*gammado323 + gamma333*(gammado133 + gammado331) + 
     gamma331*gammado333)*ginv33 + 
  0.5*((gammado131 + gammado311)*Gfromg1 + 
     (gammado132 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3 + 
     delG31*g11[ijk] + delG32*g12[ijk] + (delG11 + delG33)*g13[ijk] + 
     delG12*g23[ijk] + delG13*g33[ijk])
;

R13
=
(-0.5*deldelg1113 + gamma113*gammado111 + gamma213*gammado121 + 
     gamma313*gammado131 + gamma111*(gammado113 + gammado311) + 
     gamma211*(gammado213 + gammado321) + 
     gamma311*(gammado313 + gammado331))*ginv11 + 
  (-deldelg1213 + gamma123*gammado111 + gamma113*gammado112 + 
     gamma112*gammado113 + gamma223*gammado121 + gamma213*gammado122 + 
     gamma323*gammado131 + gamma313*gammado132 + gamma212*gammado213 + 
     gamma121*gammado311 + gamma111*(gammado123 + gammado312) + 
     gamma312*gammado313 + gamma221*gammado321 + 
     gamma211*(gammado223 + gammado322) + gamma321*gammado331 + 
     gamma311*(gammado323 + gammado332))*ginv12 + 
  (-deldelg1313 + gamma133*gammado111 + gamma233*gammado121 + 
     gamma333*gammado131 + gamma213*(gammado123 + gammado213) + 
     gamma131*gammado311 + (gamma111 + gamma313)*
      (gammado133 + gammado313) + gamma231*gammado321 + 
     gamma211*(gammado233 + gammado323) + gamma331*gammado331 + 
     2.*(gamma113*gammado113 + gamma311*gammado333))*ginv13 + 
  (-0.5*deldelg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     gamma112*gammado123 + gamma323*gammado132 + gamma212*gammado223 + 
     gamma121*gammado312 + gamma221*gammado322 + gamma312*gammado323 + 
     gamma321*gammado332)*ginv22 + 
  (-deldelg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223)*gammado123 + 
     gamma333*gammado132 + (gamma112 + gamma323)*gammado133 + 
     gamma213*gammado223 + gamma212*gammado233 + gamma131*gammado312 + 
     gamma121*gammado313 + gamma231*gammado322 + 
     (gamma221 + gamma313)*gammado323 + gamma331*gammado332 + 
     (gamma312 + gamma321)*gammado333)*ginv23 + 
  (-0.5*deldelg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma213*gammado233 + 
     gamma131*gammado313 + gamma231*gammado323 + 
     (gamma313 + gamma331)*gammado333)*ginv33 + 
  0.5*((gammado131 + gammado311)*Gfromg1 + 
     (gammado132 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3 + 
     delG31*g11[ijk] + delG32*g12[ijk] + (delG11 + delG33)*g13[ijk] + 
     delG12*g23[ijk] + delG13*g33[ijk])
;

R22
=
gammado221*Gfromg1 + gammado222*Gfromg2 + gammado223*Gfromg3 + 
  (-0.5*deldelg1122 + gamma121*gammado112 + gamma221*gammado212 + 
     2.*(gamma112*gammado211 + gamma212*gammado221 + 
        gamma312*gammado231) + gamma321*gammado312)*ginv11 + 
  (-deldelg1222 + gamma121*gammado122 + 
     gamma122*(gammado112 + 2.*gammado211) + 
     (2.*gamma112 + gamma222)*gammado212 + gamma221*gammado222 + 
     2.*(gamma222*gammado221 + gamma212*gammado222 + 
        gamma322*gammado231 + gamma312*gammado232) + gamma322*gammado312 + 
     gamma321*gammado322)*ginv12 + 
  (-deldelg1322 + gamma123*gammado112 + gamma121*gammado132 + 
     gamma223*gammado212 + gamma221*gammado232 + 
     2.*(gamma132*gammado211 + gamma112*gammado213 + 
        gamma232*gammado221 + gamma212*gammado223 + gamma332*gammado231 + 
        gamma312*gammado233) + gamma323*gammado312 + gamma321*gammado332)*
   ginv13 + (-0.5*deldelg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado232 + gammado322))*ginv22 \
+ (-deldelg2322 + gamma123*gammado122 + 
     gamma122*(gammado132 + 2.*gammado213) + gamma223*gammado222 + 
     gamma222*gammado232 + 2.*
      (gamma132*gammado212 + gamma232*gammado222 + gamma222*gammado223 + 
        gamma332*gammado232 + gamma322*gammado233) + gamma323*gammado322 + 
     gamma322*gammado332)*ginv23 + 
  (-0.5*deldelg3322 + gamma123*gammado132 + gamma223*gammado232 + 
     2.*(gamma132*gammado213 + gamma232*gammado223 + 
        gamma332*gammado233) + gamma323*gammado332)*ginv33 + 
  delG21*g12[ijk] + delG22*g22[ijk] + delG23*g23[ijk]
;

R23
=
(-0.5*deldelg1123 + gamma131*gammado112 + gamma113*gammado211 + 
     gamma231*gammado212 + gamma213*gammado221 + gamma313*gammado231 + 
     gamma112*gammado311 + gamma331*gammado312 + gamma212*gammado321 + 
     gamma312*gammado331)*ginv11 + 
  (-deldelg1223 + gamma132*gammado112 + gamma131*gammado122 + 
     gamma123*gammado211 + (gamma113 + gamma232)*gammado212 + 
     gamma223*gammado221 + (gamma213 + gamma231)*gammado222 + 
     gamma323*gammado231 + gamma313*gammado232 + gamma122*gammado311 + 
     (gamma112 + gamma332)*gammado312 + gamma222*gammado321 + 
     (gamma212 + gamma331)*gammado322 + gamma322*gammado331 + 
     gamma312*gammado332)*ginv12 + 
  (-deldelg1323 + gamma131*gammado132 + 
     gamma133*(gammado112 + gammado211) + gamma113*gammado213 + 
     gamma233*(gammado212 + gammado221) + gamma213*gammado223 + 
     gamma231*gammado232 + gamma313*gammado233 + gamma132*gammado311 + 
     gamma333*(gammado231 + gammado312) + gamma112*gammado313 + 
     gamma232*gammado321 + gamma212*gammado323 + gamma332*gammado331 + 
     gamma331*gammado332 + gamma312*gammado333)*ginv13 + 
  (-0.5*deldelg2223 + gamma132*gammado122 + gamma123*gammado212 + 
     (gamma223 + gamma232)*gammado222 + gamma323*gammado232 + 
     gamma122*gammado312 + (gamma222 + gamma332)*gammado322 + 
     gamma322*gammado332)*ginv22 + 
  (-deldelg2323 + gamma133*(gammado122 + gammado212) + 
     gamma123*gammado213 + gamma223*gammado223 + gamma323*gammado233 + 
     gamma132*(gammado132 + gammado312) + gamma122*gammado313 + 
     (gamma232 + gamma333)*(gammado232 + gammado322) + 
     gamma222*gammado323 + 2.*
      (gamma233*gammado222 + gamma332*gammado332) + gamma322*gammado333)*
   ginv23 + (-0.5*deldelg3323 + gamma133*(gammado132 + gammado213) + 
     gamma233*(gammado223 + gammado232) + gamma132*gammado313 + 
     gamma232*gammado323 + gamma333*(gammado233 + gammado332) + 
     gamma332*gammado333)*ginv33 + 
  0.5*((gammado231 + gammado321)*Gfromg1 + 
     (gammado232 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3 + 
     delG31*g12[ijk] + delG21*g13[ijk] + delG32*g22[ijk] + 
     (delG22 + delG33)*g23[ijk] + delG23*g33[ijk])
;

R23
=
(-0.5*deldelg1123 + gamma121*gammado113 + gamma113*gammado211 + 
     gamma221*gammado213 + gamma213*gammado221 + gamma313*gammado231 + 
     gamma112*gammado311 + gamma321*gammado313 + gamma212*gammado321 + 
     gamma312*gammado331)*ginv11 + 
  (-deldelg1223 + gamma121*gammado123 + gamma123*gammado211 + 
     gamma113*gammado212 + gamma223*gammado221 + gamma213*gammado222 + 
     gamma221*gammado223 + gamma323*gammado231 + gamma313*gammado232 + 
     gamma122*(gammado113 + gammado311) + gamma112*gammado312 + 
     gamma222*(gammado213 + gammado321) + gamma212*gammado322 + 
     gamma321*gammado323 + gamma322*(gammado313 + gammado331) + 
     gamma312*gammado332)*ginv12 + 
  (-deldelg1323 + gamma123*gammado113 + gamma121*gammado133 + 
     gamma133*gammado211 + (gamma113 + gamma223)*gammado213 + 
     gamma233*gammado221 + gamma213*gammado223 + gamma333*gammado231 + 
     (gamma221 + gamma313)*gammado233 + gamma132*gammado311 + 
     (gamma112 + gamma323)*gammado313 + gamma232*gammado321 + 
     gamma212*gammado323 + gamma332*gammado331 + 
     (gamma312 + gamma321)*gammado333)*ginv13 + 
  (-0.5*deldelg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     gamma323*gammado232 + gamma122*(gammado123 + gammado312) + 
     gamma222*(gammado223 + gammado322) + 
     gamma322*(gammado323 + gammado332))*ginv22 + 
  (-deldelg2323 + gamma133*gammado212 + 
     gamma123*(gammado123 + gammado213) + gamma233*gammado222 + 
     gamma333*gammado232 + gamma132*gammado312 + 
     gamma122*(gammado133 + gammado313) + gamma232*gammado322 + 
     (gamma222 + gamma323)*(gammado233 + gammado323) + 
     gamma332*gammado332 + 2.*(gamma223*gammado223 + gamma322*gammado333))*
   ginv23 + (-0.5*deldelg3323 + gamma123*gammado133 + 
     gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + gamma132*gammado313 + 
     gamma232*gammado323 + (gamma323 + gamma332)*gammado333)*ginv33 + 
  0.5*((gammado231 + gammado321)*Gfromg1 + 
     (gammado232 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3 + 
     delG31*g12[ijk] + delG21*g13[ijk] + delG32*g22[ijk] + 
     (delG22 + delG33)*g23[ijk] + delG23*g33[ijk])
;

R33
=
gammado331*Gfromg1 + gammado332*Gfromg2 + gammado333*Gfromg3 + 
  (-0.5*deldelg1133 + gamma131*gammado113 + gamma231*gammado213 + 
     gamma331*gammado313 + 2.*(gamma113*gammado311 + gamma213*gammado321 + 
        gamma313*gammado331))*ginv11 + 
  (-deldelg1233 + gamma132*gammado113 + gamma131*gammado123 + 
     gamma232*gammado213 + gamma231*gammado223 + gamma332*gammado313 + 
     gamma331*gammado323 + 2.*(gamma123*gammado311 + gamma113*gammado312 + 
        gamma223*gammado321 + gamma213*gammado322 + gamma323*gammado331 + 
        gamma313*gammado332))*ginv12 + 
  (-deldelg1333 + gamma131*gammado133 + gamma231*gammado233 + 
     gamma133*(gammado113 + 2.*gammado311) + 
     (2.*gamma113 + gamma333)*gammado313 + 
     gamma233*(gammado213 + 2.*gammado321) + gamma331*gammado333 + 
     2.*(gamma213*gammado323 + gamma333*gammado331 + gamma313*gammado333))*
   ginv13 + (-0.5*deldelg2233 + gamma132*gammado123 + 
     gamma232*gammado223 + gamma332*gammado323 + 
     2.*(gamma123*gammado312 + gamma223*gammado322 + gamma323*gammado332))*
   ginv22 + (-deldelg2333 + gamma132*gammado133 + gamma232*gammado233 + 
     gamma133*(gammado123 + 2.*gammado312) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma333*(gammado323 + 2.*gammado332) + gamma332*gammado333 + 
     2.*(gamma123*gammado313 + gamma223*gammado323 + gamma323*gammado333))*
   ginv23 + (-0.5*deldelg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33 \
+ delG31*g13[ijk] + delG32*g23[ijk] + delG33*g33[ijk]
;

f
=
phi[ijk]
;

logpsi
=
0
;



/* conditional */
if (usepsi) {

logpsi
=
log(psi[ijk])
;

f
=
f + logpsi
;

df1
=
df1 + dpop1[ijk]
;

df2
=
df2 + dpop2[ijk]
;

df3
=
df3 + dpop3[ijk]
;

ddf11
=
ddf11 + ddpop11[ijk] - pow2(dpop1[ijk])
;

ddf12
=
ddf12 + ddpop12[ijk] - dpop1[ijk]*dpop2[ijk]
;

ddf13
=
ddf13 + ddpop13[ijk] - dpop1[ijk]*dpop3[ijk]
;

ddf22
=
ddf22 + ddpop22[ijk] - pow2(dpop2[ijk])
;

ddf23
=
ddf23 + ddpop23[ijk] - dpop2[ijk]*dpop3[ijk]
;

ddf33
=
ddf33 + ddpop33[ijk] - pow2(dpop3[ijk])
;

betadf
=
betadf + beta1[ijk]*dpop1[ijk] + beta2[ijk]*dpop2[ijk] + 
  beta3[ijk]*dpop3[ijk]
;

}
/* if (usepsi) */


cddf11
=
ddf11 - df1*gamma111 - df2*gamma211 - df3*gamma311
;

cddf12
=
ddf12 - df1*gamma112 - df2*gamma212 - df3*gamma312
;

cddf12
=
ddf12 - df1*gamma121 - df2*gamma221 - df3*gamma321
;

cddf13
=
ddf13 - df1*gamma113 - df2*gamma213 - df3*gamma313
;

cddf13
=
ddf13 - df1*gamma131 - df2*gamma231 - df3*gamma331
;

cddf22
=
ddf22 - df1*gamma122 - df2*gamma222 - df3*gamma322
;

cddf23
=
ddf23 - df1*gamma123 - df2*gamma223 - df3*gamma323
;

cddf23
=
ddf23 - df1*gamma132 - df2*gamma232 - df3*gamma332
;

cddf33
=
ddf33 - df1*gamma133 - df2*gamma233 - df3*gamma333
;

trcddf
=
cddf11*ginv11 + cddf22*ginv22 + 
  2.*(cddf12*ginv12 + cddf13*ginv13 + cddf23*ginv23) + cddf33*ginv33
;

psim4
=
exp(-4.*f)
;

Rphi11
=
-2.*(cddf11 + trcddf*g11[ijk]) + (4. - 4.*ginv11*g11[ijk])*pow2(df1) - 
  g11[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi12
=
df1*df2*(4. - 8.*ginv12*g12[ijk]) - 2.*(cddf12 + trcddf*g12[ijk]) - 
  g12[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi13
=
df1*(4.*df3 - 8.*df2*ginv12*g13[ijk]) - 2.*(cddf13 + trcddf*g13[ijk]) - 
  g13[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi22
=
-2.*(cddf22 + trcddf*g22[ijk]) + (4. - 4.*ginv22*g22[ijk])*pow2(df2) - 
  g22[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv33*pow2(df3)))
;

Rphi23
=
df2*(4.*df3 - 8.*df1*ginv12*g23[ijk]) - 2.*(cddf23 + trcddf*g23[ijk]) - 
  g23[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi33
=
-2.*(cddf33 + trcddf*g33[ijk]) - 
  g33[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2))) + 
  (4. - 4.*ginv33*g33[ijk])*pow2(df3)
;

cdda11
=
dda11 - da2*gamma211 - da3*gamma311 + 
  2.*((da2*df1 + da1*df2)*ginv12 + (da3*df1 + da1*df3)*ginv13 + 
     da2*df2*ginv22 + (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g11[ijk] \
+ da1*(-4.*df1 - gamma111 + 2.*df1*ginv11*g11[ijk])
;

cdda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma112 - da2*gamma212 - 
  da3*gamma312 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12[ijk]
;

cdda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma121 - da2*gamma221 - 
  da3*gamma321 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12[ijk]
;

cdda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma113 - da2*gamma213 - 
  da3*gamma313 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13[ijk]
;

cdda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma131 - da2*gamma231 - 
  da3*gamma331 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13[ijk]
;

cdda22
=
dda22 - da1*gamma122 - da3*gamma322 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + (da3*df2 + da2*df3)*ginv23 + 
     da3*df3*ginv33)*g22[ijk] + 
  da2*(-4.*df2 - gamma222 + 2.*df2*ginv22*g22[ijk])
;

cdda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma123 - da2*gamma223 - 
  da3*gamma323 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23[ijk]
;

cdda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma132 - da2*gamma232 - 
  da3*gamma332 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23[ijk]
;

cdda33
=
dda33 - da1*gamma133 - da2*gamma233 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23)*g33[ijk] + 
  da3*(-4.*df3 - gamma333 + 2.*df3*ginv33*g33[ijk])
;

trcdda
=
(cdda11*ginv11 + cdda22*ginv22 + 
    2.*(cdda12*ginv12 + cdda13*ginv13 + cdda23*ginv23) + cdda33*ginv33)*psim4
;

K[ijk]
=
forceKzerofactor*K[ijk]
;

dK1[ijk]
=
forceKzerofactor*dK1[ijk]
;

dK2[ijk]
=
forceKzerofactor*dK2[ijk]
;

dK3[ijk]
=
forceKzerofactor*dK3[ijk]
;

AA11
=
2.*(ginv23*A12[ijk]*A13[ijk] + 
     A11[ijk]*(ginv12*A12[ijk] + ginv13*A13[ijk])) + ginv11*pow2(A11[ijk]) + 
  ginv22*pow2(A12[ijk]) + ginv33*pow2(A13[ijk])
;

AA12
=
A12[ijk]*(ginv11*A11[ijk] + ginv22*A22[ijk]) + ginv33*A13[ijk]*A23[ijk] + 
  ginv13*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  ginv23*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  ginv12*(A11[ijk]*A22[ijk] + pow2(A12[ijk]))
;

AA13
=
ginv22*A12[ijk]*A23[ijk] + ginv12*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  A13[ijk]*(ginv11*A11[ijk] + ginv33*A33[ijk]) + 
  ginv23*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv13*(A11[ijk]*A33[ijk] + pow2(A13[ijk]))
;

AA22
=
2.*(ginv23*A22[ijk]*A23[ijk] + 
     A12[ijk]*(ginv12*A22[ijk] + ginv13*A23[ijk])) + ginv11*pow2(A12[ijk]) + 
  ginv22*pow2(A22[ijk]) + ginv33*pow2(A23[ijk])
;

AA23
=
ginv11*A12[ijk]*A13[ijk] + ginv12*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  A23[ijk]*(ginv22*A22[ijk] + ginv33*A33[ijk]) + 
  ginv13*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv23*(A22[ijk]*A33[ijk] + pow2(A23[ijk]))
;

AA33
=
2.*(ginv23*A23[ijk]*A33[ijk] + 
     A13[ijk]*(ginv12*A23[ijk] + ginv13*A33[ijk])) + ginv11*pow2(A13[ijk]) + 
  ginv22*pow2(A23[ijk]) + ginv33*pow2(A33[ijk])
;

AA
=
AA11*ginv11 + AA22*ginv22 + 2.*(AA12*ginv12 + AA13*ginv13 + AA23*ginv23) + 
  AA33*ginv33
;

Ainv11
=
2.*(ginv11*(ginv12*A12[ijk] + ginv13*A13[ijk]) + ginv12*ginv13*A23[ijk]) + 
  A11[ijk]*pow2(ginv11) + A22[ijk]*pow2(ginv12) + A33[ijk]*pow2(ginv13)
;

Ainv12
=
ginv11*(ginv12*A11[ijk] + ginv22*A12[ijk] + ginv23*A13[ijk]) + 
  ginv12*(ginv13*A13[ijk] + ginv22*A22[ijk] + ginv23*A23[ijk]) + 
  ginv13*(ginv22*A23[ijk] + ginv23*A33[ijk]) + A12[ijk]*pow2(ginv12)
;

Ainv13
=
ginv11*(ginv13*A11[ijk] + ginv23*A12[ijk] + ginv33*A13[ijk]) + 
  ginv12*(ginv13*A12[ijk] + ginv23*A22[ijk] + ginv33*A23[ijk]) + 
  ginv13*(ginv23*A23[ijk] + ginv33*A33[ijk]) + A13[ijk]*pow2(ginv13)
;

Ainv22
=
2.*(ginv12*(ginv22*A12[ijk] + ginv23*A13[ijk]) + ginv22*ginv23*A23[ijk]) + 
  A11[ijk]*pow2(ginv12) + A22[ijk]*pow2(ginv22) + A33[ijk]*pow2(ginv23)
;

Ainv23
=
ginv13*(ginv22*A12[ijk] + ginv23*A13[ijk]) + 
  ginv12*(ginv13*A11[ijk] + ginv23*A12[ijk] + ginv33*A13[ijk]) + 
  ginv22*(ginv23*A22[ijk] + ginv33*A23[ijk]) + ginv23*ginv33*A33[ijk] + 
  A23[ijk]*pow2(ginv23)
;

Ainv33
=
2.*(ginv13*(ginv23*A12[ijk] + ginv33*A13[ijk]) + ginv23*ginv33*A23[ijk]) + 
  A11[ijk]*pow2(ginv13) + A22[ijk]*pow2(ginv23) + A33[ijk]*pow2(ginv33)
;

divAinv1
=
-6.*(Ainv11*df1 + Ainv12*df2 + Ainv13*df3) - Ainv11*gamma111 - 
  Ainv12*(gamma112 + gamma121) - Ainv22*gamma122 - 
  Ainv13*(gamma113 + gamma131) - Ainv23*(gamma123 + gamma132) - 
  Ainv33*gamma133 + 0.66666666666666666666666666666666666667*
   (ginv11*dK1[ijk] + ginv12*dK2[ijk] + ginv13*dK3[ijk])
;

divAinv2
=
-6.*(Ainv12*df1 + Ainv22*df2 + Ainv23*df3) - Ainv11*gamma211 - 
  Ainv12*(gamma212 + gamma221) - Ainv22*gamma222 - 
  Ainv13*(gamma213 + gamma231) - Ainv23*(gamma223 + gamma232) - 
  Ainv33*gamma233 + 0.66666666666666666666666666666666666667*
   (ginv12*dK1[ijk] + ginv22*dK2[ijk] + ginv23*dK3[ijk])
;

divAinv3
=
-6.*(Ainv13*df1 + Ainv23*df2 + Ainv33*df3) - Ainv11*gamma311 - 
  Ainv12*(gamma312 + gamma321) - Ainv22*gamma322 - 
  Ainv13*(gamma313 + gamma331) - Ainv23*(gamma323 + gamma332) - 
  Ainv33*gamma333 + 0.66666666666666666666666666666666666667*
   (ginv13*dK1[ijk] + ginv23*dK2[ijk] + ginv33*dK3[ijk])
;

R
=
AA - 0.66666666666666666666666666666666666667*pow2(K[ijk])
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666666666666666666666667*divbeta
;

ootddivbeta1
=
0.33333333333333333333333333333333333333*(ddb111 + ddb122 + ddb133)
;

ootddivbeta2
=
0.33333333333333333333333333333333333333*(ddb121 + ddb222 + ddb233)
;

ootddivbeta3
=
0.33333333333333333333333333333333333333*(ddb131 + ddb232 + ddb333)
;

lieg11
=
betadg11 + (2.*db11 - totdivbeta)*g11[ijk] + 
  2.*(db12*g12[ijk] + db13*g13[ijk])
;

lieg12
=
betadg12 + db21*g11[ijk] + (db11 + db22 - totdivbeta)*g12[ijk] + 
  db23*g13[ijk] + db12*g22[ijk] + db13*g23[ijk]
;

lieg12
=
betadg21 + db21*g11[ijk] + (db11 + db22 - totdivbeta)*g12[ijk] + 
  db23*g13[ijk] + db12*g22[ijk] + db13*g23[ijk]
;

lieg13
=
betadg13 + db31*g11[ijk] + db32*g12[ijk] + 
  (db11 + db33 - totdivbeta)*g13[ijk] + db12*g23[ijk] + db13*g33[ijk]
;

lieg13
=
betadg31 + db31*g11[ijk] + db32*g12[ijk] + 
  (db11 + db33 - totdivbeta)*g13[ijk] + db12*g23[ijk] + db13*g33[ijk]
;

lieg22
=
betadg22 - totdivbeta*g22[ijk] + 
  2.*(db21*g12[ijk] + db22*g22[ijk] + db23*g23[ijk])
;

lieg23
=
betadg23 + db31*g12[ijk] + db21*g13[ijk] + db32*g22[ijk] + 
  (db22 + db33 - totdivbeta)*g23[ijk] + db23*g33[ijk]
;

lieg23
=
betadg32 + db31*g12[ijk] + db21*g13[ijk] + db32*g22[ijk] + 
  (db22 + db33 - totdivbeta)*g23[ijk] + db23*g33[ijk]
;

lieg33
=
betadg33 - totdivbeta*g33[ijk] + 
  2.*(db31*g13[ijk] + db32*g23[ijk] + db33*g33[ijk])
;

lieA11
=
betadA11 + (2.*db11 - totdivbeta)*A11[ijk] + 
  2.*(db12*A12[ijk] + db13*A13[ijk])
;

lieA12
=
betadA12 + db21*A11[ijk] + (db11 + db22 - totdivbeta)*A12[ijk] + 
  db23*A13[ijk] + db12*A22[ijk] + db13*A23[ijk]
;

lieA12
=
betadA21 + db21*A11[ijk] + (db11 + db22 - totdivbeta)*A12[ijk] + 
  db23*A13[ijk] + db12*A22[ijk] + db13*A23[ijk]
;

lieA13
=
betadA13 + db31*A11[ijk] + db32*A12[ijk] + 
  (db11 + db33 - totdivbeta)*A13[ijk] + db12*A23[ijk] + db13*A33[ijk]
;

lieA13
=
betadA31 + db31*A11[ijk] + db32*A12[ijk] + 
  (db11 + db33 - totdivbeta)*A13[ijk] + db12*A23[ijk] + db13*A33[ijk]
;

lieA22
=
betadA22 - totdivbeta*A22[ijk] + 
  2.*(db21*A12[ijk] + db22*A22[ijk] + db23*A23[ijk])
;

lieA23
=
betadA23 + db31*A12[ijk] + db21*A13[ijk] + db32*A22[ijk] + 
  (db22 + db33 - totdivbeta)*A23[ijk] + db23*A33[ijk]
;

lieA23
=
betadA32 + db31*A12[ijk] + db21*A13[ijk] + db32*A22[ijk] + 
  (db22 + db33 - totdivbeta)*A23[ijk] + db23*A33[ijk]
;

lieA33
=
betadA33 - totdivbeta*A33[ijk] + 
  2.*(db31*A13[ijk] + db32*A23[ijk] + db33*A33[ijk])
;

lieK
=
betadK
;

liephi
=
betadf + 0.16666666666666666666666666666666666667*divbeta
;

pseudolieG1
=
betadG1 - db11*Gfromg1 - db21*Gfromg2 - db31*Gfromg3 + ddb221*ginv22 + 
  2.*ddb231*ginv23 + ddb331*ginv33 + ginv11*(ddb111 + ootddivbeta1) + 
  ginv12*(2.*ddb121 + ootddivbeta2) + ginv13*(2.*ddb131 + ootddivbeta3) + 
  Gfromg1*totdivbeta
;

pseudolieG2
=
betadG2 - db12*Gfromg1 - db22*Gfromg2 - db32*Gfromg3 + ddb112*ginv11 + 
  2.*ddb132*ginv13 + ddb332*ginv33 + ginv12*(2.*ddb122 + ootddivbeta1) + 
  ginv22*(ddb222 + ootddivbeta2) + ginv23*(2.*ddb232 + ootddivbeta3) + 
  Gfromg2*totdivbeta
;

pseudolieG3
=
betadG3 - db13*Gfromg1 - db23*Gfromg2 - db33*Gfromg3 + ddb113*ginv11 + 
  2.*ddb123*ginv12 + ddb223*ginv22 + ginv13*(2.*ddb133 + ootddivbeta1) + 
  ginv23*(2.*ddb233 + ootddivbeta2) + ginv33*(ddb333 + ootddivbeta3) + 
  Gfromg3*totdivbeta
;

rg11
=
lieg11 - 2.*A11[ijk]*alpha[ijk]
;

rg12
=
lieg12 - 2.*A12[ijk]*alpha[ijk]
;

rg13
=
lieg13 - 2.*A13[ijk]*alpha[ijk]
;

rg22
=
lieg22 - 2.*A22[ijk]*alpha[ijk]
;

rg23
=
lieg23 - 2.*A23[ijk]*alpha[ijk]
;

rg33
=
lieg33 - 2.*A33[ijk]*alpha[ijk]
;

rA11
=
lieA11 + psim4*(-cdda11 + R11*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g11[ijk] + 
  alpha[ijk]*(-2.*AA11 + psim4*Rphi11 - 
     0.33333333333333333333333333333333333333*R*g11[ijk] + A11[ijk]*K[ijk])
;

rA12
=
lieA12 + psim4*(-cdda12 + R12*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g12[ijk] + 
  alpha[ijk]*(-2.*AA12 + psim4*Rphi12 - 
     0.33333333333333333333333333333333333333*R*g12[ijk] + A12[ijk]*K[ijk])
;

rA13
=
lieA13 + psim4*(-cdda13 + R13*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g13[ijk] + 
  alpha[ijk]*(-2.*AA13 + psim4*Rphi13 - 
     0.33333333333333333333333333333333333333*R*g13[ijk] + A13[ijk]*K[ijk])
;

rA22
=
lieA22 + psim4*(-cdda22 + R22*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g22[ijk] + 
  alpha[ijk]*(-2.*AA22 + psim4*Rphi22 - 
     0.33333333333333333333333333333333333333*R*g22[ijk] + A22[ijk]*K[ijk])
;

rA23
=
lieA23 + psim4*(-cdda23 + R23*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g23[ijk] + 
  alpha[ijk]*(-2.*AA23 + psim4*Rphi23 - 
     0.33333333333333333333333333333333333333*R*g23[ijk] + A23[ijk]*K[ijk])
;

rA33
=
lieA33 + psim4*(-cdda33 + R33*alpha[ijk]) + 
  0.33333333333333333333333333333333333333*trcdda*g33[ijk] + 
  alpha[ijk]*(-2.*AA33 + psim4*Rphi33 - 
     0.33333333333333333333333333333333333333*R*g33[ijk] + A33[ijk]*K[ijk])
;

rG1
=
pseudolieG1 - 2.*(Ainv11*da1 + Ainv12*da2 + Ainv13*da3 + divAinv1*alpha[ijk])
;

rG2
=
pseudolieG2 - 2.*(Ainv12*da1 + Ainv22*da2 + Ainv23*da3 + divAinv2*alpha[ijk])
;

rG3
=
pseudolieG3 - 2.*(Ainv13*da1 + Ainv23*da2 + Ainv33*da3 + divAinv3*alpha[ijk])
;

rK
=
lieK - trcdda + alpha[ijk]*(AA + 
     0.33333333333333333333333333333333333333*pow2(K[ijk]))
;

rphi
=
liephi - 0.16666666666666666666666666666666666667*alpha[ijk]*K[ijk]
;

ralpha0
=
(6.*liephi*oplogwithshift + alpha[ijk]*(-K[ijk] + subtractK0*K0[ijk]))*
  Power(psi[ijk],lapsepsipower)
;

ralpha
=
nonconstantlapse*ralpha0*(lapseharmonicf*oploglapse + 
    (8.*oploglapse2)/(9. - 3.*alpha[ijk]) + harmoniclapse*alpha[ijk])
;



/* conditional */
if (densitizedLapse) {

E6alphaDensityWeightf
=
exp(6.*alphaDensityWeight*f)
;

ooE6alphaDensityWeightf
=
1./E6alphaDensityWeightf
;

ralphaDensity
=
ooE6alphaDensityWeightf*ralpha - 6.*alphaDensityWeight*rphi*alphaDensity[ijk]
;



/* conditional */
if (densitizedoplogWithoutShift) {

ralphaDensity
=
alphaDensity[ijk]*(alphaDensityWeight*alpha[ijk]*
     (K[ijk] - subtractK0*K0[ijk]) + 
    lapseharmonicf*(-K[ijk] + subtractK0*K0[ijk])*
     Power(psi[ijk],lapsepsipower))
;

}
/* if (densitizedoplogWithoutShift) */


ralphaDensity
=
nonconstantlapse*ralphaDensity
;

}
/* if (densitizedoplogWithoutShift) */


betaF
=
(shiftgammacoeff*Power(alpha[ijk],shiftalphapower))/
  Power(psi[ijk],shiftpsipower)
;

rbeta1
=
evolveshift*(betaF*gamma0factor + gamma2factor)*B1[ijk]
;

rbeta2
=
evolveshift*(betaF*gamma0factor + gamma2factor)*B2[ijk]
;

rbeta3
=
evolveshift*(betaF*gamma0factor + gamma2factor)*B3[ijk]
;

rB1
=
evolveshift*((gamma0factor + betaF*gamma2factor)*rG1 - shiftdriver*B1[ijk])
;

rB2
=
evolveshift*((gamma0factor + betaF*gamma2factor)*rG2 - shiftdriver*B2[ijk])
;

rB3
=
evolveshift*((gamma0factor + betaF*gamma2factor)*rG3 - shiftdriver*B3[ijk])
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

nA11[ijk]
=
dt*rA11 + pA11[ijk]
;

nA12[ijk]
=
dt*rA12 + pA12[ijk]
;

nA13[ijk]
=
dt*rA13 + pA13[ijk]
;

nA22[ijk]
=
dt*rA22 + pA22[ijk]
;

nA23[ijk]
=
dt*rA23 + pA23[ijk]
;

nA33[ijk]
=
dt*rA33 + pA33[ijk]
;

nG1[ijk]
=
dt*rG1 + pG1[ijk]
;

nG2[ijk]
=
dt*rG2 + pG2[ijk]
;

nG3[ijk]
=
dt*rG3 + pG3[ijk]
;

nK[ijk]
=
forceKzerofactor*(dt*rK + pK[ijk])
;

nphi[ijk]
=
dt*rphi + pphi[ijk]
;



/* conditional */
if (densitizedLapse) {

nalphaDensity[ijk]
=
dt*ralphaDensity + palphaDensity[ijk]
;

nf
=
logpsi + nphi[ijk]
;

nalpha[ijk]
=
exp(6.*alphaDensityWeight*nf)*nalphaDensity[ijk]
;


} else { /* if (!densitizedLapse) */

nalpha[ijk]
=
dt*ralpha + palpha[ijk]
;

}
/* if (densitizedLapse) */


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

detnginv
=
1/(2.*ng12[ijk]*ng13[ijk]*ng23[ijk] + ng11[ijk]*ng22[ijk]*ng33[ijk] - 
    ng33[ijk]*pow2(ng12[ijk]) - ng22[ijk]*pow2(ng13[ijk]) - 
    ng11[ijk]*pow2(ng23[ijk]))
;



/* conditional */
if (subtractA) {

traceA
=
detnginv*(-2.*nA23[ijk]*ng11[ijk]*ng23[ijk] + 
    nA11[ijk]*ng22[ijk]*ng33[ijk] + 
    ng11[ijk]*(nA33[ijk]*ng22[ijk] + nA22[ijk]*ng33[ijk]) + 
    2.*(ng13[ijk]*(nA23[ijk]*ng12[ijk] - nA13[ijk]*ng22[ijk] + 
          nA12[ijk]*ng23[ijk]) + 
       ng12[ijk]*(nA13[ijk]*ng23[ijk] - nA12[ijk]*ng33[ijk])) - 
    nA33[ijk]*pow2(ng12[ijk]) - nA22[ijk]*pow2(ng13[ijk]) - 
    nA11[ijk]*pow2(ng23[ijk]))
;

aux
=
-0.33333333333333333333333333333333333333*traceA
;

nA11[ijk]
=
nA11[ijk] + aux*ng11[ijk]
;

nA12[ijk]
=
nA12[ijk] + aux*ng12[ijk]
;

nA13[ijk]
=
nA13[ijk] + aux*ng13[ijk]
;

nA22[ijk]
=
nA22[ijk] + aux*ng22[ijk]
;

nA23[ijk]
=
nA23[ijk] + aux*ng23[ijk]
;

nA33[ijk]
=
nA33[ijk] + aux*ng33[ijk]
;

}
/* if (subtractA) */




/* conditional */
if (normalizedetg) {

aux
=
Power(detnginv,0.3333333333333333)
;

ng11[ijk]
=
aux*ng11[ijk]
;

ng12[ijk]
=
aux*ng12[ijk]
;

ng13[ijk]
=
aux*ng13[ijk]
;

ng22[ijk]
=
aux*ng22[ijk]
;

ng23[ijk]
=
aux*ng23[ijk]
;

ng33[ijk]
=
aux*ng33[ijk]
;

}
/* if (normalizedetg) */



} else { /* if (!normalizedetg) */

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

nA11[ijk]
=
rA11
;

nA12[ijk]
=
rA12
;

nA13[ijk]
=
rA13
;

nA22[ijk]
=
rA22
;

nA23[ijk]
=
rA23
;

nA33[ijk]
=
rA33
;

nG1[ijk]
=
rG1
;

nG2[ijk]
=
rG2
;

nG3[ijk]
=
rG3
;

nK[ijk]
=
forceKzerofactor*rK
;

nphi[ijk]
=
rphi
;



/* conditional */
if (densitizedLapse) {

nalphaDensity[ijk]
=
ralphaDensity
;

nalpha[ijk]
=
E6alphaDensityWeightf*(ralphaDensity + 
    6.*alphaDensityWeight*rphi*alphaDensity[ijk])
;


errorexit("error in nalphaDensity: we should use the new E6alphaDensityWeightf, but we are using the old one instead."); 


} else { /* if (!densitizedLapse) */

nalpha[ijk]
=
ralpha
;

}
/* if (densitizedLapse) */


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
/* if (densitizedLapse) */


} /* end of points */
} /* end of boxes */


}  /* end of function */

/* BSSN_rhs.c */
/* nvars = 215, n* = 2270,  n/ = 58,  n+ = 2066, n = 4394, O = 1 */
