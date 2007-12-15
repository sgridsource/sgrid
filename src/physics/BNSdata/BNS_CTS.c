/* BNS_CTS.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 15.12.2007 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BNS_CTS(tVarList *vlFu, tVarList *vlu,       tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs,      int nonlin)
{
tGrid *grid = vlu->grid;
int bi;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;



double *FPsi = vlldataptr(vlFu, box, 0);
double *Fbeta1 = vlldataptr(vlFu, box, 1);
double *Fbeta2 = vlldataptr(vlFu, box, 2);
double *Fbeta3 = vlldataptr(vlFu, box, 3);
double *FalphaP = vlldataptr(vlFu, box, 4);
double *FSigma = vlldataptr(vlFu, box, 5);
double *Psi = vlldataptr(vlu, box, 0);
double *beta1 = vlldataptr(vlu, box, 1);
double *beta2 = vlldataptr(vlu, box, 2);
double *beta3 = vlldataptr(vlu, box, 3);
double *alphaP = vlldataptr(vlu, box, 4);
double *Sigma = vlldataptr(vlu, box, 5);
double *FLPsi = vlldataptr(vlJdu, box, 0);
double *FLbeta1 = vlldataptr(vlJdu, box, 1);
double *FLbeta2 = vlldataptr(vlJdu, box, 2);
double *FLbeta3 = vlldataptr(vlJdu, box, 3);
double *FLalphaP = vlldataptr(vlJdu, box, 4);
double *FLSigma = vlldataptr(vlJdu, box, 5);
double *LPsi = vlldataptr(vldu, box, 0);
double *Lbeta1 = vlldataptr(vldu, box, 1);
double *Lbeta2 = vlldataptr(vldu, box, 2);
double *Lbeta3 = vlldataptr(vldu, box, 3);
double *LalphaP = vlldataptr(vldu, box, 4);
double *LSigma = vlldataptr(vldu, box, 5);
double *dLPsi1 = vlldataptr(vlduDerivs, box, 0);
double *dLPsi2 = vlldataptr(vlduDerivs, box, 1);
double *dLPsi3 = vlldataptr(vlduDerivs, box, 2);
double *ddLPsi11 = vlldataptr(vlduDerivs, box, 3);
double *ddLPsi12 = vlldataptr(vlduDerivs, box, 4);
double *ddLPsi13 = vlldataptr(vlduDerivs, box, 5);
double *ddLPsi22 = vlldataptr(vlduDerivs, box, 6);
double *ddLPsi23 = vlldataptr(vlduDerivs, box, 7);
double *ddLPsi33 = vlldataptr(vlduDerivs, box, 8);
double *dLbeta11 = vlldataptr(vlduDerivs, box, 9);
double *dLbeta12 = vlldataptr(vlduDerivs, box, 10);
double *dLbeta13 = vlldataptr(vlduDerivs, box, 11);
double *dLbeta21 = vlldataptr(vlduDerivs, box, 12);
double *dLbeta22 = vlldataptr(vlduDerivs, box, 13);
double *dLbeta23 = vlldataptr(vlduDerivs, box, 14);
double *dLbeta31 = vlldataptr(vlduDerivs, box, 15);
double *dLbeta32 = vlldataptr(vlduDerivs, box, 16);
double *dLbeta33 = vlldataptr(vlduDerivs, box, 17);
double *ddLbeta111 = vlldataptr(vlduDerivs, box, 18);
double *ddLbeta112 = vlldataptr(vlduDerivs, box, 19);
double *ddLbeta113 = vlldataptr(vlduDerivs, box, 20);
double *ddLbeta122 = vlldataptr(vlduDerivs, box, 21);
double *ddLbeta123 = vlldataptr(vlduDerivs, box, 22);
double *ddLbeta133 = vlldataptr(vlduDerivs, box, 23);
double *ddLbeta211 = vlldataptr(vlduDerivs, box, 24);
double *ddLbeta212 = vlldataptr(vlduDerivs, box, 25);
double *ddLbeta213 = vlldataptr(vlduDerivs, box, 26);
double *ddLbeta222 = vlldataptr(vlduDerivs, box, 27);
double *ddLbeta223 = vlldataptr(vlduDerivs, box, 28);
double *ddLbeta233 = vlldataptr(vlduDerivs, box, 29);
double *ddLbeta311 = vlldataptr(vlduDerivs, box, 30);
double *ddLbeta312 = vlldataptr(vlduDerivs, box, 31);
double *ddLbeta313 = vlldataptr(vlduDerivs, box, 32);
double *ddLbeta322 = vlldataptr(vlduDerivs, box, 33);
double *ddLbeta323 = vlldataptr(vlduDerivs, box, 34);
double *ddLbeta333 = vlldataptr(vlduDerivs, box, 35);
double *dLalphaP1 = vlldataptr(vlduDerivs, box, 36);
double *dLalphaP2 = vlldataptr(vlduDerivs, box, 37);
double *dLalphaP3 = vlldataptr(vlduDerivs, box, 38);
double *ddLalphaP11 = vlldataptr(vlduDerivs, box, 39);
double *ddLalphaP12 = vlldataptr(vlduDerivs, box, 40);
double *ddLalphaP13 = vlldataptr(vlduDerivs, box, 41);
double *ddLalphaP22 = vlldataptr(vlduDerivs, box, 42);
double *ddLalphaP23 = vlldataptr(vlduDerivs, box, 43);
double *ddLalphaP33 = vlldataptr(vlduDerivs, box, 44);
double *dLSigma1 = vlldataptr(vlduDerivs, box, 45);
double *dLSigma2 = vlldataptr(vlduDerivs, box, 46);
double *dLSigma3 = vlldataptr(vlduDerivs, box, 47);
double *ddLSigma11 = vlldataptr(vlduDerivs, box, 48);
double *ddLSigma12 = vlldataptr(vlduDerivs, box, 49);
double *ddLSigma13 = vlldataptr(vlduDerivs, box, 50);
double *ddLSigma22 = vlldataptr(vlduDerivs, box, 51);
double *ddLSigma23 = vlldataptr(vlduDerivs, box, 52);
double *ddLSigma33 = vlldataptr(vlduDerivs, box, 53);
int index_Psi = (vlu)->index[0];
int index_beta1 = (vlu)->index[1];
int index_beta2 = (vlu)->index[2];
int index_beta3 = (vlu)->index[3];
int index_alphaP = (vlu)->index[4];
int index_Sigma = (vlu)->index[5];
int index_LPsi = (vldu)->index[0];
int index_Lbeta1 = (vldu)->index[1];
int index_Lbeta2 = (vldu)->index[2];
int index_Lbeta3 = (vldu)->index[3];
int index_LalphaP = (vldu)->index[4];
int index_LSigma = (vldu)->index[5];
int index_dLPsi1 = (vlduDerivs)->index[0];
int index_dLPsi2 = (vlduDerivs)->index[1];
int index_dLPsi3 = (vlduDerivs)->index[2];
int index_ddLPsi11 = (vlduDerivs)->index[3];
int index_ddLPsi12 = (vlduDerivs)->index[4];
int index_ddLPsi13 = (vlduDerivs)->index[5];
int index_ddLPsi22 = (vlduDerivs)->index[6];
int index_ddLPsi23 = (vlduDerivs)->index[7];
int index_ddLPsi33 = (vlduDerivs)->index[8];
int index_dLbeta11 = (vlduDerivs)->index[9];
int index_dLbeta12 = (vlduDerivs)->index[10];
int index_dLbeta13 = (vlduDerivs)->index[11];
int index_dLbeta21 = (vlduDerivs)->index[12];
int index_dLbeta22 = (vlduDerivs)->index[13];
int index_dLbeta23 = (vlduDerivs)->index[14];
int index_dLbeta31 = (vlduDerivs)->index[15];
int index_dLbeta32 = (vlduDerivs)->index[16];
int index_dLbeta33 = (vlduDerivs)->index[17];
int index_ddLbeta111 = (vlduDerivs)->index[18];
int index_ddLbeta112 = (vlduDerivs)->index[19];
int index_ddLbeta113 = (vlduDerivs)->index[20];
int index_ddLbeta122 = (vlduDerivs)->index[21];
int index_ddLbeta123 = (vlduDerivs)->index[22];
int index_ddLbeta133 = (vlduDerivs)->index[23];
int index_ddLbeta211 = (vlduDerivs)->index[24];
int index_ddLbeta212 = (vlduDerivs)->index[25];
int index_ddLbeta213 = (vlduDerivs)->index[26];
int index_ddLbeta222 = (vlduDerivs)->index[27];
int index_ddLbeta223 = (vlduDerivs)->index[28];
int index_ddLbeta233 = (vlduDerivs)->index[29];
int index_ddLbeta311 = (vlduDerivs)->index[30];
int index_ddLbeta312 = (vlduDerivs)->index[31];
int index_ddLbeta313 = (vlduDerivs)->index[32];
int index_ddLbeta322 = (vlduDerivs)->index[33];
int index_ddLbeta323 = (vlduDerivs)->index[34];
int index_ddLbeta333 = (vlduDerivs)->index[35];
int index_dLalphaP1 = (vlduDerivs)->index[36];
int index_dLalphaP2 = (vlduDerivs)->index[37];
int index_dLalphaP3 = (vlduDerivs)->index[38];
int index_ddLalphaP11 = (vlduDerivs)->index[39];
int index_ddLalphaP12 = (vlduDerivs)->index[40];
int index_ddLalphaP13 = (vlduDerivs)->index[41];
int index_ddLalphaP22 = (vlduDerivs)->index[42];
int index_ddLalphaP23 = (vlduDerivs)->index[43];
int index_ddLalphaP33 = (vlduDerivs)->index[44];
int index_dLSigma1 = (vlduDerivs)->index[45];
int index_dLSigma2 = (vlduDerivs)->index[46];
int index_dLSigma3 = (vlduDerivs)->index[47];
int index_ddLSigma11 = (vlduDerivs)->index[48];
int index_ddLSigma12 = (vlduDerivs)->index[49];
int index_ddLSigma13 = (vlduDerivs)->index[50];
int index_ddLSigma22 = (vlduDerivs)->index[51];
int index_ddLSigma23 = (vlduDerivs)->index[52];
int index_ddLSigma33 = (vlduDerivs)->index[53];
int index_BNSdata_Psix = Ind("BNSdata_Psix");
double *dPsi1 = box->v[index_BNSdata_Psix + 0];
double *dPsi2 = box->v[index_BNSdata_Psix + 1];
double *dPsi3 = box->v[index_BNSdata_Psix + 2];
int index_BNSdata_Psixx = Ind("BNSdata_Psixx");
double *ddPsi11 = box->v[index_BNSdata_Psixx + 0];
double *ddPsi12 = box->v[index_BNSdata_Psixx + 1];
double *ddPsi13 = box->v[index_BNSdata_Psixx + 2];
double *ddPsi22 = box->v[index_BNSdata_Psixx + 3];
double *ddPsi23 = box->v[index_BNSdata_Psixx + 4];
double *ddPsi33 = box->v[index_BNSdata_Psixx + 5];
int index_BNSdata_betaxx = Ind("BNSdata_betaxx");
double *dbeta11 = box->v[index_BNSdata_betaxx + 0];
double *dbeta12 = box->v[index_BNSdata_betaxx + 1];
double *dbeta13 = box->v[index_BNSdata_betaxx + 2];
double *dbeta21 = box->v[index_BNSdata_betaxx + 3];
double *dbeta22 = box->v[index_BNSdata_betaxx + 4];
double *dbeta23 = box->v[index_BNSdata_betaxx + 5];
double *dbeta31 = box->v[index_BNSdata_betaxx + 6];
double *dbeta32 = box->v[index_BNSdata_betaxx + 7];
double *dbeta33 = box->v[index_BNSdata_betaxx + 8];
int index_BNSdata_betaxxx = Ind("BNSdata_betaxxx");
double *ddbeta111 = box->v[index_BNSdata_betaxxx + 0];
double *ddbeta112 = box->v[index_BNSdata_betaxxx + 1];
double *ddbeta113 = box->v[index_BNSdata_betaxxx + 2];
double *ddbeta122 = box->v[index_BNSdata_betaxxx + 3];
double *ddbeta123 = box->v[index_BNSdata_betaxxx + 4];
double *ddbeta133 = box->v[index_BNSdata_betaxxx + 5];
double *ddbeta211 = box->v[index_BNSdata_betaxxx + 6];
double *ddbeta212 = box->v[index_BNSdata_betaxxx + 7];
double *ddbeta213 = box->v[index_BNSdata_betaxxx + 8];
double *ddbeta222 = box->v[index_BNSdata_betaxxx + 9];
double *ddbeta223 = box->v[index_BNSdata_betaxxx + 10];
double *ddbeta233 = box->v[index_BNSdata_betaxxx + 11];
double *ddbeta311 = box->v[index_BNSdata_betaxxx + 12];
double *ddbeta312 = box->v[index_BNSdata_betaxxx + 13];
double *ddbeta313 = box->v[index_BNSdata_betaxxx + 14];
double *ddbeta322 = box->v[index_BNSdata_betaxxx + 15];
double *ddbeta323 = box->v[index_BNSdata_betaxxx + 16];
double *ddbeta333 = box->v[index_BNSdata_betaxxx + 17];
int index_BNSdata_alphaPx = Ind("BNSdata_alphaPx");
double *dalphaP1 = box->v[index_BNSdata_alphaPx + 0];
double *dalphaP2 = box->v[index_BNSdata_alphaPx + 1];
double *dalphaP3 = box->v[index_BNSdata_alphaPx + 2];
int index_BNSdata_alphaPxx = Ind("BNSdata_alphaPxx");
double *ddalphaP11 = box->v[index_BNSdata_alphaPxx + 0];
double *ddalphaP12 = box->v[index_BNSdata_alphaPxx + 1];
double *ddalphaP13 = box->v[index_BNSdata_alphaPxx + 2];
double *ddalphaP22 = box->v[index_BNSdata_alphaPxx + 3];
double *ddalphaP23 = box->v[index_BNSdata_alphaPxx + 4];
double *ddalphaP33 = box->v[index_BNSdata_alphaPxx + 5];
int index_BNSdata_Sigmax = Ind("BNSdata_Sigmax");
double *dSigma1 = box->v[index_BNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_BNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_BNSdata_Sigmax + 2];
int index_BNSdata_Sigmaxx = Ind("BNSdata_Sigmaxx");
double *ddSigma11 = box->v[index_BNSdata_Sigmaxx + 0];
double *ddSigma12 = box->v[index_BNSdata_Sigmaxx + 1];
double *ddSigma13 = box->v[index_BNSdata_Sigmaxx + 2];
double *ddSigma22 = box->v[index_BNSdata_Sigmaxx + 3];
double *ddSigma23 = box->v[index_BNSdata_Sigmaxx + 4];
double *ddSigma33 = box->v[index_BNSdata_Sigmaxx + 5];
int index_gxx = Ind("gxx");
double *g11 = box->v[index_gxx + 0];
double *g12 = box->v[index_gxx + 1];
double *g13 = box->v[index_gxx + 2];
double *g22 = box->v[index_gxx + 3];
double *g23 = box->v[index_gxx + 4];
double *g33 = box->v[index_gxx + 5];
int index_Kxx = Ind("Kxx");
double *K11 = box->v[index_Kxx + 0];
double *K12 = box->v[index_Kxx + 1];
double *K13 = box->v[index_Kxx + 2];
double *K22 = box->v[index_Kxx + 3];
double *K23 = box->v[index_Kxx + 4];
double *K33 = box->v[index_Kxx + 5];
int index_BNSdata_q = Ind("BNSdata_q");
double *q = box->v[index_BNSdata_q + 0];
int index_BNSdata_vRSx = Ind("BNSdata_vRSx");
double *vRS1 = box->v[index_BNSdata_vRSx + 0];
double *vRS2 = box->v[index_BNSdata_vRSx + 1];
double *vRS3 = box->v[index_BNSdata_vRSx + 2];





/* Jetzt geht's los! */


/* conditional */
if (nonlin) {


FirstAndSecondDerivsOf_S(box, index_Psi,     Ind("BNSdata_Psix"), Ind("BNSdata_Psixx")); 


FirstAndSecondDerivsOf_Sa(box, index_beta1,    Ind("BNSdata_betaxx"), Ind("BNSdata_betaxxx")); 


FirstAndSecondDerivsOf_S(box, index_alphaP,    Ind("BNSdata_alphaPx"), Ind("BNSdata_alphaPxx")); 


FirstAndSecondDerivsOf_S(box, index_Sigma,    Ind("BNSdata_Sigmax"), Ind("BNSdata_Sigmaxx")); 


} else { /* if (!nonlin) */


FirstAndSecondDerivsOf_S(box, index_LPsi,       index_dLPsi1, index_ddLPsi11); 


FirstAndSecondDerivsOf_Sa(box, index_Lbeta1,      index_dLbeta11, index_ddLbeta111); 


FirstAndSecondDerivsOf_S(box, index_LalphaP,      index_dLalphaP1, index_ddLalphaP11); 


FirstAndSecondDerivsOf_S(box, index_LSigma,      index_dLSigma1, index_ddLSigma11); 

}
/* if (nonlin) */



forallpoints(box, ijk) { 


         double xmax1 = grid->box[0]->x_of_X[1](                         (void *) grid->box[0], 0, 0.0,0.0,0.0);         double xmin1 = grid->box[0]->x_of_X[1](                         (void *) grid->box[0], 0, 0.0,1.0,0.0);         double xmax2 = grid->box[3]->x_of_X[1](                         (void *) grid->box[3], 0, 0.0,1.0,0.0);         double xmin2 = grid->box[3]->x_of_X[1](                         (void *) grid->box[3], 0, 0.0,0.0,0.0);         double R1  = 0.5*(xmax1-xmin1);         double R2  = 0.5*(xmax2-xmin2); double rh1 = 0.0;         double rh2 = 0.0;          if(bi==0 || bi==5)  rh1 = -3.0/(R1*R1*R1);         if(bi==3 || bi==4)  rh2 = -6.0/(R2*R2*R2);  



/* conditional */
if (nonlin) {

FPsi[ijk]
=
-rh1 + ddPsi11[ijk] + ddPsi22[ijk] + ddPsi33[ijk]
;

Fbeta1[ijk]
=
ddbeta111[ijk] + ddbeta122[ijk] + ddbeta133[ijk]
;

Fbeta2[ijk]
=
ddbeta211[ijk] + ddbeta222[ijk] + ddbeta233[ijk]
;

Fbeta3[ijk]
=
ddbeta311[ijk] + ddbeta322[ijk] + ddbeta333[ijk]
;

FalphaP[ijk]
=
ddalphaP11[ijk] + ddalphaP22[ijk] + ddalphaP33[ijk]
;

FSigma[ijk]
=
-rh2 + ddSigma11[ijk] + ddSigma22[ijk] + ddSigma33[ijk]
;


} else { /* if (!nonlin) */

FLPsi[ijk]
=
ddLPsi11[ijk] + ddLPsi22[ijk] + ddLPsi33[ijk]
;

FLbeta1[ijk]
=
ddLbeta111[ijk] + ddLbeta122[ijk] + ddLbeta133[ijk]
;

FLbeta2[ijk]
=
ddLbeta211[ijk] + ddLbeta222[ijk] + ddLbeta233[ijk]
;

FLbeta3[ijk]
=
ddLbeta311[ijk] + ddLbeta322[ijk] + ddLbeta333[ijk]
;

FLalphaP[ijk]
=
ddLalphaP11[ijk] + ddLalphaP22[ijk] + ddLalphaP33[ijk]
;

FLSigma[ijk]
=
ddLSigma11[ijk] + ddLSigma22[ijk] + ddLSigma33[ijk]
;

}
/* if (nonlin) */



} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* BNS_CTS.c */
/* nvars = 148, n* = 193,  n/ = 29,  n+ = 252, n = 474, O = 1 */
