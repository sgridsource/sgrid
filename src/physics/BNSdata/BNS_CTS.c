/* BNS_CTS.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 17.12.2007 */
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
double *FB1 = vlldataptr(vlFu, box, 1);
double *FB2 = vlldataptr(vlFu, box, 2);
double *FB3 = vlldataptr(vlFu, box, 3);
double *FalphaP = vlldataptr(vlFu, box, 4);
double *FSigma = vlldataptr(vlFu, box, 5);
double *Psi = vlldataptr(vlu, box, 0);
double *B1 = vlldataptr(vlu, box, 1);
double *B2 = vlldataptr(vlu, box, 2);
double *B3 = vlldataptr(vlu, box, 3);
double *alphaP = vlldataptr(vlu, box, 4);
double *Sigma = vlldataptr(vlu, box, 5);
double *FLPsi = vlldataptr(vlJdu, box, 0);
double *FLB1 = vlldataptr(vlJdu, box, 1);
double *FLB2 = vlldataptr(vlJdu, box, 2);
double *FLB3 = vlldataptr(vlJdu, box, 3);
double *FLalphaP = vlldataptr(vlJdu, box, 4);
double *FLSigma = vlldataptr(vlJdu, box, 5);
double *LPsi = vlldataptr(vldu, box, 0);
double *LB1 = vlldataptr(vldu, box, 1);
double *LB2 = vlldataptr(vldu, box, 2);
double *LB3 = vlldataptr(vldu, box, 3);
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
double *dLB11 = vlldataptr(vlduDerivs, box, 9);
double *dLB12 = vlldataptr(vlduDerivs, box, 10);
double *dLB13 = vlldataptr(vlduDerivs, box, 11);
double *dLB21 = vlldataptr(vlduDerivs, box, 12);
double *dLB22 = vlldataptr(vlduDerivs, box, 13);
double *dLB23 = vlldataptr(vlduDerivs, box, 14);
double *dLB31 = vlldataptr(vlduDerivs, box, 15);
double *dLB32 = vlldataptr(vlduDerivs, box, 16);
double *dLB33 = vlldataptr(vlduDerivs, box, 17);
double *ddLB111 = vlldataptr(vlduDerivs, box, 18);
double *ddLB112 = vlldataptr(vlduDerivs, box, 19);
double *ddLB113 = vlldataptr(vlduDerivs, box, 20);
double *ddLB122 = vlldataptr(vlduDerivs, box, 21);
double *ddLB123 = vlldataptr(vlduDerivs, box, 22);
double *ddLB133 = vlldataptr(vlduDerivs, box, 23);
double *ddLB211 = vlldataptr(vlduDerivs, box, 24);
double *ddLB212 = vlldataptr(vlduDerivs, box, 25);
double *ddLB213 = vlldataptr(vlduDerivs, box, 26);
double *ddLB222 = vlldataptr(vlduDerivs, box, 27);
double *ddLB223 = vlldataptr(vlduDerivs, box, 28);
double *ddLB233 = vlldataptr(vlduDerivs, box, 29);
double *ddLB311 = vlldataptr(vlduDerivs, box, 30);
double *ddLB312 = vlldataptr(vlduDerivs, box, 31);
double *ddLB313 = vlldataptr(vlduDerivs, box, 32);
double *ddLB322 = vlldataptr(vlduDerivs, box, 33);
double *ddLB323 = vlldataptr(vlduDerivs, box, 34);
double *ddLB333 = vlldataptr(vlduDerivs, box, 35);
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
int index_B1 = (vlu)->index[1];
int index_B2 = (vlu)->index[2];
int index_B3 = (vlu)->index[3];
int index_alphaP = (vlu)->index[4];
int index_Sigma = (vlu)->index[5];
int index_LPsi = (vldu)->index[0];
int index_LB1 = (vldu)->index[1];
int index_LB2 = (vldu)->index[2];
int index_LB3 = (vldu)->index[3];
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
int index_dLB11 = (vlduDerivs)->index[9];
int index_dLB12 = (vlduDerivs)->index[10];
int index_dLB13 = (vlduDerivs)->index[11];
int index_dLB21 = (vlduDerivs)->index[12];
int index_dLB22 = (vlduDerivs)->index[13];
int index_dLB23 = (vlduDerivs)->index[14];
int index_dLB31 = (vlduDerivs)->index[15];
int index_dLB32 = (vlduDerivs)->index[16];
int index_dLB33 = (vlduDerivs)->index[17];
int index_ddLB111 = (vlduDerivs)->index[18];
int index_ddLB112 = (vlduDerivs)->index[19];
int index_ddLB113 = (vlduDerivs)->index[20];
int index_ddLB122 = (vlduDerivs)->index[21];
int index_ddLB123 = (vlduDerivs)->index[22];
int index_ddLB133 = (vlduDerivs)->index[23];
int index_ddLB211 = (vlduDerivs)->index[24];
int index_ddLB212 = (vlduDerivs)->index[25];
int index_ddLB213 = (vlduDerivs)->index[26];
int index_ddLB222 = (vlduDerivs)->index[27];
int index_ddLB223 = (vlduDerivs)->index[28];
int index_ddLB233 = (vlduDerivs)->index[29];
int index_ddLB311 = (vlduDerivs)->index[30];
int index_ddLB312 = (vlduDerivs)->index[31];
int index_ddLB313 = (vlduDerivs)->index[32];
int index_ddLB322 = (vlduDerivs)->index[33];
int index_ddLB323 = (vlduDerivs)->index[34];
int index_ddLB333 = (vlduDerivs)->index[35];
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
int index_BNSdata_Bxx = Ind("BNSdata_Bxx");
double *dB11 = box->v[index_BNSdata_Bxx + 0];
double *dB12 = box->v[index_BNSdata_Bxx + 1];
double *dB13 = box->v[index_BNSdata_Bxx + 2];
double *dB21 = box->v[index_BNSdata_Bxx + 3];
double *dB22 = box->v[index_BNSdata_Bxx + 4];
double *dB23 = box->v[index_BNSdata_Bxx + 5];
double *dB31 = box->v[index_BNSdata_Bxx + 6];
double *dB32 = box->v[index_BNSdata_Bxx + 7];
double *dB33 = box->v[index_BNSdata_Bxx + 8];
int index_BNSdata_Bxxx = Ind("BNSdata_Bxxx");
double *ddB111 = box->v[index_BNSdata_Bxxx + 0];
double *ddB112 = box->v[index_BNSdata_Bxxx + 1];
double *ddB113 = box->v[index_BNSdata_Bxxx + 2];
double *ddB122 = box->v[index_BNSdata_Bxxx + 3];
double *ddB123 = box->v[index_BNSdata_Bxxx + 4];
double *ddB133 = box->v[index_BNSdata_Bxxx + 5];
double *ddB211 = box->v[index_BNSdata_Bxxx + 6];
double *ddB212 = box->v[index_BNSdata_Bxxx + 7];
double *ddB213 = box->v[index_BNSdata_Bxxx + 8];
double *ddB222 = box->v[index_BNSdata_Bxxx + 9];
double *ddB223 = box->v[index_BNSdata_Bxxx + 10];
double *ddB233 = box->v[index_BNSdata_Bxxx + 11];
double *ddB311 = box->v[index_BNSdata_Bxxx + 12];
double *ddB312 = box->v[index_BNSdata_Bxxx + 13];
double *ddB313 = box->v[index_BNSdata_Bxxx + 14];
double *ddB322 = box->v[index_BNSdata_Bxxx + 15];
double *ddB323 = box->v[index_BNSdata_Bxxx + 16];
double *ddB333 = box->v[index_BNSdata_Bxxx + 17];
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


FirstAndSecondDerivsOf_Sa(box, index_B1,    Ind("BNSdata_Bxx"), Ind("BNSdata_Bxxx")); 


FirstAndSecondDerivsOf_S(box, index_alphaP,    Ind("BNSdata_alphaPx"), Ind("BNSdata_alphaPxx")); 


FirstAndSecondDerivsOf_S(box, index_Sigma,    Ind("BNSdata_Sigmax"), Ind("BNSdata_Sigmaxx")); 


} else { /* if (!nonlin) */


FirstAndSecondDerivsOf_S(box, index_LPsi,       index_dLPsi1, index_ddLPsi11); 


FirstAndSecondDerivsOf_Sa(box, index_LB1,      index_dLB11, index_ddLB111); 


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

FB1[ijk]
=
ddB111[ijk] + ddB122[ijk] + ddB133[ijk]
;

FB2[ijk]
=
ddB211[ijk] + ddB222[ijk] + ddB233[ijk]
;

FB3[ijk]
=
ddB311[ijk] + ddB322[ijk] + ddB333[ijk]
;

FalphaP[ijk]
=
-rh2 + ddalphaP11[ijk] + ddalphaP22[ijk] + ddalphaP33[ijk]
;

FSigma[ijk]
=
ddSigma11[ijk] + ddSigma22[ijk] + ddSigma33[ijk]
;


} else { /* if (!nonlin) */

FLPsi[ijk]
=
ddLPsi11[ijk] + ddLPsi22[ijk] + ddLPsi33[ijk]
;

FLB1[ijk]
=
ddLB111[ijk] + ddLB122[ijk] + ddLB133[ijk]
;

FLB2[ijk]
=
ddLB211[ijk] + ddLB222[ijk] + ddLB233[ijk]
;

FLB3[ijk]
=
ddLB311[ijk] + ddLB322[ijk] + ddLB333[ijk]
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
