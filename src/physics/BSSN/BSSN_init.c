/* BSSN_init.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 9.5.2005 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BSSN.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BSSN_init(tGrid *grid, int igb, int iK, int ipsi, int igt, int iAt, int iG, int itrK, int iphi, int i_alpha, int i_alphaDensity)
{
int bi;

double alphaDensityWeight = Getd("BSSN_alphaDensityWeight");
int usepsi = 1;
for(bi = 0; bi < grid->nboxes; bi++)
{
tBox *box = grid->box[bi];
int ijk;


double *gb11 = box->v[igb+0];
double *gb12 = box->v[igb+1];
double *gb13 = box->v[igb+2];
double *gb22 = box->v[igb+3];
double *gb23 = box->v[igb+4];
double *gb33 = box->v[igb+5];
double *K11 = box->v[iK+0];
double *K12 = box->v[iK+1];
double *K13 = box->v[iK+2];
double *K22 = box->v[iK+3];
double *K23 = box->v[iK+4];
double *K33 = box->v[iK+5];
double *psi = box->v[ipsi+0];
double *gt11 = box->v[igt+0];
double *gt12 = box->v[igt+1];
double *gt13 = box->v[igt+2];
double *gt22 = box->v[igt+3];
double *gt23 = box->v[igt+4];
double *gt33 = box->v[igt+5];
double *At11 = box->v[iAt+0];
double *At12 = box->v[iAt+1];
double *At13 = box->v[iAt+2];
double *At22 = box->v[iAt+3];
double *At23 = box->v[iAt+4];
double *At33 = box->v[iAt+5];
double *G1 = box->v[iG+0];
double *G2 = box->v[iG+1];
double *G3 = box->v[iG+2];
double *K = box->v[itrK+0];
double *phi = box->v[iphi+0];
double *alpha = box->v[i_alpha+0];
double *alphaDensity = box->v[i_alphaDensity+0];
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


double detgb;
double detgtinv;
double p;
double p0;
double Gtdown1;
double Gtdown2;
double Gtdown3;
double gtinv11;
double gtinv12;
double gtinv13;
double gtinv22;
double gtinv23;
double gtinv33;
double Kt11;
double Kt12;
double Kt13;
double Kt22;
double Kt23;
double Kt33;



/* Jetzt geht's los! */

forallpoints(box, ijk) { 

detgb
=
2.*gb12[ijk]*gb13[ijk]*gb23[ijk] + gb11[ijk]*gb22[ijk]*gb33[ijk] - 
  gb33[ijk]*pow2(gb12[ijk]) - gb22[ijk]*pow2(gb13[ijk]) - 
  gb11[ijk]*pow2(gb23[ijk])
;

p
=
Power(detgb,-0.3333333333333333)
;

p0
=
Cal(usepsi,Power(psi[ijk],-4),1.)
;

gt11[ijk]
=
p*gb11[ijk]
;

gt12[ijk]
=
p*gb12[ijk]
;

gt13[ijk]
=
p*gb13[ijk]
;

gt22[ijk]
=
p*gb22[ijk]
;

gt23[ijk]
=
p*gb23[ijk]
;

gt33[ijk]
=
p*gb33[ijk]
;

Kt11
=
p*p0*K11[ijk]
;

Kt12
=
p*p0*K12[ijk]
;

Kt13
=
p*p0*K13[ijk]
;

Kt22
=
p*p0*K22[ijk]
;

Kt23
=
p*p0*K23[ijk]
;

Kt33
=
p*p0*K33[ijk]
;

alphaDensity[ijk]
=
Power(p*p0,1.5*alphaDensityWeight)*alpha[ijk]
;

detgtinv
=
1/(2.*gt12[ijk]*gt13[ijk]*gt23[ijk] + gt11[ijk]*gt22[ijk]*gt33[ijk] - 
    gt33[ijk]*pow2(gt12[ijk]) - gt22[ijk]*pow2(gt13[ijk]) - 
    gt11[ijk]*pow2(gt23[ijk]))
;

gtinv11
=
detgtinv*(gt22[ijk]*gt33[ijk] - pow2(gt23[ijk]))
;

gtinv12
=
detgtinv*(gt13[ijk]*gt23[ijk] - gt12[ijk]*gt33[ijk])
;

gtinv13
=
detgtinv*(-(gt13[ijk]*gt22[ijk]) + gt12[ijk]*gt23[ijk])
;

gtinv22
=
detgtinv*(gt11[ijk]*gt33[ijk] - pow2(gt13[ijk]))
;

gtinv23
=
detgtinv*(gt12[ijk]*gt13[ijk] - gt11[ijk]*gt23[ijk])
;

gtinv33
=
detgtinv*(gt11[ijk]*gt22[ijk] - pow2(gt12[ijk]))
;

phi[ijk]
=
-0.25*log(p)
;

K[ijk]
=
gtinv11*Kt11 + gtinv22*Kt22 + 2.*
   (gtinv12*Kt12 + gtinv13*Kt13 + gtinv23*Kt23) + gtinv33*Kt33
;

At11[ijk]
=
Kt11 - 0.33333333333333333333333333333333333333*gt11[ijk]*K[ijk]
;

At12[ijk]
=
Kt12 - 0.33333333333333333333333333333333333333*gt12[ijk]*K[ijk]
;

At13[ijk]
=
Kt13 - 0.33333333333333333333333333333333333333*gt13[ijk]*K[ijk]
;

At22[ijk]
=
Kt22 - 0.33333333333333333333333333333333333333*gt22[ijk]*K[ijk]
;

At23[ijk]
=
Kt23 - 0.33333333333333333333333333333333333333*gt23[ijk]*K[ijk]
;

At33[ijk]
=
Kt33 - 0.33333333333333333333333333333333333333*gt33[ijk]*K[ijk]
;


} /* end 1st forallpoints loop */ 


FirstDerivsOf_Sab(box, Ind("BSSN_g"),                     Ind("ADMvars_dgxxx")); 


forallpoints(box, ijk) { 

detgtinv
=
1/(2.*gt12[ijk]*gt13[ijk]*gt23[ijk] + gt11[ijk]*gt22[ijk]*gt33[ijk] - 
    gt33[ijk]*pow2(gt12[ijk]) - gt22[ijk]*pow2(gt13[ijk]) - 
    gt11[ijk]*pow2(gt23[ijk]))
;

gtinv11
=
detgtinv*(gt22[ijk]*gt33[ijk] - pow2(gt23[ijk]))
;

gtinv12
=
detgtinv*(gt13[ijk]*gt23[ijk] - gt12[ijk]*gt33[ijk])
;

gtinv13
=
detgtinv*(-(gt13[ijk]*gt22[ijk]) + gt12[ijk]*gt23[ijk])
;

gtinv22
=
detgtinv*(gt11[ijk]*gt33[ijk] - pow2(gt13[ijk]))
;

gtinv23
=
detgtinv*(gt12[ijk]*gt13[ijk] - gt11[ijk]*gt23[ijk])
;

gtinv33
=
detgtinv*(gt11[ijk]*gt22[ijk] - pow2(gt12[ijk]))
;

Gtdown1
=
gtinv11*dgt111[ijk] + gtinv12*(dgt112[ijk] + dgt121[ijk]) + 
  gtinv22*dgt122[ijk] + gtinv13*(dgt113[ijk] + dgt131[ijk]) + 
  gtinv23*(dgt123[ijk] + dgt132[ijk]) + gtinv33*dgt133[ijk]
;

Gtdown2
=
gtinv11*dgt121[ijk] + gtinv12*(dgt122[ijk] + dgt221[ijk]) + 
  gtinv22*dgt222[ijk] + gtinv13*(dgt123[ijk] + dgt231[ijk]) + 
  gtinv23*(dgt223[ijk] + dgt232[ijk]) + gtinv33*dgt233[ijk]
;

Gtdown3
=
gtinv11*dgt131[ijk] + gtinv12*(dgt132[ijk] + dgt231[ijk]) + 
  gtinv22*dgt232[ijk] + gtinv13*(dgt133[ijk] + dgt331[ijk]) + 
  gtinv23*(dgt233[ijk] + dgt332[ijk]) + gtinv33*dgt333[ijk]
;

G1[ijk]
=
Gtdown1*gtinv11 + Gtdown2*gtinv12 + Gtdown3*gtinv13
;

G2[ijk]
=
Gtdown1*gtinv12 + Gtdown2*gtinv22 + Gtdown3*gtinv23
;

G3[ijk]
=
Gtdown1*gtinv13 + Gtdown2*gtinv23 + Gtdown3*gtinv33
;

} /* end of points */
} /* end of boxes */


}  /* end of function */

/* BSSN_init.c */
/* nvars = 50, n* = 192,  n/ = 19,  n+ = 174, n = 385, O = 1 */
