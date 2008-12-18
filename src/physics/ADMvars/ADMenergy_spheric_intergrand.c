/* ADMenergy_spheric_intergrand.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 18.12.2008 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "ADMvars.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void ADMenergy_spheric_intergrand(tVarList *u)
{
tGrid *grid = u->grid;
int bi;

for(bi = 0; bi < grid->nboxes; bi++)
{
tBox *box = grid->box[bi];
int ijk;

FirstDerivsOf_Sab(box, Ind("gxx"), Ind("ADMvars_dgxxx"));

forallpoints(box, ijk)
{
double *Eadm = vlldataptr(u, box, 0);
double *gb11 = vlldataptr(u, box, 1);
double *gb12 = vlldataptr(u, box, 2);
double *gb13 = vlldataptr(u, box, 3);
double *gb22 = vlldataptr(u, box, 4);
double *gb23 = vlldataptr(u, box, 5);
double *gb33 = vlldataptr(u, box, 6);
double *x1 = vlldataptr(u, box, 7);
double *x2 = vlldataptr(u, box, 8);
double *x3 = vlldataptr(u, box, 9);
double *psi = vlldataptr(u, box, 10);
double *dpop1 = vlldataptr(u, box, 11);
double *dpop2 = vlldataptr(u, box, 12);
double *dpop3 = vlldataptr(u, box, 13);
int index_ADMvars_dgxxx = Ind("ADMvars_dgxxx");
double *dgb111 = box->v[index_ADMvars_dgxxx + 0];
double *dgb112 = box->v[index_ADMvars_dgxxx + 1];
double *dgb113 = box->v[index_ADMvars_dgxxx + 2];
double *dgb121 = box->v[index_ADMvars_dgxxx + 3];
double *dgb122 = box->v[index_ADMvars_dgxxx + 4];
double *dgb123 = box->v[index_ADMvars_dgxxx + 5];
double *dgb131 = box->v[index_ADMvars_dgxxx + 6];
double *dgb132 = box->v[index_ADMvars_dgxxx + 7];
double *dgb133 = box->v[index_ADMvars_dgxxx + 8];
double *dgb221 = box->v[index_ADMvars_dgxxx + 9];
double *dgb222 = box->v[index_ADMvars_dgxxx + 10];
double *dgb223 = box->v[index_ADMvars_dgxxx + 11];
double *dgb231 = box->v[index_ADMvars_dgxxx + 12];
double *dgb232 = box->v[index_ADMvars_dgxxx + 13];
double *dgb233 = box->v[index_ADMvars_dgxxx + 14];
double *dgb331 = box->v[index_ADMvars_dgxxx + 15];
double *dgb332 = box->v[index_ADMvars_dgxxx + 16];
double *dgb333 = box->v[index_ADMvars_dgxxx + 17];

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
double detg;
double detgb;
double f;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double integrand;
double n1;
double n2;
double n3;
double r2;



/* Jetzt geht's los! */
delg111
=
dgb111[ijk]
;

delg112
=
dgb121[ijk]
;

delg113
=
dgb131[ijk]
;

delg122
=
dgb221[ijk]
;

delg123
=
dgb231[ijk]
;

delg133
=
dgb331[ijk]
;

delg211
=
dgb112[ijk]
;

delg212
=
dgb122[ijk]
;

delg213
=
dgb132[ijk]
;

delg222
=
dgb222[ijk]
;

delg223
=
dgb232[ijk]
;

delg233
=
dgb332[ijk]
;

delg311
=
dgb113[ijk]
;

delg312
=
dgb123[ijk]
;

delg313
=
dgb133[ijk]
;

delg322
=
dgb223[ijk]
;

delg323
=
dgb233[ijk]
;

delg333
=
dgb333[ijk]
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

delg111
=
delg111*f + delf1*gb11[ijk]
;

delg112
=
delg112*f + delf1*gb12[ijk]
;

delg113
=
delg113*f + delf1*gb13[ijk]
;

delg122
=
delg122*f + delf1*gb22[ijk]
;

delg123
=
delg123*f + delf1*gb23[ijk]
;

delg133
=
delg133*f + delf1*gb33[ijk]
;

delg211
=
delg211*f + delf2*gb11[ijk]
;

delg212
=
delg212*f + delf2*gb12[ijk]
;

delg213
=
delg213*f + delf2*gb13[ijk]
;

delg222
=
delg222*f + delf2*gb22[ijk]
;

delg223
=
delg223*f + delf2*gb23[ijk]
;

delg233
=
delg233*f + delf2*gb33[ijk]
;

delg311
=
delg311*f + delf3*gb11[ijk]
;

delg312
=
delg312*f + delf3*gb12[ijk]
;

delg313
=
delg313*f + delf3*gb13[ijk]
;

delg322
=
delg322*f + delf3*gb22[ijk]
;

delg323
=
delg323*f + delf3*gb23[ijk]
;

delg333
=
delg333*f + delf3*gb33[ijk]
;

detgb
=
2.*gb12[ijk]*gb13[ijk]*gb23[ijk] + gb11[ijk]*gb22[ijk]*gb33[ijk] - 
  gb33[ijk]*pow2(gb12[ijk]) - gb22[ijk]*pow2(gb13[ijk]) - 
  gb11[ijk]*pow2(gb23[ijk])
;

detg
=
detgb*Power(f,3)
;

ginv11
=
(gb22[ijk]*gb33[ijk] - pow2(gb23[ijk]))/(detgb*f)
;

ginv12
=
(gb13[ijk]*gb23[ijk] - gb12[ijk]*gb33[ijk])/(detgb*f)
;

ginv13
=
(-(gb13[ijk]*gb22[ijk]) + gb12[ijk]*gb23[ijk])/(detgb*f)
;

ginv22
=
(gb11[ijk]*gb33[ijk] - pow2(gb13[ijk]))/(detgb*f)
;

ginv23
=
(gb12[ijk]*gb13[ijk] - gb11[ijk]*gb23[ijk])/(detgb*f)
;

ginv33
=
(gb11[ijk]*gb22[ijk] - pow2(gb12[ijk]))/(detgb*f)
;

r2
=
pow2(x1[ijk]) + pow2(x2[ijk]) + pow2(x3[ijk])
;



/* conditional */
if (r2) {

n1
=
x1[ijk]/sqrt(r2)
;

n2
=
x2[ijk]/sqrt(r2)
;

n3
=
x3[ijk]/sqrt(r2)
;

integrand
=
(detg*(ginv22*(0.*ginv23*((delg223 + delg322)*n2 + delg222*n3) + 
         0.0625*ginv33*((-delg233 + delg323)*n2 + (delg223 - delg322)*n3)) \
+ ginv12*(0.*ginv22*(delg222*n1 + (delg122 + delg212)*n2) + 
         ginv33*(0.0625*((-delg233 + delg323)*n1 + 
               (-delg133 + delg313)*n2) + 
            (0.0625*(delg123 + delg213) - 0.125*delg312)*n3) + 
         ginv13*((0.125*delg123 - 0.0625*(delg213 + delg312))*n1 + 
            0.0625*((-delg113 + delg311)*n2 + (-delg112 + delg211)*n3)) + 
         ginv23*((0.125*delg213 - 0.0625*(delg123 + delg312))*n2 + 
            0.0625*((-delg223 + delg322)*n1 + (delg122 - delg212)*n3))) + 
      ginv13*(ginv23*(0.0625*((delg233 - delg323)*n1 + 
               (delg133 - delg313)*n2) + 
            (-0.0625*(delg123 + delg213) + 0.125*delg312)*n3) + 
         0.*ginv33*(delg333*n1 + (delg133 + delg313)*n3) + 
         ginv22*((-0.125*delg213 + 0.0625*(delg123 + delg312))*n2 + 
            0.0625*((delg223 - delg322)*n1 + (-delg122 + delg212)*n3))) + 
      ginv11*(0.*(ginv12*((delg112 + delg211)*n1 + delg111*n2) + 
            ginv13*((delg113 + delg311)*n1 + delg111*n3)) + 
         ginv23*((-0.125*delg123 + 0.0625*(delg213 + delg312))*n1 + 
            0.0625*((delg113 - delg311)*n2 + (delg112 - delg211)*n3)) + 
         0.0625*(ginv22*((-delg122 + delg212)*n1 + 
               (delg112 - delg211)*n2) + 
            ginv33*((-delg133 + delg313)*n1 + (delg113 - delg311)*n3))) + 
      0.0625*(((delg122 - delg212)*n1 + (-delg112 + delg211)*n2)*
          pow2(ginv12) + ((delg133 - delg313)*n1 + 
            (-delg113 + delg311)*n3)*pow2(ginv13) + 
         ((delg233 - delg323)*n2 + (-delg223 + delg322)*n3)*pow2(ginv23)) + 
      0.*(ginv23*ginv33*(delg333*n2 + (delg233 + delg323)*n3) + 
         delg111*n1*pow2(ginv11) + delg222*n2*pow2(ginv22) + 
         delg333*n3*pow2(ginv33))))/PI
;


} else { /* if (!r2) */

integrand
=
0.
;

}
/* if (r2) */


Eadm[ijk]
=
integrand
;

} /* end of points */

} /* end of boxes */


}  /* end of function */

/* ADMenergy_spheric_intergrand.c */
/* nvars = 32, n* = 221,  n/ = 31,  n+ = 174, n = 426, O = 1 */
