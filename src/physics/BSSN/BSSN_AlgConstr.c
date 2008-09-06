/* BSSN_AlgConstr.c */
/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, 5.9.2008 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "BSSN.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void BSSN_AlgConstr(tVarList *ucur)
{
tGrid *grid = ucur->grid;
int bi;

const int subtractA     = Getv("BSSN_subtractA", "yes");
const int normalizedetg = Getv("BSSN_normalizedetg", "yes");
const double oothree = 1.0/3.0;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;



double *g11 = vlldataptr(ucur, box, 0);
double *g12 = vlldataptr(ucur, box, 1);
double *g13 = vlldataptr(ucur, box, 2);
double *g22 = vlldataptr(ucur, box, 3);
double *g23 = vlldataptr(ucur, box, 4);
double *g33 = vlldataptr(ucur, box, 5);
double *A11 = vlldataptr(ucur, box, 6);
double *A12 = vlldataptr(ucur, box, 7);
double *A13 = vlldataptr(ucur, box, 8);
double *A22 = vlldataptr(ucur, box, 9);
double *A23 = vlldataptr(ucur, box, 10);
double *A33 = vlldataptr(ucur, box, 11);


double aux;
double detg;
double detginv;
double trA;



/* Jetzt geht's los! */

forallpoints(box, ijk) { 

detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;



/* conditional */
if (normalizedetg) {

aux
=
Power(detg,-oothree)
;

g11[ijk]
=
aux*g11[ijk]
;

g12[ijk]
=
aux*g12[ijk]
;

g13[ijk]
=
aux*g13[ijk]
;

g22[ijk]
=
aux*g22[ijk]
;

g23[ijk]
=
aux*g23[ijk]
;

g33[ijk]
=
aux*g33[ijk]
;

detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;

}
/* if (normalizedetg) */


detginv
=
1/detg
;

trA
=
detginv*(-2.*A23[ijk]*g11[ijk]*g23[ijk] + A11[ijk]*g22[ijk]*g33[ijk] + 
    g11[ijk]*(A33[ijk]*g22[ijk] + A22[ijk]*g33[ijk]) + 
    2.*(g13[ijk]*(A23[ijk]*g12[ijk] - A13[ijk]*g22[ijk] + 
          A12[ijk]*g23[ijk]) + 
       g12[ijk]*(A13[ijk]*g23[ijk] - A12[ijk]*g33[ijk])) - 
    A33[ijk]*pow2(g12[ijk]) - A22[ijk]*pow2(g13[ijk]) - 
    A11[ijk]*pow2(g23[ijk]))
;



/* conditional */
if (subtractA) {

aux
=
-(oothree*trA)
;

A11[ijk]
=
A11[ijk] + aux*g11[ijk]
;

A12[ijk]
=
A12[ijk] + aux*g12[ijk]
;

A13[ijk]
=
A13[ijk] + aux*g13[ijk]
;

A22[ijk]
=
A22[ijk] + aux*g22[ijk]
;

A23[ijk]
=
A23[ijk] + aux*g23[ijk]
;

A33[ijk]
=
A33[ijk] + aux*g33[ijk]
;

}
/* if (subtractA) */



/* no storage at the moment */ 

} /* end of points */
} /* end of boxes */


}  /* end of function */

/* BSSN_AlgConstr.c */
/* nvars = 14, n* = 90,  n/ = 27,  n+ = 30, n = 147, O = 1 */
