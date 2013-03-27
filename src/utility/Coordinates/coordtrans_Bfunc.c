/* coordtrans_Bfunc.c */
/* Copyright (C) 2013 Wolfgang Tichy, 27.3.2013 */
/* Produced with Mathematica */

#include "sgrid.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     (log((double) (x)))
#define Sqrt(x)    (sqrt((double) (x)))
#define Sin(x)     (sin((double) (x)))
#define Cos(x)     (cos((double) (x)))
#define Csc(x)     (1.0/sin((double) (x)))
#define Cot(x)     (1.0/tan((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))


extern int Coordinates_AnsorgNS_b_ParIndex;


/* func with trafo and 1st and 2nd derivs */
void ABphi_dABphi_ddABphi_of_Bfunc_AfBfPhi(tBox *box, int ind, double Af, double Bf, double Phi, double ABphi[4],double dABphi[4][4],double ddABphi[4][4][4])
{
double ss = 0.5;
double aa = 2.0*ss-2.0;
double bb = 3.0-3.0*ss;
double es = 0.05;

/* coord transforms */
ABphi[1] = Af;
ABphi[2] = aa*Power(Bf,3) + Bf*ss + bb*pow2(Bf);
ABphi[3] = Phi;

/* 1st derivs */
dABphi[1][1] = 1.;
dABphi[1][2] = 0;
dABphi[1][3] = 0;
dABphi[2][1] = 0;
dABphi[2][2] = 2.*bb*Bf + ss + 3.*aa*pow2(Bf);
dABphi[2][3] = 0;
dABphi[3][1] = 0;
dABphi[3][2] = 0;
dABphi[3][3] = 1.;

/* 2nd derivs */
ddABphi[1][1][1] = 0;
ddABphi[1][1][2] = 0;
ddABphi[1][1][3] = 0;
ddABphi[1][2][2] = 0;
ddABphi[1][2][3] = 0;
ddABphi[1][3][3] = 0;
ddABphi[2][1][1] = 0;
ddABphi[2][1][2] = 0;
ddABphi[2][1][3] = 0;
ddABphi[2][2][2] = 2.*(bb + 3.*aa*Bf);
ddABphi[2][2][3] = 0;
ddABphi[2][3][3] = 0;
ddABphi[3][1][1] = 0;
ddABphi[3][1][2] = 0;
ddABphi[3][1][3] = 0;
ddABphi[3][2][2] = 0;
ddABphi[3][2][3] = 0;
ddABphi[3][3][3] = 0;
}

/* end of: coordtrans_Bfunc.c */
/* nvars = 0, n* = 29,  n/ = 17,  n+ = 7, n = 53, O = 0 */
