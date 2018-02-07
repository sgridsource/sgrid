/* coordtrans_AnsorgXRphi.c */
/* Copyright (C) 2013 Wolfgang Tichy, 22.1.2018 */
/* Produced with Mathematica */

#include "sgrid.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     (log((double) (x)))
#define Sin(x)     (sin((double) (x)))
#define Cos(x)     (cos((double) (x)))
#define Csc(x)     (1.0/sin((double) (x)))
#define Cot(x)     (1.0/tan((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))


extern int Coordinates_AnsorgNS_b_ParIndex;


/* func with trafo and 1st and 2nd derivs */
void xyz_dxyz_ddxyz_of_AnsorgXRphi_XRphi(tBox *box, int ind, double X, double R, double phi, double xyz[4],double dxyz[4][4],double ddxyz[4][4][4])
{
double b = GetCachedNumValByParIndex(Coordinates_AnsorgNS_b_ParIndex);
double XsqrPlusRsqr = pow2(R) + pow2(X);
double DXsqrPlusRsqrDX = 2*X;
double DXsqrPlusRsqrDR = 2*R;
double DDXsqrPlusRsqrDXDX = 2;
double DDXsqrPlusRsqrDRDR = 2;
if(XsqrPlusRsqr==0.0) XsqrPlusRsqr=dequaleps*dequaleps*dequaleps;

/* coord transforms */
xyz[1] = 0.5*b*(-pow2(R) + pow2(X))*(1. + pow2inv(XsqrPlusRsqr));
xyz[2] = b*R*X*Cos(phi)*(-1. + pow2inv(XsqrPlusRsqr));
xyz[3] = b*R*X*(-1. + pow2inv(XsqrPlusRsqr))*Sin(phi);

/* 1st derivs */
dxyz[1][1] = (b*DXsqrPlusRsqrDX*(pow2(R) - 1.*pow2(X)))/Power(XsqrPlusRsqr,3) + b*X*(1. + pow2inv(XsqrPlusRsqr));
dxyz[1][2] = (b*DXsqrPlusRsqrDR*(pow2(R) - 1.*pow2(X)))/Power(XsqrPlusRsqr,3) - 1.*b*R*(1. + pow2inv(XsqrPlusRsqr));
dxyz[1][3] = 0;
dxyz[2][1] = b*R*Cos(phi)*(-1. - (2.*DXsqrPlusRsqrDX*X)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr));
dxyz[2][2] = b*X*Cos(phi)*(-1. - (2.*DXsqrPlusRsqrDR*R)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr));
dxyz[2][3] = -1.*b*R*X*(-1. + pow2inv(XsqrPlusRsqr))*Sin(phi);
dxyz[3][1] = b*R*(-1. - (2.*DXsqrPlusRsqrDX*X)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr))*Sin(phi);
dxyz[3][2] = b*X*(-1. - (2.*DXsqrPlusRsqrDR*R)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr))*Sin(phi);
dxyz[3][3] = b*R*X*Cos(phi)*(-1. + pow2inv(XsqrPlusRsqr));

/* 2nd derivs */
ddxyz[1][1][1] = (b*(-4.*DXsqrPlusRsqrDX*X*XsqrPlusRsqr + Power(XsqrPlusRsqr,4) + (DDXsqrPlusRsqrDXDX*XsqrPlusRsqr - 3.*pow2(DXsqrPlusRsqrDX))*pow2(R) + (-1.*DDXsqrPlusRsqrDXDX*XsqrPlusRsqr + 3.*pow2(DXsqrPlusRsqrDX))*pow2(X) + Power(XsqrPlusRsqr,4)*pow2inv(XsqrPlusRsqr)))/Power(XsqrPlusRsqr,4);
ddxyz[1][1][2] = (b*(2.*DXsqrPlusRsqrDX*R*XsqrPlusRsqr - 2.*DXsqrPlusRsqrDR*X*XsqrPlusRsqr - 3.*DXsqrPlusRsqrDR*DXsqrPlusRsqrDX*pow2(R) + 3.*DXsqrPlusRsqrDR*DXsqrPlusRsqrDX*pow2(X)))/Power(XsqrPlusRsqr,4);
ddxyz[1][1][3] = 0;
ddxyz[1][2][2] = (-1.*b*(-4.*DXsqrPlusRsqrDR*R*XsqrPlusRsqr + Power(XsqrPlusRsqr,4) + (-1.*DDXsqrPlusRsqrDRDR*XsqrPlusRsqr + 3.*pow2(DXsqrPlusRsqrDR))*pow2(R) + (DDXsqrPlusRsqrDRDR*XsqrPlusRsqr - 3.*pow2(DXsqrPlusRsqrDR))*pow2(X) + Power(XsqrPlusRsqr,4)*pow2inv(XsqrPlusRsqr)))/Power(XsqrPlusRsqr,4);
ddxyz[1][2][3] = 0;
ddxyz[1][3][3] = 0;
ddxyz[2][1][1] = (-2.*b*R*Cos(phi)*(2.*DXsqrPlusRsqrDX*XsqrPlusRsqr + DDXsqrPlusRsqrDXDX*X*XsqrPlusRsqr - 3.*X*pow2(DXsqrPlusRsqrDX)))/Power(XsqrPlusRsqr,4);
ddxyz[2][1][2] = (b*Cos(phi)*(6.*DXsqrPlusRsqrDR*DXsqrPlusRsqrDX*R*X - 2.*DXsqrPlusRsqrDR*R*XsqrPlusRsqr - 2.*DXsqrPlusRsqrDX*X*XsqrPlusRsqr - 1.*Power(XsqrPlusRsqr,4) + Power(XsqrPlusRsqr,4)*pow2inv(XsqrPlusRsqr)))/Power(XsqrPlusRsqr,4);
ddxyz[2][1][3] = b*R*(1. + (2.*DXsqrPlusRsqrDX*X)/Power(XsqrPlusRsqr,3) - 1.*pow2inv(XsqrPlusRsqr))*Sin(phi);
ddxyz[2][2][2] = (-2.*b*X*Cos(phi)*(2.*DXsqrPlusRsqrDR*XsqrPlusRsqr + DDXsqrPlusRsqrDRDR*R*XsqrPlusRsqr - 3.*R*pow2(DXsqrPlusRsqrDR)))/Power(XsqrPlusRsqr,4);
ddxyz[2][2][3] = b*X*(1. + (2.*DXsqrPlusRsqrDR*R)/Power(XsqrPlusRsqr,3) - 1.*pow2inv(XsqrPlusRsqr))*Sin(phi);
ddxyz[2][3][3] = -1.*b*R*X*Cos(phi)*(-1. + pow2inv(XsqrPlusRsqr));
ddxyz[3][1][1] = (-2.*b*R*(2.*DXsqrPlusRsqrDX*XsqrPlusRsqr + DDXsqrPlusRsqrDXDX*X*XsqrPlusRsqr - 3.*X*pow2(DXsqrPlusRsqrDX))*Sin(phi))/Power(XsqrPlusRsqr,4);
ddxyz[3][1][2] = (b*(6.*DXsqrPlusRsqrDR*DXsqrPlusRsqrDX*R*X - 2.*DXsqrPlusRsqrDR*R*XsqrPlusRsqr - 2.*DXsqrPlusRsqrDX*X*XsqrPlusRsqr - 1.*Power(XsqrPlusRsqr,4) + Power(XsqrPlusRsqr,4)*pow2inv(XsqrPlusRsqr))*Sin(phi))/Power(XsqrPlusRsqr,4);
ddxyz[3][1][3] = b*R*Cos(phi)*(-1. - (2.*DXsqrPlusRsqrDX*X)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr));
ddxyz[3][2][2] = (-2.*b*X*(2.*DXsqrPlusRsqrDR*XsqrPlusRsqr + DDXsqrPlusRsqrDRDR*R*XsqrPlusRsqr - 3.*R*pow2(DXsqrPlusRsqrDR))*Sin(phi))/Power(XsqrPlusRsqr,4);
ddxyz[3][2][3] = b*X*Cos(phi)*(-1. - (2.*DXsqrPlusRsqrDR*R)/Power(XsqrPlusRsqr,3) + pow2inv(XsqrPlusRsqr));
ddxyz[3][3][3] = -1.*b*R*X*(-1. + pow2inv(XsqrPlusRsqr))*Sin(phi);
}

/* end of: coordtrans_AnsorgXRphi.c */
/* nvars = 0, n* = 210,  n/ = 36,  n+ = 87, n = 333, O = 0 */
