/* coordtrans_AnsorgNS3.c */
/* Copyright (C) 2013 Wolfgang Tichy, 26.3.2013 */
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


extern double (*Coordinates_AnsorgNS_sigmam)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_dsigmam_dphi)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_ddsigmam_dBdB)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_ddsigmam_dBdphi)(tBox *box, int ind, double B, double phi);
extern double (*Coordinates_AnsorgNS_ddsigmam_dphidphi)(tBox *box, int ind, double B, double phi);


/* func with trafo and 1st and 2nd derivs */
void XRphi_dXRphi_ddXRphi_of_AnsorgNS3_ABphi(tBox *box, int ind, double A, double B, double phi, double XRphi[4],double dXRphi[4][4],double ddXRphi[4][4][4])
{
double sigOfBphi = Coordinates_AnsorgNS_sigmam(box, ind, B, phi);
double sigOf1phi = Coordinates_AnsorgNS_sigmam(box, -1, 1.0, phi);
double d1sigOfBphi = Coordinates_AnsorgNS_dsigmam_dB(box, ind, B, phi);
double d2sigOfBphi = Coordinates_AnsorgNS_dsigmam_dphi(box, ind, B, phi);
double d1d1sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dBdB(box, ind, B, phi);
double d1d2sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dBdphi(box, ind, B, phi);
double d2d2sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dphidphi(box, ind, B, phi);
double AbstanhWsBphi = Abstanh(0.25*sigOfBphi, 0.25*PI*B);
double d1AbstanhWsBphi = dAbstanhdx(0.25*sigOfBphi, 0.25*PI*B);
double d2AbstanhWsBphi = dAbstanhdy(0.25*sigOfBphi, 0.25*PI*B);
double d1d1AbstanhWsBphi = ddAbstanhdxdx(0.25*sigOfBphi, 0.25*PI*B);
double d1d2AbstanhWsBphi = ddAbstanhdxdy(0.25*sigOfBphi, 0.25*PI*B);
double d2d2AbstanhWsBphi = ddAbstanhdydy(0.25*sigOfBphi, 0.25*PI*B);
double ArgtanhWsBphi = Argtanh(0.25*sigOfBphi, 0.25*PI*B);
double d1ArgtanhWsBphi = dArgtanhdx(0.25*sigOfBphi, 0.25*PI*B);
double d2ArgtanhWsBphi = dArgtanhdy(0.25*sigOfBphi, 0.25*PI*B);
double d1d1ArgtanhWsBphi = ddArgtanhdxdx(0.25*sigOfBphi, 0.25*PI*B);
double d1d2ArgtanhWsBphi = ddArgtanhdxdy(0.25*sigOfBphi, 0.25*PI*B);
double d2d2ArgtanhWsBphi = ddArgtanhdydy(0.25*sigOfBphi, 0.25*PI*B);
double AbstanhWs1phi = Abstanh(0.25*sigOf1phi, 0.25*PI);
double ArgtanhWs1phi = Argtanh(0.25*sigOf1phi, 0.25*PI);
double AbsCpOfBphi = Sqrt(AbstanhWsBphi);
double ArgCpOfBphi = ArgtanhWsBphi/2.;
double ReCpOfBphi = AbsCpOfBphi*Cos(ArgCpOfBphi);
double ImCpOfBphi = AbsCpOfBphi*Sin(ArgCpOfBphi);
double AbsCpOf1phi = Sqrt(AbstanhWs1phi);
double ArgCpOf1phi = ArgtanhWs1phi/2.;
double ReCpOf1phi = AbsCpOf1phi*Cos(ArgCpOf1phi);
double ImCpOf1phi = AbsCpOf1phi*Sin(ArgCpOf1phi);
double d1AbsCpOfBphi = ((d1AbstanhWsBphi*d1sigOfBphi)/4. + (d2AbstanhWsBphi*PI)/4.)/(2.*Sqrt(AbstanhWsBphi));
double d1ArgCpOfBphi = ((d1ArgtanhWsBphi*d1sigOfBphi)/4. + (d2ArgtanhWsBphi*PI)/4.)/2.;
double d1ReCpOfBphi = d1AbsCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*d1ArgCpOfBphi*Sin(ArgCpOfBphi);
double d1ImCpOfBphi = AbsCpOfBphi*d1ArgCpOfBphi*Cos(ArgCpOfBphi) + d1AbsCpOfBphi*Sin(ArgCpOfBphi);
double d2AbsCpOfBphi = (d1AbstanhWsBphi*d2sigOfBphi)/(8.*Sqrt(AbstanhWsBphi));
double d2ArgCpOfBphi = (d1ArgtanhWsBphi*d2sigOfBphi)/8.;
double d2ReCpOfBphi = d2AbsCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*d2ArgCpOfBphi*Sin(ArgCpOfBphi);
double d2ImCpOfBphi = AbsCpOfBphi*d2ArgCpOfBphi*Cos(ArgCpOfBphi) + d2AbsCpOfBphi*Sin(ArgCpOfBphi);
double d1d1AbsCpOfBphi = ((d1AbstanhWsBphi*d1d1sigOfBphi)/4. + (d1sigOfBphi*((d1d1AbstanhWsBphi*d1sigOfBphi)/4. + (d1d2AbstanhWsBphi*PI)/4.))/4. + (PI*((d1d2AbstanhWsBphi*d1sigOfBphi)/4. + (d2d2AbstanhWsBphi*PI)/4.))/4.)/(2.*Sqrt(AbstanhWsBphi)) - pow2((d1AbstanhWsBphi*d1sigOfBphi)/4. + (d2AbstanhWsBphi*PI)/4.)/(4.*Power(AbstanhWsBphi,1.5));
double d1d1ArgCpOfBphi = ((d1ArgtanhWsBphi*d1d1sigOfBphi)/4. + (d1sigOfBphi*((d1d1ArgtanhWsBphi*d1sigOfBphi)/4. + (d1d2ArgtanhWsBphi*PI)/4.))/4. + (PI*((d1d2ArgtanhWsBphi*d1sigOfBphi)/4. + (d2d2ArgtanhWsBphi*PI)/4.))/4.)/2.;
double d1d1ReCpOfBphi = d1d1AbsCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*Cos(ArgCpOfBphi)*pow2(d1ArgCpOfBphi) - 2*d1AbsCpOfBphi*d1ArgCpOfBphi*Sin(ArgCpOfBphi) - AbsCpOfBphi*d1d1ArgCpOfBphi*Sin(ArgCpOfBphi);
double d1d1ImCpOfBphi = 2*d1AbsCpOfBphi*d1ArgCpOfBphi*Cos(ArgCpOfBphi) + AbsCpOfBphi*d1d1ArgCpOfBphi*Cos(ArgCpOfBphi) + d1d1AbsCpOfBphi*Sin(ArgCpOfBphi) - AbsCpOfBphi*pow2(d1ArgCpOfBphi)*Sin(ArgCpOfBphi);
double d1d2AbsCpOfBphi = -(d1AbstanhWsBphi*d2sigOfBphi*((d1AbstanhWsBphi*d1sigOfBphi)/4. + (d2AbstanhWsBphi*PI)/4.))/(16.*Power(AbstanhWsBphi,1.5)) + ((d1AbstanhWsBphi*d1d2sigOfBphi)/4. + (d1d1AbstanhWsBphi*d1sigOfBphi*d2sigOfBphi)/16. + (d1d2AbstanhWsBphi*d2sigOfBphi*PI)/16.)/(2.*Sqrt(AbstanhWsBphi));
double d1d2ArgCpOfBphi = ((d1ArgtanhWsBphi*d1d2sigOfBphi)/4. + (d1d1ArgtanhWsBphi*d1sigOfBphi*d2sigOfBphi)/16. + (d1d2ArgtanhWsBphi*d2sigOfBphi*PI)/16.)/2.;
double d1d2ReCpOfBphi = d1d2AbsCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*d1ArgCpOfBphi*d2ArgCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*d1d2ArgCpOfBphi*Sin(ArgCpOfBphi) - d1ArgCpOfBphi*d2AbsCpOfBphi*Sin(ArgCpOfBphi) - d1AbsCpOfBphi*d2ArgCpOfBphi*Sin(ArgCpOfBphi);
double d1d2ImCpOfBphi = AbsCpOfBphi*d1d2ArgCpOfBphi*Cos(ArgCpOfBphi) + d1ArgCpOfBphi*d2AbsCpOfBphi*Cos(ArgCpOfBphi) + d1AbsCpOfBphi*d2ArgCpOfBphi*Cos(ArgCpOfBphi) + d1d2AbsCpOfBphi*Sin(ArgCpOfBphi) - AbsCpOfBphi*d1ArgCpOfBphi*d2ArgCpOfBphi*Sin(ArgCpOfBphi);
double d2d2AbsCpOfBphi = (d1AbstanhWsBphi*d2d2sigOfBphi)/(8.*Sqrt(AbstanhWsBphi)) + (d1d1AbstanhWsBphi*pow2(d2sigOfBphi))/(32.*Sqrt(AbstanhWsBphi)) - (pow2(d1AbstanhWsBphi)*pow2(d2sigOfBphi))/(64.*Power(AbstanhWsBphi,1.5));
double d2d2ArgCpOfBphi = (d1ArgtanhWsBphi*d2d2sigOfBphi)/8. + (d1d1ArgtanhWsBphi*pow2(d2sigOfBphi))/32.;
double d2d2ReCpOfBphi = d2d2AbsCpOfBphi*Cos(ArgCpOfBphi) - AbsCpOfBphi*Cos(ArgCpOfBphi)*pow2(d2ArgCpOfBphi) - 2*d2AbsCpOfBphi*d2ArgCpOfBphi*Sin(ArgCpOfBphi) - AbsCpOfBphi*d2d2ArgCpOfBphi*Sin(ArgCpOfBphi);
double d2d2ImCpOfBphi = 2*d2AbsCpOfBphi*d2ArgCpOfBphi*Cos(ArgCpOfBphi) + AbsCpOfBphi*d2d2ArgCpOfBphi*Cos(ArgCpOfBphi) + d2d2AbsCpOfBphi*Sin(ArgCpOfBphi) - AbsCpOfBphi*pow2(d2ArgCpOfBphi)*Sin(ArgCpOfBphi);
/* end var defs */

/* coord transforms */
XRphi[1] = (1. - A)*(-(B*ReCpOf1phi) + ReCpOfBphi) + B*Cos((1. - A)*ArgCpOf1phi + A*PIh);
XRphi[2] = A*(1. - B) + (1. - A)*(-(B*ImCpOf1phi) + ImCpOfBphi) + B*Sin((1. - A)*ArgCpOf1phi + A*PIh);
XRphi[3] = phi;

/* 1st derivs */
dXRphi[1][1] = B*ReCpOf1phi - 1.*ReCpOfBphi + B*(ArgCpOf1phi - 1.*PIh)*Sin(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
dXRphi[1][2] = -1.*(-1. + A)*(d1ReCpOfBphi - 1.*ReCpOf1phi) + Cos(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
dXRphi[1][3] = d2ReCpOfBphi - 1.*A*d2ReCpOfBphi;
dXRphi[2][1] = 1. - 1.*B + B*ImCpOf1phi - 1.*ImCpOfBphi + B*(-1.*ArgCpOf1phi + PIh)*Cos(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
dXRphi[2][2] = d1ImCpOfBphi - 1.*ImCpOf1phi + A*(-1. - 1.*d1ImCpOfBphi + ImCpOf1phi) + Sin(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
dXRphi[2][3] = d2ImCpOfBphi - 1.*A*d2ImCpOfBphi;
dXRphi[3][1] = 0;
dXRphi[3][2] = 0;
dXRphi[3][3] = 1.;

/* 2nd derivs */
ddXRphi[1][1][1] = -1.*B*Cos(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh)*pow2(-1.*ArgCpOf1phi + PIh);
ddXRphi[1][1][2] = -1.*d1ReCpOfBphi + ReCpOf1phi + (ArgCpOf1phi - 1.*PIh)*Sin(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
ddXRphi[1][1][3] = -1.*d2ReCpOfBphi;
ddXRphi[1][2][2] = d1d1ReCpOfBphi - 1.*A*d1d1ReCpOfBphi;
ddXRphi[1][2][3] = d1d2ReCpOfBphi - 1.*A*d1d2ReCpOfBphi;
ddXRphi[1][3][3] = d2d2ReCpOfBphi - 1.*A*d2d2ReCpOfBphi;
ddXRphi[2][1][1] = -1.*B*pow2(-1.*ArgCpOf1phi + PIh)*Sin(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
ddXRphi[2][1][2] = -1. - 1.*d1ImCpOfBphi + ImCpOf1phi + (-1.*ArgCpOf1phi + PIh)*Cos(ArgCpOf1phi - 1.*A*ArgCpOf1phi + A*PIh);
ddXRphi[2][1][3] = -1.*d2ImCpOfBphi;
ddXRphi[2][2][2] = d1d1ImCpOfBphi - 1.*A*d1d1ImCpOfBphi;
ddXRphi[2][2][3] = d1d2ImCpOfBphi - 1.*A*d1d2ImCpOfBphi;
ddXRphi[2][3][3] = d2d2ImCpOfBphi - 1.*A*d2d2ImCpOfBphi;
ddXRphi[3][1][1] = 0;
ddXRphi[3][1][2] = 0;
ddXRphi[3][1][3] = 0;
ddXRphi[3][2][2] = 0;
ddXRphi[3][2][3] = 0;
ddXRphi[3][3][3] = 0;
}

/* end of: coordtrans_AnsorgNS3.c */
/* nvars = 0, n* = 273,  n/ = 64,  n+ = 123, n = 460, O = 0 */
