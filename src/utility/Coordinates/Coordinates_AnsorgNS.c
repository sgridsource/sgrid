/* Coordinates.c */
/* Wolfgang Tichy 4/2005 */

#include "sgrid.h"
#include "Coordinates.h"

/* define some often used code pieces */
#define CALL_ALL_AnsorgNS(dom) \
  tBox *box = (tBox *) aux;\
  double xyz[4];\
  double dABphi[4][4];\
  double ddABphi[4][4][4];\
  xyz_dABphi_ddABphi_of_AnsorgNS(box, ind, dom,A,B,phi, xyz, dABphi, ddABphi);

#define RETURN_x(dom) CALL_ALL_AnsorgNS(dom); return xyz[1]
#define RETURN_y(dom) CALL_ALL_AnsorgNS(dom); return xyz[2]
#define RETURN_z(dom) CALL_ALL_AnsorgNS(dom); return xyz[3]
#define RETURN_dAdx(dom) CALL_ALL_AnsorgNS(dom); return dABphi[1][1]
#define RETURN_dAdy(dom) CALL_ALL_AnsorgNS(dom); return dABphi[1][2]
#define RETURN_dAdz(dom) CALL_ALL_AnsorgNS(dom); return dABphi[1][3]

#define RETURN_dBdx(dom) CALL_ALL_AnsorgNS(dom); return dABphi[2][1]
#define RETURN_dBdy(dom) CALL_ALL_AnsorgNS(dom); return dABphi[2][2]
#define RETURN_dBdz(dom) CALL_ALL_AnsorgNS(dom); return dABphi[2][3]

#define RETURN_dphidx(dom) CALL_ALL_AnsorgNS(dom); return dABphi[3][1]
#define RETURN_dphidy(dom) CALL_ALL_AnsorgNS(dom); return dABphi[3][2]
#define RETURN_dphidz(dom) CALL_ALL_AnsorgNS(dom); return dABphi[3][3]

#define RETURN_ddAdxdx(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][1][1]
#define RETURN_ddAdxdy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][1][2]
#define RETURN_ddAdxdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][1][3]
#define RETURN_ddAdydy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][2][2]
#define RETURN_ddAdydz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][2][3]
#define RETURN_ddAdzdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[1][3][3]

#define RETURN_ddBdxdx(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][1][1]
#define RETURN_ddBdxdy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][1][2]
#define RETURN_ddBdxdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][1][3]
#define RETURN_ddBdydy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][2][2]
#define RETURN_ddBdydz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][2][3]
#define RETURN_ddBdzdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[2][3][3]

#define RETURN_ddphidxdx(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][1][1]
#define RETURN_ddphidxdy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][1][2]
#define RETURN_ddphidxdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][1][3]
#define RETURN_ddphidydy(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][2][2]
#define RETURN_ddphidydz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][2][3]
#define RETURN_ddphidzdz(dom) CALL_ALL_AnsorgNS(dom); return ddABphi[3][3][3]

/* code pieces to get dx/dX */
#define CALL_INV_AnsorgNS(dom) \
  tBox *box = (tBox *) aux;\
  double xyz[4];\
  double dxyz[4][4];\
  double ddxyz[4][4][4];\
  xyz_dxyz_ddxyz_of_AnsorgNS_ABphi(box, ind, dom,A,B,phi, xyz, dxyz, ddxyz);

#define RETINV_x(dom) CALL_INV_AnsorgNS(dom); return xyz[1]
#define RETINV_y(dom) CALL_INV_AnsorgNS(dom); return xyz[2]
#define RETINV_z(dom) CALL_INV_AnsorgNS(dom); return xyz[3]
#define RETINV_dxdA(dom) CALL_INV_AnsorgNS(dom);   return dxyz[1][1]
#define RETINV_dxdB(dom) CALL_INV_AnsorgNS(dom);   return dxyz[1][2]
#define RETINV_dxdphi(dom) CALL_INV_AnsorgNS(dom); return dxyz[1][3]

#define RETINV_dydA(dom) CALL_INV_AnsorgNS(dom);   return dxyz[2][1]
#define RETINV_dydB(dom) CALL_INV_AnsorgNS(dom);   return dxyz[2][2]
#define RETINV_dydphi(dom) CALL_INV_AnsorgNS(dom); return dxyz[2][3]

#define RETINV_dzdA(dom) CALL_INV_AnsorgNS(dom);   return dxyz[3][1]
#define RETINV_dzdB(dom) CALL_INV_AnsorgNS(dom);   return dxyz[3][2]
#define RETINV_dzdphi(dom) CALL_INV_AnsorgNS(dom); return dxyz[3][3]

#define RETINV_ddxdAdA(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[1][1][1]
#define RETINV_ddxdAdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[1][1][2]
#define RETINV_ddxdAdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[1][1][3]
#define RETINV_ddxdBdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[1][2][2]
#define RETINV_ddxdBdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[1][2][3]
#define RETINV_ddxdphidphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[1][3][3]

#define RETINV_ddydAdA(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[2][1][1]
#define RETINV_ddydAdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[2][1][2]
#define RETINV_ddydAdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[2][1][3]
#define RETINV_ddydBdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[2][2][2]
#define RETINV_ddydBdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[2][2][3]
#define RETINV_ddydphidphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[2][3][3]

#define RETINV_ddzdAdA(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[3][1][1]
#define RETINV_ddzdAdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[3][1][2]
#define RETINV_ddzdAdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[3][1][3]
#define RETINV_ddzdBdB(dom) CALL_INV_AnsorgNS(dom);   return ddxyz[3][2][2]
#define RETINV_ddzdBdphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[3][2][3]
#define RETINV_ddzdphidphi(dom) CALL_INV_AnsorgNS(dom); return ddxyz[3][3][3]


/* Coord. trafo and its derivs for domain 0,1,2,3 */
void xyz_dxyz_ddxyz_of_AnsorgNS_ABphi(tBox *box, int ind, int domain, 
                      double A, double B, double phi,
                      double xyz[4], double dxyz[4][4], double ddxyz[4][4][4])
{
  double XRphi[4];
  double dXRphi_dABphi[4][4];
  double ddXRphi_ddABphi[4][4][4];
  double dxyz_dXRphi[4][4];
  double ddxyz_ddXRphi[4][4][4];

  /* get XRphi(A,B,phi) */
  if(domain==0)
  {
    if(B==0.0)      B = dequaleps;
    else if(B==1.0) B-= dequaleps;
    if(A==1.0) A-= dequaleps;
    XRphi_dXRphi_ddXRphi_of_AnsorgNS0_ABphi(box, ind, A,B,phi, XRphi,
                                            dXRphi_dABphi, ddXRphi_ddABphi);
  }
  else if(domain==1)
  {
    if(B==0.0)      B = dequaleps;
    else if(B==1.0) B-= dequaleps;
    XRphi_dXRphi_ddXRphi_of_AnsorgNS1_ABphi(box, ind, A,B,phi, XRphi,
                                            dXRphi_dABphi, ddXRphi_ddABphi);
  }
  else if(domain==2)
  {
    if(B==0.0)      B = dequaleps;
    else if(B==1.0) B-= dequaleps;
    XRphi_dXRphi_ddXRphi_of_AnsorgNS2_ABphi(box, ind, A,B,phi, XRphi,
                                            dXRphi_dABphi, ddXRphi_ddABphi);
  }
  else if(domain==3)
  {
    if(B==0.0)      B = dequaleps;
    else if(B==1.0) B-= dequaleps;
    if(A==1.0) A-= dequaleps;
    XRphi_dXRphi_ddXRphi_of_AnsorgNS3_ABphi(box, ind, A,B,phi, XRphi,
                                            dXRphi_dABphi, ddXRphi_ddABphi);
  }
  else
    errorexit("xyz_dxyz_ddxyz_of_AnsorgNS0_ABphi: domain must be 0,1,2,3");

  /* get xyz(X,R,phi) */
  xyz_dxyz_ddxyz_of_AnsorgXRphi_XRphi(box, ind, XRphi[1],XRphi[2],XRphi[3],
                                      xyz, dxyz_dXRphi, ddxyz_ddXRphi);
  /* get dxyz=dxyz_dABphi */
  dxdX_from_dxdU_dUdX(dxyz, dxyz_dXRphi, dXRphi_dABphi);

  /* get ddxyz=ddxyz_ddABphi */
  ddxdXdX_from_dxdU_dUdX_ddxdUdU_ddUdXdX(ddxyz, dxyz_dXRphi, dXRphi_dABphi,
                                         ddxyz_ddXRphi, ddXRphi_ddABphi);
}

/* Coord. trafo and its derivs for domain 0,1,2,3, but with extra trafo
   ABphi_dABphi_ddABphi_of_Bfunc_AfBfPhi in coordtrans_Bfunc.m */
void xyz_dxyz_ddxyz_of_AnsorgNS_AfBfPhi(tBox *box, int ind, int domain, 
                      double Af, double Bf, double Phi,
                      double xyz[4], double dxyz[4][4], double ddxyz[4][4][4])
{
  double A,B,phi;
  double ABphi[4];
  double dABphi_dAfBfPhi[4][4];
  double ddABphi_ddAfBfPhi[4][4][4];
  double dxyz_dABphi[4][4];
  double ddxyz_ddABphi[4][4][4];

  /* get ABphi(Af,Bf,Phi) */
  ABphi_dABphi_ddABphi_of_Bfunc_AfBfPhi(box, ind, Af,Bf,Phi, ABphi,
                                        dABphi_dAfBfPhi,ddABphi_ddAfBfPhi);
  A = ABphi[1];
  B = ABphi[2];
  phi = ABphi[3];

  /* get xyz(A,B,phi) */
  xyz_dxyz_ddxyz_of_AnsorgNS_ABphi(box, ind, domain, A,B,phi,
                                   xyz, dxyz_dABphi, ddxyz_ddABphi);
  /* get dxyz=dxyz_dAfBfPhi */
  dxdX_from_dxdU_dUdX(dxyz, dxyz_dABphi, dABphi_dAfBfPhi);

  /* get ddxyz=ddxyz_ddAfBfPhi */
  ddxdXdX_from_dxdU_dUdX_ddxdUdU_ddUdXdX(ddxyz, dxyz_dABphi, dABphi_dAfBfPhi,
                                         ddxyz_ddABphi, ddABphi_ddAfBfPhi);
}


/* Coord. trafo with inverted derivs */
void xyz_dABphi_ddABphi_of_AnsorgNS(tBox *box, int ind, int domain,
                  double A, double B, double phi,
                  double xyz[4], 
                  double dABphi_dxyz[4][4], double ddABphi_ddxyz[4][4][4])
{
  double dxyz_dABphi[4][4];
  double ddxyz_ddABphi[4][4][4];

  /* get xyz, dxyz_dABphi, ddxyz_ddABphi */
  xyz_dxyz_ddxyz_of_AnsorgNS_ABphi(box, ind, domain, A,B,phi,
                                   xyz, dxyz_dABphi, ddxyz_ddABphi);
  /* we can use xyz_dxyz_ddxyz_of_AnsorgNS_AfBfPhi instead of
     xyz_dxyz_ddxyz_of_AnsorgNS_ABphi */
//  xyz_dxyz_ddxyz_of_AnsorgNS_AfBfPhi(box, ind, domain, A,B,phi,
//                                   xyz, dxyz_dABphi, ddxyz_ddABphi);
  /* invert derivs */
  dXdx_from_dxdX(dABphi_dxyz, dxyz_dABphi);
  ddXdxdx_from_dXdx_ddxdXdX(ddABphi_ddxyz, dABphi_dxyz, ddxyz_ddABphi);
}

/*************************************************************************/
/* Start NAnsorgNS? functions                                            */
/*************************************************************************/
/* Coord. trafos for domain 0 */
double   x_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_x(0); }
double   y_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_y(0); }
double   z_of_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_z(0); }

/* first coord. derivs for domain 0 */
double  dA_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdx(0); }
double  dA_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdy(0); }
double  dA_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdz(0); }

double  dB_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdx(0); }
double  dB_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdy(0); }
double  dB_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdz(0); }

double  dphi_dx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidx(0); }
double  dphi_dy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidy(0); }
double  dphi_dz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidz(0); }

/* second coord. derivs for domain 0 */
double  ddA_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdx(0); }
double  ddA_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdy(0); }
double  ddA_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdz(0); }
double  ddA_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydy(0); }
double  ddA_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydz(0); }
double  ddA_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdzdz(0); }

double  ddB_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdx(0); }
double  ddB_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdy(0); }
double  ddB_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdz(0); }
double  ddB_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydy(0); }
double  ddB_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydz(0); }
double  ddB_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdzdz(0); }

double  ddphi_dxdx_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdx(0); }
double  ddphi_dxdy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdy(0); }
double  ddphi_dxdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdz(0); }
double  ddphi_dydy_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydy(0); }
double  ddphi_dydz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydz(0); }
double  ddphi_dzdz_NAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidzdz(0); }

/* Coord. trafos for domain 1 */
double   x_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_x(1); }
double   y_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_y(1); }
double   z_of_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_z(1); }

/* first coord. derivs for domain 1 */
double  dA_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdx(1); }
double  dA_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdy(1); }
double  dA_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdz(1); }

double  dB_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdx(1); }
double  dB_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdy(1); }
double  dB_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdz(1); }

double  dphi_dx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidx(1); }
double  dphi_dy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidy(1); }
double  dphi_dz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidz(1); }

/* second coord. derivs for domain 1 */
double  ddA_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdx(1); }
double  ddA_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdy(1); }
double  ddA_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdz(1); }
double  ddA_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydy(1); }
double  ddA_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydz(1); }
double  ddA_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdzdz(1); }

double  ddB_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdx(1); }
double  ddB_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdy(1); }
double  ddB_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdz(1); }
double  ddB_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydy(1); }
double  ddB_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydz(1); }
double  ddB_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdzdz(1); }

double  ddphi_dxdx_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdx(1); }
double  ddphi_dxdy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdy(1); }
double  ddphi_dxdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdz(1); }
double  ddphi_dydy_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydy(1); }
double  ddphi_dydz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydz(1); }
double  ddphi_dzdz_NAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidzdz(1); }

/* Coord. trafos for domain 2 */
double   x_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_x(2); }
double   y_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_y(2); }
double   z_of_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_z(2); }

/* first coord. derivs for domain 2 */
double  dA_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdx(2); }
double  dA_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdy(2); }
double  dA_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdz(2); }

double  dB_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdx(2); }
double  dB_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdy(2); }
double  dB_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdz(2); }

double  dphi_dx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidx(2); }
double  dphi_dy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidy(2); }
double  dphi_dz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidz(2); }

/* second coord. derivs for domain 2 */
double  ddA_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdx(2); }
double  ddA_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdy(2); }
double  ddA_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdz(2); }
double  ddA_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydy(2); }
double  ddA_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydz(2); }
double  ddA_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdzdz(2); }

double  ddB_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdx(2); }
double  ddB_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdy(2); }
double  ddB_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdz(2); }
double  ddB_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydy(2); }
double  ddB_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydz(2); }
double  ddB_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdzdz(2); }

double  ddphi_dxdx_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdx(2); }
double  ddphi_dxdy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdy(2); }
double  ddphi_dxdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdz(2); }
double  ddphi_dydy_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydy(2); }
double  ddphi_dydz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydz(2); }
double  ddphi_dzdz_NAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidzdz(2); }


/* Coord. trafos for domain 3 */
double   x_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_x(3); }
double   y_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_y(3); }
double   z_of_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_z(3); }

/* first coord. derivs for domain 3 */
double  dA_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdx(3); }
double  dA_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdy(3); }
double  dA_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dAdz(3); }

double  dB_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdx(3); }
double  dB_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdy(3); }
double  dB_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dBdz(3); }

double  dphi_dx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidx(3); }
double  dphi_dy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidy(3); }
double  dphi_dz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_dphidz(3); }

/* second coord. derivs for domain 3 */
double  ddA_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdx(3); }
double  ddA_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdy(3); }
double  ddA_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdxdz(3); }
double  ddA_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydy(3); }
double  ddA_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdydz(3); }
double  ddA_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddAdzdz(3); }

double  ddB_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdx(3); }
double  ddB_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdy(3); }
double  ddB_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdxdz(3); }
double  ddB_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydy(3); }
double  ddB_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdydz(3); }
double  ddB_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddBdzdz(3); }

double  ddphi_dxdx_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdx(3); }
double  ddphi_dxdy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdy(3); }
double  ddphi_dxdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidxdz(3); }
double  ddphi_dydy_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydy(3); }
double  ddphi_dydz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidydz(3); }
double  ddphi_dzdz_NAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETURN_ddphidzdz(3); }

/*************************************************************************/
/* Start IAnsorgNS? functions                                            */
/*************************************************************************/
/* Coord. trafos for domain 0 */
double   x_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_x(0); }
double   y_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_y(0); }
double   z_of_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_z(0); }

/* first coord. derivs for domain 0 */
double  dx_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdA(0); }
double  dx_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdB(0); }
double  dx_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdphi(0); }

double  dy_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydA(0); }
double  dy_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydB(0); }
double  dy_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydphi(0); }

double  dz_dA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdA(0); }
double  dz_dB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdB(0); }
double  dz_dp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdphi(0); }

/* second coord. derivs for domain 0 */
double  ddx_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdA(0); }
double  ddx_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdB(0); }
double  ddx_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdphi(0); }
double  ddx_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdB(0); }
double  ddx_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdphi(0); }
double  ddx_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdphidphi(0); }

double  ddy_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdA(0); }
double  ddy_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdB(0); }
double  ddy_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdphi(0); }
double  ddy_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdB(0); }
double  ddy_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdphi(0); }
double  ddy_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydphidphi(0); }

double  ddz_dAdA_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdA(0); }
double  ddz_dAdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdB(0); }
double  ddz_dAdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdphi(0); }
double  ddz_dBdB_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdB(0); }
double  ddz_dBdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdphi(0); }
double  ddz_dpdp_IAnsorgNS0(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdphidphi(0); }

/* Coord. trafos for domain 1 */
double   x_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_x(1); }
double   y_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_y(1); }
double   z_of_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_z(1); }

/* first coord. derivs for domain 1 */
double  dx_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdA(1); }
double  dx_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdB(1); }
double  dx_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdphi(1); }

double  dy_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydA(1); }
double  dy_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydB(1); }
double  dy_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydphi(1); }

double  dz_dA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdA(1); }
double  dz_dB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdB(1); }
double  dz_dp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdphi(1); }

/* second coord. derivs for domain 1 */
double  ddx_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdA(1); }
double  ddx_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdB(1); }
double  ddx_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdphi(1); }
double  ddx_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdB(1); }
double  ddx_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdphi(1); }
double  ddx_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdphidphi(1); }

double  ddy_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdA(1); }
double  ddy_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdB(1); }
double  ddy_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdphi(1); }
double  ddy_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdB(1); }
double  ddy_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdphi(1); }
double  ddy_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydphidphi(1); }

double  ddz_dAdA_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdA(1); }
double  ddz_dAdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdB(1); }
double  ddz_dAdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdphi(1); }
double  ddz_dBdB_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdB(1); }
double  ddz_dBdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdphi(1); }
double  ddz_dpdp_IAnsorgNS1(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdphidphi(1); }

/* Coord. trafos for domain 2 */
double   x_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_x(2); }
double   y_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_y(2); }
double   z_of_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_z(2); }

/* first coord. derivs for domain 2 */
double  dx_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdA(2); }
double  dx_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdB(2); }
double  dx_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdphi(2); }

double  dy_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydA(2); }
double  dy_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydB(2); }
double  dy_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydphi(2); }

double  dz_dA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdA(2); }
double  dz_dB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdB(2); }
double  dz_dp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdphi(2); }

/* second coord. derivs for domain 2 */
double  ddx_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdA(2); }
double  ddx_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdB(2); }
double  ddx_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdphi(2); }
double  ddx_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdB(2); }
double  ddx_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdphi(2); }
double  ddx_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdphidphi(2); }

double  ddy_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdA(2); }
double  ddy_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdB(2); }
double  ddy_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdphi(2); }
double  ddy_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdB(2); }
double  ddy_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdphi(2); }
double  ddy_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydphidphi(2); }

double  ddz_dAdA_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdA(2); }
double  ddz_dAdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdB(2); }
double  ddz_dAdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdphi(2); }
double  ddz_dBdB_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdB(2); }
double  ddz_dBdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdphi(2); }
double  ddz_dpdp_IAnsorgNS2(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdphidphi(2); }


/* Coord. trafos for domain 3 */
double   x_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_x(3); }
double   y_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_y(3); }
double   z_of_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_z(3); }

/* first coord. derivs for domain 3 */
double  dx_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdA(3); }
double  dx_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdB(3); }
double  dx_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dxdphi(3); }

double  dy_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydA(3); }
double  dy_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydB(3); }
double  dy_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dydphi(3); }

double  dz_dA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdA(3); }
double  dz_dB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdB(3); }
double  dz_dp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_dzdphi(3); }

/* second coord. derivs for domain 3 */
double  ddx_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdA(3); }
double  ddx_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdB(3); }
double  ddx_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdAdphi(3); }
double  ddx_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdB(3); }
double  ddx_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdBdphi(3); }
double  ddx_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddxdphidphi(3); }

double  ddy_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdA(3); }
double  ddy_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdB(3); }
double  ddy_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydAdphi(3); }
double  ddy_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdB(3); }
double  ddy_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydBdphi(3); }
double  ddy_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddydphidphi(3); }

double  ddz_dAdA_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdA(3); }
double  ddz_dAdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdB(3); }
double  ddz_dAdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdAdphi(3); }
double  ddz_dBdB_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdB(3); }
double  ddz_dBdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdBdphi(3); }
double  ddz_dpdp_IAnsorgNS3(void *aux, int ind, double A, double B, double phi)
{ RETINV_ddzdphidphi(3); }
