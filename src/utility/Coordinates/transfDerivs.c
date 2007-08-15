/* transfDerivs.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Coordinates.h"

/* compute dX^i/dx^n from dx^n/dX^l and
   d^2 X^i /(dx^j dx^k) from d^2 x^n/(dX^m dX^l) */
/* IDEA:
(dX^i / dx^n) (dx^n / dX^l) = delta^i_l
==>
d/dx^j ((dX^i/dx^n)(dx^n/dX^l)) = 0
[d^2 X^i /(dx^j dx^n)](dx^n/dX^l) + (dX^i/dx^n) d/dx^j (dx^n/dX^l) = 0
[d^2 X^i /(dx^j dx^n)](dx^n/dX^l) + (dX^i/dx^n) (dX^m/dx^j)d/dX^m (dx^n/dX^l) = 0  | *dX^l/dx^k
d^2 X^i /(dx^j dx^k) + (dX^i/dx^n)(dX^m/dx^j)(dX^l/dx^k) d^2 x^n/(dX^m dX^l)=0
*/


/* compute dX^i/dx^n from dx^n/dX^l */
void dXdx_from_dxdX(double dXdx[][], double dxdX[][])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}} 
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33
     Inverse[M]*Det[M] =
       {{-(m23*m32) + m22*m33, m13*m32 - m12*m33, -(m13*m22) + m12*m23},
        {m23*m31 - m21*m33, -(m13*m31) + m11*m33, m13*m21 - m11*m23},
        {-(m22*m31) + m21*m32, m12*m31 - m11*m32, -(m12*m21) + m11*m22}}  */
  double DetM = dxdX[1][1]*dxdX[2][2]*dxdX[3][3] - 
                dxdX[1][3]*dxdX[2][2]*dxdX[3][1] + 
                dxdX[1][2]*dxdX[2][3]*dxdX[3][1] +
                dxdX[1][3]*dxdX[2][1]*dxdX[3][2] - 
                dxdX[1][1]*dxdX[2][3]*dxdX[3][2] - 
                dxdX[1][2]*dxdX[2][1]*dxdX[3][3];
  dXdx[1][1] = (-(dxdX[2][3]*dxdX[3][2]) + dxdX[2][2]*dxdX[3][3])/DetM;
  dXdx[1][2] = (  dxdX[1][3]*dxdX[3][2]  - dxdX[1][2]*dxdX[3][3])/DetM;
  dXdx[1][3] = (-(dxdX[1][3]*dxdX[2][2]) + dxdX[1][2]*dxdX[2][3])/DetM;
  dXdx[2][1] = (  dxdX[2][3]*dxdX[3][1]  - dxdX[2][1]*dxdX[3][3])/DetM;
  dXdx[2][2] = (-(dxdX[1][3]*dxdX[3][1]) + dxdX[1][1]*dxdX[3][3])/DetM;
  dXdx[2][3] = (  dxdX[1][3]*dxdX[2][1]  - dxdX[1][1]*dxdX[2][3])/DetM;
  dXdx[3][1] = (-(dxdX[2][2]*dxdX[3][1]) + dxdX[2][1]*dxdX[3][2])/DetM;
  dXdx[3][2] = (  dxdX[1][2]*dxdX[3][1]  - dxdX[1][1]*dxdX[3][2])/DetM;
  dXdx[3][3] = (-(dxdX[1][2]*dxdX[2][1]) + dxdX[1][1]*dxdX[2][2])/DetM;
}

/* compute d^2 X^i /(dx^j dx^k) from dx^n/dX^l and d^2 x^n/(dX^m dX^l) */
void ddXdxdx_from_dXdx_ddxdXdX(double ddXdxdx[][][],
                               double dXdx[][], double ddxdXdX[][][])
{
  for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
  for(k=j;k<=3;k++)
    for(n=1;n<=3;n++)
    for(m=1;m<=3;m++)
    for(l=m;l<=3;l++)
      ddXdxdx[i][j][k] = -dXdx[i][n]*dXdx[m][j]*dXdx[l][k]*ddxdXdX[n][m][l];
}
