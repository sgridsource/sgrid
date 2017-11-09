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
void dXdx_from_dxdX(double dXdx[4][4], double dxdX[4][4])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}} 
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33
     Inverse[M]*Det[M] =
       {{-(m23*m32) + m22*m33, m13*m32 - m12*m33, -(m13*m22) + m12*m23},
        {m23*m31 - m21*m33, -(m13*m31) + m11*m33, m13*m21 - m11*m23},
        {-(m22*m31) + m21*m32, m12*m31 - m11*m32, -(m12*m21) + m11*m22}}  */
  LDOUBLE DetM =dxdX[1][1]*dxdX[2][2]*dxdX[3][3] - 
                dxdX[1][3]*dxdX[2][2]*dxdX[3][1] + 
                dxdX[1][2]*dxdX[2][3]*dxdX[3][1] +
                dxdX[1][3]*dxdX[2][1]*dxdX[3][2] - 
                dxdX[1][1]*dxdX[2][3]*dxdX[3][2] - 
                dxdX[1][2]*dxdX[2][1]*dxdX[3][3];
  LDOUBLE sum;
  if(DetM==0.0) DetM=dequaleps*dequaleps*dequaleps;
  sum        = (-(dxdX[2][3]*dxdX[3][2]) + dxdX[2][2]*dxdX[3][3])/DetM;
  dXdx[1][1] = sum;
  sum        = (  dxdX[1][3]*dxdX[3][2]  - dxdX[1][2]*dxdX[3][3])/DetM;
  dXdx[1][2] = sum;
  sum        = (-(dxdX[1][3]*dxdX[2][2]) + dxdX[1][2]*dxdX[2][3])/DetM;
  dXdx[1][3] = sum;
  sum        = (  dxdX[2][3]*dxdX[3][1]  - dxdX[2][1]*dxdX[3][3])/DetM;
  dXdx[2][1] = sum;
  sum        = (-(dxdX[1][3]*dxdX[3][1]) + dxdX[1][1]*dxdX[3][3])/DetM;
  dXdx[2][2] = sum;
  sum        = (  dxdX[1][3]*dxdX[2][1]  - dxdX[1][1]*dxdX[2][3])/DetM;
  dXdx[2][3] = sum;
  sum        = (-(dxdX[2][2]*dxdX[3][1]) + dxdX[2][1]*dxdX[3][2])/DetM;
  dXdx[3][1] = sum;
  sum        = (  dxdX[1][2]*dxdX[3][1]  - dxdX[1][1]*dxdX[3][2])/DetM;
  dXdx[3][2] = sum;
  sum        = (-(dxdX[1][2]*dxdX[2][1]) + dxdX[1][1]*dxdX[2][2])/DetM;
  dXdx[3][3] = sum;
}
/* compute dx^i/dX^n from dX^n/dx^l */
void dxdX_from_dXdx(double dxdX[4][4], double dXdx[4][4])
{
  dXdx_from_dxdX(dxdX, dXdx);
}

/* compute d^2 X^i /(dx^j dx^k) from dx^n/dX^l and d^2 x^n/(dX^m dX^l) */
void ddXdxdx_from_dXdx_ddxdXdX(double ddXdxdx[4][4][4],
                               double dXdx[4][4], double ddxdXdX[4][4][4])
{
  int i,j,k, n,m,l;
  for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
  for(k=j;k<=3;k++)
  {
    LDOUBLE ddXdxdx_ijk = 0.0;
    for(n=1;n<=3;n++)
    for(m=1;m<=3;m++)
    for(l=1;l<=3;l++)
    {
      LDOUBLE ddxdXdX_nml;
      if(l>=m)  ddxdXdX_nml = ddxdXdX[n][m][l];
      else      ddxdXdX_nml = ddxdXdX[n][l][m];
      ddXdxdx_ijk -= dXdx[i][n]*dXdx[m][j]*dXdx[l][k]*ddxdXdX_nml;
    }
    ddXdxdx[i][j][k] = ddXdxdx_ijk;
  }
}


/* 2 coord trafos: X -> U = g(X)  and  U -> x = f(U):
    x^i = f^i(U^m)
    U^j = g^j(X^n)
compute dx^i/dX^m and d^2 x^i/(dX^m dX^n) from 
   dx^i/dU^n , dU^n/dX^m and d^2 x^i/(dU^j dU^k) , d^2 U^j/(dX^m dX^n) */
/* IDEA:
dx^i/dX^m = dx^i/dU^n dU^n/dX^m
d^2 x^i/(dX^m dX^n) = d/dX^m (dU^k/dX^n dx^i/dU^k)
 = d^2 dU^k/(dX^m dX^n) dx^i/dU^k + dU^k/dX^n d^2 x^i/(dU^j dU^k) dU^j/dX^m
 = (dx^i/dU^j) d^2 dU^j/(dX^m dX^n) + dU^j/dX^m dU^k/dX^n d^2 x^i/(dU^j dU^k)
*/

/* dx^i/dX^m = dx^i/dU^n dU^n/dX^m */
void dxdX_from_dxdU_dUdX(double dxdX[4][4],
                         double dxdU[4][4], double dUdX[4][4])
{
  int i,m,n;

  for(i=1;i<=3;i++)
  for(m=1;m<=3;m++)
  {
    LDOUBLE sum=0.0;

    for(n=1;n<=3;n++) sum += dxdU[i][n] * dUdX[n][m];
    dxdX[i][m] = sum;
  }
}

/* d^2 x^i/(dX^m dX^n) = 
  (dx^i/dU^j) d^2 dU^j/(dX^m dX^n) + dU^j/dX^m dU^k/dX^n d^2 x^i/(dU^j dU^k) */
void ddxdXdX_from_dxdU_dUdX_ddxdUdU_ddUdXdX(double ddxdXdX[4][4][4],
                         double dxdU[4][4], double dUdX[4][4],
                         double ddxdUdU[4][4][4], double ddUdXdX[4][4][4])
{
  int i,m,n, j,k;

  for(i=1;i<=3;i++)
  for(m=1;m<=3;m++)
  for(n=m;n<=3;n++)
  {
    LDOUBLE sum1=0.0;
    LDOUBLE sum2=0.0;

    for(j=1;j<=3;j++)
    {
      sum1 += dxdU[i][j] * ddUdXdX[j][m][n];
      for(k=1;k<=3;k++)
      {
        LDOUBLE ddxdUdU_ijk;
        if(k>=j)  ddxdUdU_ijk = ddxdUdU[i][j][k];
        else      ddxdUdU_ijk = ddxdUdU[i][k][j];
        sum2 += dUdX[j][m] * dUdX[k][n] * ddxdUdU_ijk;
      }
    }
    ddxdXdX[i][m][n] = sum1 + sum2;
  }
}


/* check box->dx_dX[k][l] at point in box */
double check_box_dx_dX(tBox *box, double X, double Y, double Z)
{
  double sum, temp, norm;
  int k, l;
  double dXdx[4][4];
  double dxdX[4][4];

  if(box->dx_dX[1][1]==NULL)
  {
    printf("box->dx_dX[1][1]==NULL\n");
    return -1.0;
  }

  /* get dxdX */
  for(k=1; k<=3; k++)
  for(l=1; l<=3; l++)
    dxdX[k][l] = box->dx_dX[k][l](box, -1, X,Y,Z);

  /* compute inverse dXdx from dxdX */
  dXdx_from_dxdX(dXdx, dxdX);

  /* compare dXdx and box->dX_dx[k][l] */
  printf("box->dX_dx[k][l]:\n");
  for(k=1; k<=3; k++)
  {
    for(l=1; l<=3; l++)
      printf("  %+.16e", box->dX_dx[k][l](box, -1, X,Y,Z));
    printf("\n");
  }
  printf("dXdx[k][l], i.e. inverse of box->dx_dX[k][l]:\n");
  for(k=1; k<=3; k++)
  {
    for(l=1; l<=3; l++)
      printf("  %+.16e", dXdx[k][l]);
    printf("\n");
  }

  sum=0.0;
  for(k=1; k<=3; k++)
  for(l=1; l<=3; l++)
  {
    temp = dXdx[k][l] - box->dX_dx[k][l](box, -1, X,Y,Z);
    sum += temp*temp;
  }
  norm = sqrt(sum);
  printf("L2-diff = %g\n", norm);
  return norm;
}
