/* cartesianDerivs.c */
/* Wolfgang Tichy 4/2005 */

#include "sgrid.h"
#include "Coordinates.h"

/* compute du/dx_m  where m=1,2,3 */
/* Note dX_dx[l][m] = dX_l/dx_m where x_m is the mth cartesian coord,
   i.e x_m=(x,y,z) and X_l is e.g. X_l=(rho,phi,z)*/
void cart_partials_old(tBox *box, double *u, double *u1, double *u2, double *u3)
{
  spec_Deriv1(box, 1, u, u1);
  spec_Deriv1(box, 2, u, u2);
  spec_Deriv1(box, 3, u, u3);

  if( box->x_of_X[1] != NULL )
  {
    int i,m, ind;
    double *Xp = box->v[Ind("X")];
    double *Yp = box->v[Ind("Y")];
    double *Zp = box->v[Ind("Z")];
    double *du[4];
    double dv[4];

    du[1] = u1;        du[2] = u2;        du[3] = u3;

    /* compute cartesian derivs at all points */
    forallpoints(box,ind)
    {
      double X = Xp[ind];
      double Y = Yp[ind];
      double Z = Zp[ind];
    
      /* Transform derivs to cartesian coords */
      for(m=1; m<=3; m++)
      {
        dv[m] = 0.0;
          
        for(i=1; i<=3; i++)
          dv[m] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z)  *  du[i][ind];
      }
      /* copy dv into du */
      for(m=1; m<=3; m++)
        du[m][ind] = dv[m];
    } /* end loop over all points */

    /* set derivs at singular points */
    for(m=1; m<=3; m++)
    {
      if( box->Sing_d_dx[m] != NULL )
        box->Sing_d_dx[m]((void *) box,   (void *) u,
                          (void *) du[1], (void *) du[2], (void *) du[3]);
    }
  }
}


/* compute du/dx_m and d^2u/(dx_m dx_n) where m,n=1,2,3 */
/* Note ddX_dxdx[l][m][n] = dX_l/(dx_m dx_n) where x_m is the mth 
   cartesian coord, i.e x_m=(x,y,z) and X_l is e.g. X_l=(rho,phi,z)*/
void cart_partial_all_old(tBox *box, double *u, double *u1, double *u2, double *u3,
                      double *u11, double *u12, double *u13,
                      double *u22, double *u23, double *u33 )
{
  /* compute spectral derivs and use u1,u2,u3, u11,u12,u13,u22,u23,u33
     as temporal storage */
  spec_allDerivs(box, u, u1,u2,u3, u11,u12,u13,u22,u23,u33);

  /* Note: the pointers x_of_X[i] are NULL if X_l = x_l */
  if( box->x_of_X[1] != NULL )
  {
    int i,j,m,n, ind;
    double *Xp = box->v[Ind("X")];
    double *Yp = box->v[Ind("Y")];
    double *Zp = box->v[Ind("Z")];
    double *du[4];
    double *ddu[4][4];
    double dv[4];
    double ddv[4][4];

    du[1] = u1;        du[2] = u2;        du[3] = u3;
    ddu[1][1] = u11;   ddu[1][2] = u12;   ddu[1][3] = u13;  
    ddu[2][1] = u12;   ddu[2][2] = u22;   ddu[2][3] = u23;  
    ddu[3][1] = u13;   ddu[3][2] = u23;   ddu[3][3] = u33;  

    if( box->ddX_dxdx[1][1][1] == NULL )
      errorexit("cart_partial_all: box->ddX_dxdx[i][m][n] = NULL.");

    /* loop over all points */
    forallpoints(box,ind)
    {
      double X = Xp[ind];
      double Y = Yp[ind];
      double Z = Zp[ind];
      
      /* Transform derivs to cartesian coords */
      /* Note:
      d u / dx^m = (dX^i / dx^m)  (d u / dX^i)
      d^2 u / dx^m dx^n =  (dX^i / dx^m dx^n) (d u / dX^i)
                          + (dX^i / dx^m) (d u / dX^i dx^n)
                        =   (dX^i / dx^m dx^n) (d u / dX^i)
                          + (dX^i / dx^m) (dX^j / dx^n) (d u / dX^i dX^j)  */
      for(m=1; m<=3; m++)
      {
        dv[m] = 0.0;
        
        for(i=1; i<=3; i++)
          dv[m] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z)  *  du[i][ind];

        for(n=m; n<=3; n++)
        {
          ddv[m][n] = 0.0;

          for(i=1; i<=3; i++)
            ddv[m][n] += box->ddX_dxdx[i][m][n]((void *) box, ind, X,Y,Z) * du[i][ind];

          for(i=1; i<=3; i++)
            for(j=1; j<=3; j++)
              ddv[m][n] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z) *
                           box->dX_dx[j][n]((void *) box, ind, X,Y,Z) *
                           ddu[i][j][ind];
        }
      }

      /* copy dv and ddv into du and ddu */
      for(m=1; m<=3; m++)
      {
        du[m][ind] = dv[m];
        
        for(n=m; n<=3; n++)
          ddu[m][n][ind] = ddv[m][n];
      }
    } /* end loop over all points */
  } /* end if */
}


/* corrdinate trafos used in cart_partials and cart_partial_all */
double dX_dx(int i, int m, int dXd,
             tBox *box, int ind, double X, double Y, double Z)
{
  if(box->v[dXd])
    return box->v[dXd + 3*(i-1) + m-1][ind];
  else
    return box->dX_dx[i][m]((void *) box, ind, X,Y,Z);
}

double ddX_dxdx(int i, int m, int n, int ddXdd,
                tBox *box, int ind, double X, double Y, double Z)
{
  if(box->v[ddXdd])
  {
    int vi;

    if(m==1)  vi=ddXdd +6*(i-1) +(n-1);
    else      vi=ddXdd +6*(i-1) +m+n-1;

    return box->v[vi][ind];
  }
  else
    return box->ddX_dxdx[i][m][n]((void *) box, ind, X,Y,Z);
}


/* compute du/dx_m  where m=1,2,3 */
/* Note dX_dx[l][m] = dX_l/dx_m where x_m is the mth cartesian coord,
   i.e x_m=(x,y,z) and X_l is e.g. X_l=(rho,phi,z)*/
void cart_partials(tBox *box, double *u, double *u1, double *u2, double *u3)
{
  int dXd = Ind("dXdx");

  spec_Deriv1(box, 1, u, u1);
  spec_Deriv1(box, 2, u, u2);
  spec_Deriv1(box, 3, u, u3);

  if( box->x_of_X[1] != NULL )
  {
    int ind, k;
    double *du[4];

    du[1] = u1;        du[2] = u2;        du[3] = u3;

    /* compute cartesian derivs at all points */
    if(box->v[dXd]) /* use stored coord trafos */
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forallpoints(box,ind)
      {
        int i,m,m_1;
        double dv[4];

        /* Transform derivs to cartesian coords */
        for(m=1; m<=3; m++)
        {
          dv[m] = 0.0;
          m_1 = m-1;
            
          for(i=1; i<=3; i++)
            /* dv[m] += dX_dx(i,m, dXd, box, ind, X,Y,Z) * du[i][ind]; */
            dv[m] += box->v[dXd + 3*(i-1) + m_1][ind] * du[i][ind];
        }
        /* copy dv into du */
        for(m=1; m<=3; m++)
          du[m][ind] = dv[m];
      }
    }
    else  /* compute  coord trafos */
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forallpoints(box,ind)
      {
        int i,m;
        double dv[4];
        double *Xp = box->v[Ind("X")];
        double *Yp = box->v[Ind("Y")];
        double *Zp = box->v[Ind("Z")];
        double X = Xp[ind];
        double Y = Yp[ind];
        double Z = Zp[ind];
      
        /* Transform derivs to cartesian coords */
        for(m=1; m<=3; m++)
        {
          dv[m] = 0.0;
            
          for(i=1; i<=3; i++)
            /* dv[m] += dX_dx(i,m, dXd, box, ind, X,Y,Z) * du[i][ind]; */
            dv[m] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z) * du[i][ind];
        }
        /* copy dv into du */
        for(m=1; m<=3; m++)
          du[m][ind] = dv[m];
      } /* end cartesian derivs at all points */
    }

    /* set derivs at singular points */
    for(k=1; k<=3; k++)
    {
      if( box->Sing_d_dx[k] != NULL )
        box->Sing_d_dx[k]((void *) box,   (void *) u,
                          (void *) du[1], (void *) du[2], (void *) du[3]);
    }
  }
}


/* compute du/dx_m and d^2u/(dx_m dx_n) where m,n=1,2,3 */
/* Note ddX_dxdx[l][m][n] = dX_l/(dx_m dx_n) where x_m is the mth 
   cartesian coord, i.e x_m=(x,y,z) and X_l is e.g. X_l=(rho,phi,z)*/
void cart_partial_all(tBox *box, double *u, double *u1, double *u2, double *u3,
                      double *u11, double *u12, double *u13,
                      double *u22, double *u23, double *u33 )
{
  int dXd = Ind("dXdx");
  int ddXdd = Ind("ddXddxx");

  /* compute spectral derivs and use u1,u2,u3, u11,u12,u13,u22,u23,u33
     as temporal storage */
  spec_allDerivs(box, u, u1,u2,u3, u11,u12,u13,u22,u23,u33);

  /* Note: the pointers x_of_X[i] are NULL if X_l = x_l */
  if( box->x_of_X[1] != NULL )
  {
    int ind;
    double *du[4];
    double *ddu[4][4];

    du[1] = u1;        du[2] = u2;        du[3] = u3;
    ddu[1][1] = u11;   ddu[1][2] = u12;   ddu[1][3] = u13;  
    ddu[2][1] = u12;   ddu[2][2] = u22;   ddu[2][3] = u23;  
    ddu[3][1] = u13;   ddu[3][2] = u23;   ddu[3][3] = u33;  

    if( box->ddX_dxdx[1][1][1] == NULL )
      errorexit("cart_partial_all: box->ddX_dxdx[i][m][n] = NULL.");

    /* compute cartesian derivs at all points */
    if(box->v[ddXdd] && box->v[dXd]) /* use stored coord trafos */
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forallpoints(box,ind) /* loop over all points */
      {
        int i,j, m,n, m_1,n_1, vi;
        double dv[4];
        double ddv[4][4];
        /* Transform derivs to cartesian coords */
        /* Note:
        d u / dx^m = (dX^i / dx^m)  (d u / dX^i)
        d^2 u / dx^m dx^n =  (dX^i / dx^m dx^n) (d u / dX^i)
                            + (dX^i / dx^m) (d u / dX^i dx^n)
                          =   (dX^i / dx^m dx^n) (d u / dX^i)
                            + (dX^i / dx^m) (dX^j / dx^n) (d u / dX^i dX^j)  */
        for(m=1; m<=3; m++)
        {
          dv[m] = 0.0;
          m_1 = m-1;

          for(i=1; i<=3; i++)
            dv[m] += box->v[dXd + 3*(i-1) + m_1][ind] * du[i][ind];
  
          for(n=m; n<=3; n++)
          {
            ddv[m][n] = 0.0;
            n_1 = n-1;
  
            for(i=1; i<=3; i++)
            {
              if(m==1)  vi=ddXdd +6*(i-1) +(n_1);
              else      vi=ddXdd +6*(i-1) +m+n_1;
              /* ddv[m][n] += ddX_dxdx(i,m,n, ddXdd, box, ind, X,Y,Z) * du[i][ind]; */
              ddv[m][n] += box->v[vi][ind] * du[i][ind];
            }
            for(i=1; i<=3; i++)
              for(j=1; j<=3; j++)
                /* ddv[m][n] += dX_dx(i,m, dXd, box, ind, X,Y,Z) *
                             dX_dx(j,n, dXd, box, ind, X,Y,Z) *
                             ddu[i][j][ind]; */
                ddv[m][n] += box->v[dXd + 3*(i-1) + m_1][ind] * 
                             box->v[dXd + 3*(j-1) + n_1][ind] *
                             ddu[i][j][ind];
          }
        }
        /* copy dv and ddv into du and ddu */
        for(m=1; m<=3; m++)
        {
          du[m][ind] = dv[m];
          
          for(n=m; n<=3; n++)
            ddu[m][n][ind] = ddv[m][n];
        }
      } /* end loop over all points */
    }
    else /* compute coord trafos */
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forallpoints(box,ind) /* loop over all points */
      {
        int i,j, m,n;
        double dv[4];
        double ddv[4][4];
        double *Xp = box->v[Ind("X")];
        double *Yp = box->v[Ind("Y")];
        double *Zp = box->v[Ind("Z")];
        double X = Xp[ind];
        double Y = Yp[ind];
        double Z = Zp[ind];
        
        /* Transform derivs to cartesian coords */
        /* Note:
        d u / dx^m = (dX^i / dx^m)  (d u / dX^i)
        d^2 u / dx^m dx^n =  (dX^i / dx^m dx^n) (d u / dX^i)
                            + (dX^i / dx^m) (d u / dX^i dx^n)
                          =   (dX^i / dx^m dx^n) (d u / dX^i)
                            + (dX^i / dx^m) (dX^j / dx^n) (d u / dX^i dX^j)  */
        for(m=1; m<=3; m++)
        {
          dv[m] = 0.0;

          for(i=1; i<=3; i++)
            dv[m] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z) * du[i][ind];
              
          for(n=m; n<=3; n++)
          {
            ddv[m][n] = 0.0;
  
            for(i=1; i<=3; i++)
              /* ddv[m][n] += ddX_dxdx(i,m,n, ddXdd, box, ind, X,Y,Z) * du[i][ind]; */
              ddv[m][n] += box->ddX_dxdx[i][m][n]((void *) box, ind, X,Y,Z) *
                           du[i][ind];

            for(i=1; i<=3; i++)
              for(j=1; j<=3; j++)
                /* ddv[m][n] += dX_dx(i,m, dXd, box, ind, X,Y,Z) *
                             dX_dx(j,n, dXd, box, ind, X,Y,Z) *
                             ddu[i][j][ind]; */
                ddv[m][n] += box->dX_dx[i][m]((void *) box, ind, X,Y,Z) *
                             box->dX_dx[j][n]((void *) box, ind, X,Y,Z) *
                             ddu[i][j][ind];
          }
        }
        /* copy dv and ddv into du and ddu */
        for(m=1; m<=3; m++)
        {
          du[m][ind] = dv[m];
          
          for(n=m; n<=3; n++)
            ddu[m][n][ind] = ddv[m][n];
        }
      } /* end compute cartesian derivs at all points */
    }
  } /* end if */
}


/************************************************************************/
/*  Helper functions to compute Cartesian derivs of certain tensors     */
/************************************************************************/ 

/* compute all first derivs S_{ab,c} of a symmetric tensor S_{ab} in a box */
void FirstDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc)
{
  /* tensor from index */
  double *S11 = box->v[i_Sab+0];
  double *S12 = box->v[i_Sab+1];
  double *S13 = box->v[i_Sab+2];
  double *S22 = box->v[i_Sab+3];
  double *S23 = box->v[i_Sab+4];
  double *S33 = box->v[i_Sab+5];

  /* tensor derivs from index */
  double *dS111 = box->v[i_dSabc+0];
  double *dS112 = box->v[i_dSabc+1];
  double *dS113 = box->v[i_dSabc+2];
  double *dS121 = box->v[i_dSabc+3];
  double *dS122 = box->v[i_dSabc+4];
  double *dS123 = box->v[i_dSabc+5];
  double *dS131 = box->v[i_dSabc+6];
  double *dS132 = box->v[i_dSabc+7];
  double *dS133 = box->v[i_dSabc+8];
  double *dS221 = box->v[i_dSabc+9];
  double *dS222 = box->v[i_dSabc+10];
  double *dS223 = box->v[i_dSabc+11];
  double *dS231 = box->v[i_dSabc+12];
  double *dS232 = box->v[i_dSabc+13];
  double *dS233 = box->v[i_dSabc+14];
  double *dS331 = box->v[i_dSabc+15];
  double *dS332 = box->v[i_dSabc+16];
  double *dS333 = box->v[i_dSabc+17];

  /* compute all derivs in box */
  cart_partials(box, S11, dS111,dS112,dS113);
  cart_partials(box, S12, dS121,dS122,dS123);
  cart_partials(box, S13, dS131,dS132,dS133);
  cart_partials(box, S22, dS221,dS222,dS223);
  cart_partials(box, S23, dS231,dS232,dS233);
  cart_partials(box, S33, dS331,dS332,dS333);
} 

/* compute all first and second derivs S_{ab,c} and S_{ab,cd} 
   of a symmetric tensor S_{ab} in a box */
void FirstAndSecondDerivsOf_Sab(tBox *box, int i_Sab, 
                                int i_dSabc, int i_ddSabcd)
{
  /* first tensor derivs from index */
  double *dS111 = box->v[i_dSabc+0];
  double *dS112 = box->v[i_dSabc+1];
  double *dS113 = box->v[i_dSabc+2];
  double *dS121 = box->v[i_dSabc+3];
  double *dS122 = box->v[i_dSabc+4];
  double *dS123 = box->v[i_dSabc+5];
  double *dS131 = box->v[i_dSabc+6];
  double *dS132 = box->v[i_dSabc+7];
  double *dS133 = box->v[i_dSabc+8];
  double *dS221 = box->v[i_dSabc+9];
  double *dS222 = box->v[i_dSabc+10];
  double *dS223 = box->v[i_dSabc+11];
  double *dS231 = box->v[i_dSabc+12];
  double *dS232 = box->v[i_dSabc+13];
  double *dS233 = box->v[i_dSabc+14];
  double *dS331 = box->v[i_dSabc+15];
  double *dS332 = box->v[i_dSabc+16];
  double *dS333 = box->v[i_dSabc+17];

  /* second tensor derivs from index */
  double *ddS1111 = box->v[i_ddSabcd+0];
  double *ddS1112 = box->v[i_ddSabcd+1];
  double *ddS1113 = box->v[i_ddSabcd+2];
  double *ddS1122 = box->v[i_ddSabcd+3];
  double *ddS1123 = box->v[i_ddSabcd+4];
  double *ddS1133 = box->v[i_ddSabcd+5];

  double *ddS1211 = box->v[i_ddSabcd+6];
  double *ddS1212 = box->v[i_ddSabcd+7];
  double *ddS1213 = box->v[i_ddSabcd+8];
  double *ddS1222 = box->v[i_ddSabcd+9];
  double *ddS1223 = box->v[i_ddSabcd+10];
  double *ddS1233 = box->v[i_ddSabcd+11];

  double *ddS1311 = box->v[i_ddSabcd+12];
  double *ddS1312 = box->v[i_ddSabcd+13];
  double *ddS1313 = box->v[i_ddSabcd+14];
  double *ddS1322 = box->v[i_ddSabcd+15];
  double *ddS1323 = box->v[i_ddSabcd+16];
  double *ddS1333 = box->v[i_ddSabcd+17];

  double *ddS2211 = box->v[i_ddSabcd+18];
  double *ddS2212 = box->v[i_ddSabcd+19];
  double *ddS2213 = box->v[i_ddSabcd+20];
  double *ddS2222 = box->v[i_ddSabcd+21];
  double *ddS2223 = box->v[i_ddSabcd+22];
  double *ddS2233 = box->v[i_ddSabcd+23];

  double *ddS2311 = box->v[i_ddSabcd+24];
  double *ddS2312 = box->v[i_ddSabcd+25];
  double *ddS2313 = box->v[i_ddSabcd+26];
  double *ddS2322 = box->v[i_ddSabcd+27];
  double *ddS2323 = box->v[i_ddSabcd+28];
  double *ddS2333 = box->v[i_ddSabcd+29];

  double *ddS3311 = box->v[i_ddSabcd+30];
  double *ddS3312 = box->v[i_ddSabcd+31];
  double *ddS3313 = box->v[i_ddSabcd+32];
  double *ddS3322 = box->v[i_ddSabcd+33];
  double *ddS3323 = box->v[i_ddSabcd+34];
  double *ddS3333 = box->v[i_ddSabcd+35];

  FirstDerivsOf_Sab(box, i_Sab, i_dSabc);
  
  /* compute all second derivs in box */
  cart_partials(box, dS111, ddS1111,ddS1112,ddS1113);
  cart_partials(box, dS112, ddS1112,ddS1122,ddS1123);
  cart_partials(box, dS113, ddS1113,ddS1123,ddS1133);

  cart_partials(box, dS121, ddS1211,ddS1212,ddS1213);
  cart_partials(box, dS122, ddS1212,ddS1222,ddS1223);
  cart_partials(box, dS123, ddS1213,ddS1223,ddS1233);

  cart_partials(box, dS131, ddS1311,ddS1312,ddS1313);
  cart_partials(box, dS132, ddS1312,ddS1322,ddS1323);
  cart_partials(box, dS133, ddS1313,ddS1323,ddS1333);

  cart_partials(box, dS221, ddS2211,ddS2212,ddS2213);
  cart_partials(box, dS222, ddS2212,ddS2222,ddS2223);
  cart_partials(box, dS223, ddS2213,ddS2223,ddS2233);

  cart_partials(box, dS231, ddS2311,ddS2312,ddS2313);
  cart_partials(box, dS232, ddS2312,ddS2322,ddS2323);
  cart_partials(box, dS233, ddS2313,ddS2323,ddS2333);

  cart_partials(box, dS331, ddS3311,ddS3312,ddS3313);
  cart_partials(box, dS332, ddS3312,ddS3322,ddS3323);
  cart_partials(box, dS333, ddS3313,ddS3323,ddS3333);
} 


/* compute all first derivs S_{a,b} of a tensor S_{a} in a box */
void FirstDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab)
{
  /* tensor from index */
  double *S1 = box->v[i_Sa+0];
  double *S2 = box->v[i_Sa+1];
  double *S3 = box->v[i_Sa+2];

  /* tensor derivs from index */
  double *dS11 = box->v[i_dSab+0];
  double *dS12 = box->v[i_dSab+1];
  double *dS13 = box->v[i_dSab+2];
  double *dS21 = box->v[i_dSab+3];
  double *dS22 = box->v[i_dSab+4];
  double *dS23 = box->v[i_dSab+5];
  double *dS31 = box->v[i_dSab+6];
  double *dS32 = box->v[i_dSab+7];
  double *dS33 = box->v[i_dSab+8];

  /* compute all derivs in box */
  cart_partials(box, S1, dS11,dS12,dS13);
  cart_partials(box, S2, dS21,dS22,dS23);
  cart_partials(box, S3, dS31,dS32,dS33);
} 

/* compute all first and second derivs S_{a,b} and S_{a,bc} 
   of a tensor S_{a} in a box */
void FirstAndSecondDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc)
{
  /* first tensor derivs from index */
  double *dS11 = box->v[i_dSab+0];
  double *dS12 = box->v[i_dSab+1];
  double *dS13 = box->v[i_dSab+2];
  double *dS21 = box->v[i_dSab+3];
  double *dS22 = box->v[i_dSab+4];
  double *dS23 = box->v[i_dSab+5];
  double *dS31 = box->v[i_dSab+6];
  double *dS32 = box->v[i_dSab+7];
  double *dS33 = box->v[i_dSab+8];

  /* second tensor derivs from index */
  double *ddS111 = box->v[i_ddSabc+0];
  double *ddS112 = box->v[i_ddSabc+1];
  double *ddS113 = box->v[i_ddSabc+2];
  double *ddS122 = box->v[i_ddSabc+3];
  double *ddS123 = box->v[i_ddSabc+4];
  double *ddS133 = box->v[i_ddSabc+5];

  double *ddS211 = box->v[i_ddSabc+6];
  double *ddS212 = box->v[i_ddSabc+7];
  double *ddS213 = box->v[i_ddSabc+8];
  double *ddS222 = box->v[i_ddSabc+9];
  double *ddS223 = box->v[i_ddSabc+10];
  double *ddS233 = box->v[i_ddSabc+11];

  double *ddS311 = box->v[i_ddSabc+12];
  double *ddS312 = box->v[i_ddSabc+13];
  double *ddS313 = box->v[i_ddSabc+14];
  double *ddS322 = box->v[i_ddSabc+15];
  double *ddS323 = box->v[i_ddSabc+16];
  double *ddS333 = box->v[i_ddSabc+17];

  FirstDerivsOf_Sa(box, i_Sa, i_dSab);
  
  /* compute all second derivs in box */
  cart_partials(box, dS11, ddS111,ddS112,ddS113);
  cart_partials(box, dS12, ddS112,ddS122,ddS123);
  cart_partials(box, dS13, ddS113,ddS123,ddS133);

  cart_partials(box, dS21, ddS211,ddS212,ddS213);
  cart_partials(box, dS22, ddS212,ddS222,ddS223);
  cart_partials(box, dS23, ddS213,ddS223,ddS233);

  cart_partials(box, dS31, ddS311,ddS312,ddS313);
  cart_partials(box, dS32, ddS312,ddS322,ddS323);
  cart_partials(box, dS33, ddS313,ddS323,ddS333);
} 


/* compute all first derivs S_{,a} of a tensor S in a box */
void FirstDerivsOf_S(tBox *box, int i_S, int i_dSa)
{
  /* tensor from index */
  double *S = box->v[i_S];

  /* tensor derivs from index */
  double *dS1 = box->v[i_dSa+0];
  double *dS2 = box->v[i_dSa+1];
  double *dS3 = box->v[i_dSa+2];

  /* compute all derivs in box */
  cart_partials(box, S, dS1,dS2,dS3);
} 

/* compute all first and second derivs S_{,a} and S_{,ab} 
   of a tensor S_{a} in a box */
void FirstAndSecondDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab)
{
  /* first tensor derivs from index */
  double *dS1 = box->v[i_dSa+0];
  double *dS2 = box->v[i_dSa+1];
  double *dS3 = box->v[i_dSa+2];

  /* second tensor derivs from index */
  double *ddS11 = box->v[i_ddSab+0];
  double *ddS12 = box->v[i_ddSab+1];
  double *ddS13 = box->v[i_ddSab+2];
  double *ddS22 = box->v[i_ddSab+3];
  double *ddS23 = box->v[i_ddSab+4];
  double *ddS33 = box->v[i_ddSab+5];

  FirstDerivsOf_Sa(box, i_S, i_dSa);
  
  /* compute all second derivs in box */
  cart_partials(box, dS1, ddS11,ddS12,ddS13);
  cart_partials(box, dS2, ddS12,ddS22,ddS23);
  cart_partials(box, dS3, ddS13,ddS23,ddS33);
} 



/* compute all first and second derivs S_{ab,c} and S_{ab,cd} 
   of a symmetric tensor S_{ab} in a box. In this routine the second
   derivs are computed directly from the DD matrix!                  */
void allDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc, int i_ddSabcd)
{
  /* tensor from index */
  double *S11 = box->v[i_Sab+0];
  double *S12 = box->v[i_Sab+1];
  double *S13 = box->v[i_Sab+2];
  double *S22 = box->v[i_Sab+3];
  double *S23 = box->v[i_Sab+4];
  double *S33 = box->v[i_Sab+5];

  /* first tensor derivs from index */
  double *dS111 = box->v[i_dSabc+0];
  double *dS112 = box->v[i_dSabc+1];
  double *dS113 = box->v[i_dSabc+2];
  double *dS121 = box->v[i_dSabc+3];
  double *dS122 = box->v[i_dSabc+4];
  double *dS123 = box->v[i_dSabc+5];
  double *dS131 = box->v[i_dSabc+6];
  double *dS132 = box->v[i_dSabc+7];
  double *dS133 = box->v[i_dSabc+8];
  double *dS221 = box->v[i_dSabc+9];
  double *dS222 = box->v[i_dSabc+10];
  double *dS223 = box->v[i_dSabc+11];
  double *dS231 = box->v[i_dSabc+12];
  double *dS232 = box->v[i_dSabc+13];
  double *dS233 = box->v[i_dSabc+14];
  double *dS331 = box->v[i_dSabc+15];
  double *dS332 = box->v[i_dSabc+16];
  double *dS333 = box->v[i_dSabc+17];

  /* second tensor derivs from index */
  double *ddS1111 = box->v[i_ddSabcd+0];
  double *ddS1112 = box->v[i_ddSabcd+1];
  double *ddS1113 = box->v[i_ddSabcd+2];
  double *ddS1122 = box->v[i_ddSabcd+3];
  double *ddS1123 = box->v[i_ddSabcd+4];
  double *ddS1133 = box->v[i_ddSabcd+5];

  double *ddS1211 = box->v[i_ddSabcd+6];
  double *ddS1212 = box->v[i_ddSabcd+7];
  double *ddS1213 = box->v[i_ddSabcd+8];
  double *ddS1222 = box->v[i_ddSabcd+9];
  double *ddS1223 = box->v[i_ddSabcd+10];
  double *ddS1233 = box->v[i_ddSabcd+11];

  double *ddS1311 = box->v[i_ddSabcd+12];
  double *ddS1312 = box->v[i_ddSabcd+13];
  double *ddS1313 = box->v[i_ddSabcd+14];
  double *ddS1322 = box->v[i_ddSabcd+15];
  double *ddS1323 = box->v[i_ddSabcd+16];
  double *ddS1333 = box->v[i_ddSabcd+17];

  double *ddS2211 = box->v[i_ddSabcd+18];
  double *ddS2212 = box->v[i_ddSabcd+19];
  double *ddS2213 = box->v[i_ddSabcd+20];
  double *ddS2222 = box->v[i_ddSabcd+21];
  double *ddS2223 = box->v[i_ddSabcd+22];
  double *ddS2233 = box->v[i_ddSabcd+23];

  double *ddS2311 = box->v[i_ddSabcd+24];
  double *ddS2312 = box->v[i_ddSabcd+25];
  double *ddS2313 = box->v[i_ddSabcd+26];
  double *ddS2322 = box->v[i_ddSabcd+27];
  double *ddS2323 = box->v[i_ddSabcd+28];
  double *ddS2333 = box->v[i_ddSabcd+29];

  double *ddS3311 = box->v[i_ddSabcd+30];
  double *ddS3312 = box->v[i_ddSabcd+31];
  double *ddS3313 = box->v[i_ddSabcd+32];
  double *ddS3322 = box->v[i_ddSabcd+33];
  double *ddS3323 = box->v[i_ddSabcd+34];
  double *ddS3333 = box->v[i_ddSabcd+35];

  cart_partial_all(box, S11,  dS111,dS112,dS113,
                   ddS1111,ddS1112,ddS1113,ddS1122,ddS1123,ddS1133);
  cart_partial_all(box, S12,  dS121,dS122,dS123,
                   ddS1211,ddS1212,ddS1213,ddS1222,ddS1223,ddS1233);
  cart_partial_all(box, S13,  dS131,dS132,dS133,
                   ddS1311,ddS1312,ddS1313,ddS1322,ddS1323,ddS1333);
  cart_partial_all(box, S22,  dS221,dS222,dS223,
                   ddS2211,ddS2212,ddS2213,ddS2222,ddS2223,ddS2233);
  cart_partial_all(box, S23,  dS231,dS232,dS233,
                   ddS2311,ddS2312,ddS2313,ddS2322,ddS2323,ddS2333);
  cart_partial_all(box, S33,  dS331,dS332,dS333,
                   ddS3311,ddS3312,ddS3313,ddS3322,ddS3323,ddS3333);
} 

/* compute all first and second derivs S_{a,b} and S_{a,bc}
   of a tensor S_{a} in a box. In this routine the second
   derivs are computed directly from the DD matrix!          */
void allDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc)
{
  /* tensor from index */
  double *S1 = box->v[i_Sa+0];
  double *S2 = box->v[i_Sa+1];
  double *S3 = box->v[i_Sa+2];

  /* first tensor derivs from index */
  double *dS11 = box->v[i_dSab+0];
  double *dS12 = box->v[i_dSab+1];
  double *dS13 = box->v[i_dSab+2];
  double *dS21 = box->v[i_dSab+3];
  double *dS22 = box->v[i_dSab+4];
  double *dS23 = box->v[i_dSab+5];
  double *dS31 = box->v[i_dSab+6];
  double *dS32 = box->v[i_dSab+7];
  double *dS33 = box->v[i_dSab+8];

  /* second tensor derivs from index */
  double *ddS111 = box->v[i_ddSabc+0];
  double *ddS112 = box->v[i_ddSabc+1];
  double *ddS113 = box->v[i_ddSabc+2];
  double *ddS122 = box->v[i_ddSabc+3];
  double *ddS123 = box->v[i_ddSabc+4];
  double *ddS133 = box->v[i_ddSabc+5];

  double *ddS211 = box->v[i_ddSabc+6];
  double *ddS212 = box->v[i_ddSabc+7];
  double *ddS213 = box->v[i_ddSabc+8];
  double *ddS222 = box->v[i_ddSabc+9];
  double *ddS223 = box->v[i_ddSabc+10];
  double *ddS233 = box->v[i_ddSabc+11];

  double *ddS311 = box->v[i_ddSabc+12];
  double *ddS312 = box->v[i_ddSabc+13];
  double *ddS313 = box->v[i_ddSabc+14];
  double *ddS322 = box->v[i_ddSabc+15];
  double *ddS323 = box->v[i_ddSabc+16];
  double *ddS333 = box->v[i_ddSabc+17];

  cart_partial_all(box, S1,  dS11,dS12,dS13,
                   ddS111,ddS112,ddS113,ddS122,ddS123,ddS133);
  cart_partial_all(box, S2,  dS21,dS22,dS23,
                   ddS211,ddS212,ddS213,ddS222,ddS223,ddS233);
  cart_partial_all(box, S3,  dS31,dS32,dS33,
                   ddS311,ddS312,ddS313,ddS322,ddS323,ddS333);
} 

/* compute all first and second derivs S_{,a} and S_{,ab}
   of a scalar S in a box. In this routine the second
   derivs are computed directly from the DD matrix!          */
void allDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab)
{
  /* scalar from index */
  double *S = box->v[i_S];

  /* first tensor derivs from index */
  double *dS1 = box->v[i_dSa+0];
  double *dS2 = box->v[i_dSa+1];
  double *dS3 = box->v[i_dSa+2];

  /* second tensor derivs from index */
  double *ddS11 = box->v[i_ddSab+0];
  double *ddS12 = box->v[i_ddSab+1];
  double *ddS13 = box->v[i_ddSab+2];
  double *ddS22 = box->v[i_ddSab+3];
  double *ddS23 = box->v[i_ddSab+4];
  double *ddS33 = box->v[i_ddSab+5];

  cart_partial_all(box, S,  dS1,dS2,dS3,
                   ddS11,ddS12,ddS13,ddS22,ddS23,ddS33);
} 


/* wrapper for 1st and 2nd derivs of Sab */
void D_and_DD_of_Sab(tBox *box, int i_Sab, int i_dSabc, int i_ddSabcd)
{
  if(Getv("Coordinates_useDD","no"))
    FirstAndSecondDerivsOf_Sab(box, i_Sab, i_dSabc, i_ddSabcd);
  else
    allDerivsOf_Sab(box, i_Sab, i_dSabc, i_ddSabcd);
}
/* wrapper for 1st and 2nd derivs of Sa */
void D_and_DD_of_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc)
{
  if(Getv("Coordinates_useDD","no"))
    FirstAndSecondDerivsOf_Sa(box, i_Sa, i_dSab, i_ddSabc);
  else
    allDerivsOf_Sa(box, i_Sa, i_dSab, i_ddSabc);
}
/* wrapper for 1st and 2nd derivs of S */
void D_and_DD_of_S(tBox *box, int i_S, int i_dSa, int i_ddSab)
{
  if(Getv("Coordinates_useDD","no"))
    FirstAndSecondDerivsOf_S(box, i_S, i_dSa, i_ddSab);
  else
    allDerivsOf_S(box, i_S, i_dSa, i_ddSab);
}
