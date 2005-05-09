/* utility.c */
/* Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "ADMvars.h"


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
  double *dS111 = box->v[i_Sabc+0];
  double *dS112 = box->v[i_Sabc+1];
  double *dS113 = box->v[i_Sabc+2];
  double *dS121 = box->v[i_Sabc+3];
  double *dS122 = box->v[i_Sabc+4];
  double *dS123 = box->v[i_Sabc+5];
  double *dS131 = box->v[i_Sabc+6];
  double *dS132 = box->v[i_Sabc+7];
  double *dS133 = box->v[i_Sabc+8];
  double *dS221 = box->v[i_Sabc+9];
  double *dS222 = box->v[i_Sabc+10];
  double *dS223 = box->v[i_Sabc+11];
  double *dS231 = box->v[i_Sabc+12];
  double *dS232 = box->v[i_Sabc+13];
  double *dS233 = box->v[i_Sabc+14];
  double *dS331 = box->v[i_Sabc+15];
  double *dS332 = box->v[i_Sabc+16];
  double *dS333 = box->v[i_Sabc+17];

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
  FirstDerivsOf_Sab(box, i_Sab, i_dSabc);
  
  /* first tensor derivs from index */
  double *dS111 = box->v[i_Sabc+0];
  double *dS112 = box->v[i_Sabc+1];
  double *dS113 = box->v[i_Sabc+2];
  double *dS121 = box->v[i_Sabc+3];
  double *dS122 = box->v[i_Sabc+4];
  double *dS123 = box->v[i_Sabc+5];
  double *dS131 = box->v[i_Sabc+6];
  double *dS132 = box->v[i_Sabc+7];
  double *dS133 = box->v[i_Sabc+8];
  double *dS221 = box->v[i_Sabc+9];
  double *dS222 = box->v[i_Sabc+10];
  double *dS223 = box->v[i_Sabc+11];
  double *dS231 = box->v[i_Sabc+12];
  double *dS232 = box->v[i_Sabc+13];
  double *dS233 = box->v[i_Sabc+14];
  double *dS331 = box->v[i_Sabc+15];
  double *dS332 = box->v[i_Sabc+16];
  double *dS333 = box->v[i_Sabc+17];

  /* second tensor derivs from index */
  double *dS1111 = box->v[i_Sabcd+0];
  double *dS1112 = box->v[i_Sabcd+1];
  double *dS1113 = box->v[i_Sabcd+2];
  double *dS1122 = box->v[i_Sabcd+3];
  double *dS1123 = box->v[i_Sabcd+4];
  double *dS1133 = box->v[i_Sabcd+5];

  double *dS1211 = box->v[i_Sabcd+6];
  double *dS1212 = box->v[i_Sabcd+7];
  double *dS1213 = box->v[i_Sabcd+8];
  double *dS1222 = box->v[i_Sabcd+9];
  double *dS1223 = box->v[i_Sabcd+10];
  double *dS1233 = box->v[i_Sabcd+11];

  double *dS1311 = box->v[i_Sabcd+12];
  double *dS1312 = box->v[i_Sabcd+13];
  double *dS1313 = box->v[i_Sabcd+14];
  double *dS1322 = box->v[i_Sabcd+15];
  double *dS1323 = box->v[i_Sabcd+16];
  double *dS1333 = box->v[i_Sabcd+17];

  double *dS2211 = box->v[i_Sabcd+18];
  double *dS2212 = box->v[i_Sabcd+19];
  double *dS2213 = box->v[i_Sabcd+20];
  double *dS2222 = box->v[i_Sabcd+21];
  double *dS2223 = box->v[i_Sabcd+22];
  double *dS2233 = box->v[i_Sabcd+23];

  double *dS2311 = box->v[i_Sabcd+24];
  double *dS2312 = box->v[i_Sabcd+25];
  double *dS2313 = box->v[i_Sabcd+26];
  double *dS2322 = box->v[i_Sabcd+27];
  double *dS2323 = box->v[i_Sabcd+28];
  double *dS2333 = box->v[i_Sabcd+29];

  double *dS3311 = box->v[i_Sabcd+30];
  double *dS3312 = box->v[i_Sabcd+31];
  double *dS3313 = box->v[i_Sabcd+32];
  double *dS3322 = box->v[i_Sabcd+33];
  double *dS3323 = box->v[i_Sabcd+34];
  double *dS3333 = box->v[i_Sabcd+35];

  /* compute all second derivs in box */
  cart_partials(box, dS111, ddS1111,dS1112,dS1113);
  cart_partials(box, dS112, ddS1112,dS1122,dS1123);
  cart_partials(box, dS113, ddS1113,dS1123,dS1133);

  cart_partials(box, dS121, ddS1211,dS1212,dS1213);
  cart_partials(box, dS122, ddS1212,dS1222,dS1223);
  cart_partials(box, dS123, ddS1213,dS1223,dS1233);

  cart_partials(box, dS131, ddS1311,dS1312,dS1313);
  cart_partials(box, dS132, ddS1312,dS1322,dS1323);
  cart_partials(box, dS133, ddS1313,dS1323,dS1333);

  cart_partials(box, dS221, ddS2211,dS2212,dS2213);
  cart_partials(box, dS222, ddS2212,dS2222,dS2223);
  cart_partials(box, dS223, ddS2213,dS2223,dS2233);

  cart_partials(box, dS231, ddS2311,dS2312,dS2313);
  cart_partials(box, dS232, ddS2312,dS2322,dS2323);
  cart_partials(box, dS233, ddS2313,dS2323,dS2333);

  cart_partials(box, dS331, ddS3311,dS3312,dS3313);
  cart_partials(box, dS332, ddS3312,dS3322,dS3323);
  cart_partials(box, dS333, ddS3313,dS3323,dS3333);
} 



double detg(double g11, double g12, double g13, 
	    double g22, double g23, double g33)
{
  double gginv11, gginv12, gginv13;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  return g11*gginv11 + g12*gginv12 + g13*gginv13;
}



double invg(double g11, double g12, double g13, 
	    double g22, double g23, double g33,
	    double *i11, double *i12, double *i13, 
	    double *i22, double *i23, double *i33)
{
  double detg, gginv11, gginv12, gginv13, gginv22, gginv23, gginv33;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  gginv22 = g11*g33 - g13*g13;
  gginv23 = g12*g13 - g11*g23;
  gginv33 = g11*g22 - g12*g12;
  detg = g11*gginv11 + g12*gginv12 + g13*gginv13;
  *i11 = gginv11/detg;
  *i12 = gginv12/detg;
  *i13 = gginv13/detg;
  *i22 = gginv22/detg;
  *i23 = gginv23/detg;
  *i33 = gginv33/detg;
  return detg;
}




/* compute and save initial value of the trace of K */
int set_K_initial(tGrid *grid)
{
  double *gxx = Ptr(level, "gxx");
  double *gxy = Ptr(level, "gxy");
  double *gxz = Ptr(level, "gxz");
  double *gyy = Ptr(level, "gyy");
  double *gyz = Ptr(level, "gyz");
  double *gzz = Ptr(level, "gzz");

  double *Kxx = Ptr(level, "Kxx");
  double *Kxy = Ptr(level, "Kxy");
  double *Kxz = Ptr(level, "Kxz");
  double *Kyy = Ptr(level, "Kyy");
  double *Kyz = Ptr(level, "Kyz");
  double *Kzz = Ptr(level, "Kzz");

  double *psi = Ptr(level, "psi");
  double *K_initial = PtrEnable(level, "K_initial");

  double ixx, ixy, ixz, iyy, iyz, izz;
  int i,b;

  /* K_initial is zero upon enabling, return if that is what we want */
  if (!Getv("physics", "KerrSchild")) 
    return;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
      
    forallpoints(box, i)
    {
      invg(gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
	   &ixx, &ixy, &ixz, &iyy, &iyz, &izz);
    
      K_initial[i] = pow(psi[i],-4.0) * 
        (ixx*Kxx[i] + iyy*Kyy[i] + izz*Kzz[i] 
         + 2*(ixy*Kxy[i] + ixz*Kxz[i] + iyz*Kyz[i]));
    }
  }
}
