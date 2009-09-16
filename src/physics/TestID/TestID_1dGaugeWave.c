/* TestID_1dGaugeWave.c */
/* Copyright (C) 2009 Wolfgang Tichy */

#include "sgrid.h"
#include "TestID.h"


void TestID_1dGaugeWave(tGrid *grid, int i_x, int i_gb, int i_K, int i_psi, 
                       int i_dpsiopsi, int i_ddpsiopsi, 
                       int i_alpha, int i_beta)
{
  int bi;
  double t = 0.0;
  double A = Getd("TestID_amplitude");

  forallboxes(grid, bi)
  {
    tBox *box = grid->box[bi];
    int i;
    double Lx = box->bbox[1] - box->bbox[0];

    forallpoints(box, i)
    {
      double *x1 = box->v[i_x+0];
      double *x2 = box->v[i_x+1];
      double *x3 = box->v[i_x+2];
      double *gb11 = box->v[i_gb+0];
      double *gb12 = box->v[i_gb+1];
      double *gb13 = box->v[i_gb+2];
      double *gb22 = box->v[i_gb+3];
      double *gb23 = box->v[i_gb+4];
      double *gb33 = box->v[i_gb+5];
      double *K11 = box->v[i_K+0];
      double *K12 = box->v[i_K+1];
      double *K13 = box->v[i_K+2];
      double *K22 = box->v[i_K+3];
      double *K23 = box->v[i_K+4];
      double *K33 = box->v[i_K+5];
      double *psi = box->v[i_psi+0];
      double *dpsiopsi1 = box->v[i_dpsiopsi+0];
      double *dpsiopsi2 = box->v[i_dpsiopsi+1];
      double *dpsiopsi3 = box->v[i_dpsiopsi+2];
      double *ddpsiopsi11 = box->v[i_ddpsiopsi+0];
      double *ddpsiopsi12 = box->v[i_ddpsiopsi+1];
      double *ddpsiopsi13 = box->v[i_ddpsiopsi+2];
      double *ddpsiopsi22 = box->v[i_ddpsiopsi+3];
      double *ddpsiopsi23 = box->v[i_ddpsiopsi+4];
      double *ddpsiopsi33 = box->v[i_ddpsiopsi+5];
      double *alpha = box->v[i_alpha+0];
      double *beta1 = box->v[i_beta+0];
      double *beta2 = box->v[i_beta+1];
      double *beta3 = box->v[i_beta+2];
      /*
      double r2 = x1[i] * x1[i] + x2[i] * x2[i] + x3[i] * x3[i];
      double r = sqrt(r2);
      */
      double arg = 2.0*PI*(x1[i]-t)/Lx;
      double Ac = A*cos(arg);
      double As = A*sin(arg);

      /* conformal factor */
      psi[i] = 1.0;
      dpsiopsi1[i] = dpsiopsi2[i] = dpsiopsi3[i] = 0.0;
      ddpsiopsi11[i] = ddpsiopsi12[i] = ddpsiopsi13[i] = 0.0; 
      ddpsiopsi22[i] = ddpsiopsi23[i] = ddpsiopsi33[i] = 0.0;

      /**************************************************/
      /* set data as in gr-qc/0310100 and gr-qc/0609047 */
      /**************************************************/
      /* conformal metric */
      gb11[i] = 1.0 + As;
      gb22[i] = gb33[i] = 1.0;
      gb12[i] = gb13[i] = gb23[i] = 0.0;

      /* extrinsic curv */
      K12[i] = K13[i] = K22[i] = K23[i] = K33[i] = 0.0;      
      K11[i] = +PI*Ac/( Lx * sqrt(1.0+As) );
      /* contrary to gr-qc/0609047, we choose + in K11 to make it agree with 
         ds^2 = -(1+a)dt^2 + (1+a)dx^2 + dy^2 + dz^2,   a = A sin[2(x-t)/L]
         which moves to the right */

      /* lapse and shift */
      alpha[i] = sqrt(1.0+As);
      beta1[i] = beta2[i] = beta3[i] = 0.0;

//printf("x1[i]=%g Lx=%g gb11[i]=%g K11[i]=%g alpha[i]=%g\n",
//x1[i],Lx,gb11[i],K11[i],alpha[i]);

    } /* end of points */
  } /* end of boxes */
}  /* end of function */
