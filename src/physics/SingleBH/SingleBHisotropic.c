/* SingleBHisotropic.c */
/* Copyright (C) 2005 Wolfgang Tichy */

#include "sgrid.h"
#include "SingleBH.h"


void SingleBHisotropic(tGrid *grid, int i_x, int i_gb, int i_K, int i_psi, 
                       int i_dpsiopsi, int i_ddpsiopsi, 
                       int i_alpha, int i_beta)
{
  int bi;

  double M = Getd("BHmass1");
  int ConformalFactor = Getv("SingleBH_ConformalFactor","yes");

  for(bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int i;

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

      double r2 = x1[i] * x1[i] + x2[i] * x2[i] + x3[i] * x3[i];
      double r = sqrt(r2);
      double psi4;

      /* conformal factor */
      psi[i] = 1 + M/(2.0*r);
   
      /* conformal metric */
      gb12[i] = gb13[i] = gb23[i] = 0.0;
      gb11[i] = gb22[i] = gb33[i] = 1.0;

      /* extrinsic K is zero */
      K11[i] = K12[i] = K13[i] = K22[i] = K23[i] = K33[i] = 0.0;      
      
      /* lapse and shift */
      alpha[i] = ( 1.0 - M/(2.0*r) )/( 1.0 + M/(2.0*r) );
      beta1[i] = beta2[i] = beta3[i] = 0.0;

      if(ConformalFactor)
      {
        double r3 = r2*r;

        dpsiopsi1[i] = (-M/(2.0*r2) * x1[i]/r)/psi[i];
        dpsiopsi2[i] = (-M/(2.0*r2) * x2[i]/r)/psi[i];
        dpsiopsi3[i] = (-M/(2.0*r2) * x3[i]/r)/psi[i];
        
        /* ddpsiopsi[a,b] 
            = 0.25*( (6M/r^3) (x[a]/r) (x[b]/r) -(2M/r^3)delta[a,b] )/psi */
        ddpsiopsi11[i] =( (3*M/(2*r3))*(x1[i]/r)*(x1[i]/r) -(M/(2*r3)) )/psi[i];
        ddpsiopsi12[i] =( (3*M/(2*r3))*(x1[i]/r)*(x2[i]/r) )/psi[i];
        ddpsiopsi13[i] =( (3*M/(2*r3))*(x1[i]/r)*(x3[i]/r) )/psi[i];
        ddpsiopsi22[i] =( (3*M/(2*r3))*(x2[i]/r)*(x2[i]/r) -(M/(2*r3)) )/psi[i];
        ddpsiopsi23[i] =( (3*M/(2*r3))*(x2[i]/r)*(x3[i]/r) )/psi[i];
        ddpsiopsi33[i] =( (3*M/(2*r3))*(x3[i]/r)*(x3[i]/r) -(M/(2*r3)) )/psi[i];
      } 
      else
      { /* if (!ConformalFactor) */
        psi4 = pow(psi[i],4.0);
        gb12[i] = gb13[i] = gb23[i] = 0.0;
        gb11[i] = gb22[i] = gb33[i] = psi4;
              
        psi[i] = 1.0;
        dpsiopsi1[i] = dpsiopsi2[i] = dpsiopsi3[i] = 0.0;
        ddpsiopsi11[i] = ddpsiopsi12[i] = ddpsiopsi13[i] = 0.0;
        ddpsiopsi22[i] = ddpsiopsi23[i] = ddpsiopsi33[i] = 0.0;
      }
      
    } /* end of points */
  } /* end of boxes */
}  /* end of function */
