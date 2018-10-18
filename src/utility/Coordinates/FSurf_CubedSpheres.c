/* FSurf_CubedSpheres.c */
/* Wolfgang Tichy, Oct 2018 */
/* functions to to get sigma on surface of cubed sphere */

#include "sgrid.h"
#include "Coordinates.h"


/* global vars in this file */
int isigma01_co; /* index of Coordinates_CubedSphere_sigma01_co */
int lmax;        /* max l in Ylm expansion */

/* funcs in this file */
double FSurf_CubSph_sigma01_func(tBox *box, int si, double A, double B);



/* return value of surface function sigma01 */
double FSurf_CubSph_sigma01_func(tBox *box, int si, double A, double B)
{
  int N = box->nnodes;
  int offset = si ? N/2 : 0;  /* offset for coeffs in co */
  int ijk, l,m;
  double *co = box->v[isigma01_co];
  double fv, Theta,Phi;
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);

  /* get Theta,Phi from A,B */
  ThetaPhi_of_AB_CubSph(box, A,B, &Theta,&Phi);

  /* make tables of Ylm at Theta,Phi */
  set_YlmTabs(lmax, Theta,Phi, ReYtab, ImYtab);

  /* get func val fv at Theta,Phi */
  fv = 0.;
  ijk=offset;
  for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      double Rco, Ico, Re_Ylm, Im_Ylm;
      /* get real and imag part of coeffs in co */
      Rco = co[ijk++];
      Ico = co[ijk++];

      /* get Ylm at Theta,Phi */
      Ylm_from_Tabs(lmax, ReYtab, ImYtab, l,m, &Re_Ylm,&Im_Ylm);

      /* There is a choice of sign here: define the inner product by
         (f,g) = int f^* g
         and define
         psi_mode = (Ylm, psi)  */
      /* fv = \sum c_lm Ylm
            = \sum (  Re_c_lm Re_Ylm -   Im_c_lm Im_Ylm +
                    i Re_c_lm Im_Ylm + i Im_c_lm Re_Ylm   )
         We assume that fv is real: */
      fv += Rco*Re_Ylm - Ico*Im_Ylm;
    }

  free(ImYtab);
  free(ReYtab);
  return fv;
}


/* put the sYlm in two variables on the grid */
/*
If we have output at the points of SphericalDF in sgrid I can use
We can use spec_SurfaceIntegral to compute surface integrals.
For each mode coeff we need, we could add one radial point.
E.g. make vars for Re and Im where at
i=0     we put rY_0^{0}  + iY_0^{0}
i=1     we put rY_1^{-1} + iY_1^{-1}
i=2     we put rY_1^{0}  + iY_1^{0}
i=3     we put rY_1^{1}  + iY_1^{1}
i=4     we put rY_2^{-2} + iY_2^{-2}
i=5     we put rY_2^{-1} + iY_2^{-1}
i=6     we put rY_2^{0}  + iY_2^{0}
i=7     we put rY_2^{1}  + iY_2^{1}
Then I can use spec_SurfaceIntegral over these vars to compute the real part
of all coeffs.
*/
/* NOTE: Re_Ylmp,Im_Ylmp have size N1*n2*n3 */                         
int FSurf_CubSph_set_Ylm(tBox *box, int N1, double *Re_Ylmp, double *Im_Ylmp,
                         int lmax)
{
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *Yp = box->v[Ind("Y")];
  double *Zp = box->v[Ind("Z")];
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);
  int nYs = (lmax*(lmax+1))/2 + lmax+1; /* number of Ylm's we use */
  int l,m, i,j,k, ijk, Ijk;

  if(nYs>N1) errorexit("size of Re_Ylmp,Im_Ylmp arrays is too small");

  /* loop over l and m and set Ylm over surface. Put each Ylm at a 
     different radial coord for each l,m */
  //printf("setting Ylm in box%d\n", box->b);
  for(k=0; k<n3; k++)
  for(j=0; j<n2; j++)
  {
    double A,B, theta,phi, Re_Ylm,Im_Ylm;

    /* point at i=0 */
    i=0;
    ijk=Index(i,j,k);

    /* get A,B at point ijk */
    A = Yp[ijk];
    B = Zp[ijk];

    /* get theta,phi from A,B */
    ThetaPhi_of_AB_CubSph(box, A,B, &theta,&phi);

    /* make tables of Ylm at Theta,Phi */
    set_YlmTabs(lmax, theta,phi, ReYtab, ImYtab);

    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      /* get Ylm at theta,phi */
      Ylm_from_Tabs(lmax, ReYtab, ImYtab, l,m, &Re_Ylm,&Im_Ylm);

      /* set spherical harmonic Ylm at point ijk. 
         NOTE: Re_Ylmp and Im_Ylmp may not be on the grid and thus have a
         different range for the index i */
      Ijk=Ind_n1n2(i,j,k, N1,n2);
      Re_Ylmp[Ijk] = Re_Ylm;
      Im_Ylmp[Ijk] = Im_Ylm;
      i++;
    }
  }
    
  free(ImYtab);
  free(ReYtab);
  /* return total number of coeffs up to l=lmax */
  return i;
}


/* Compute coeffs of var that has real part at varindex Re_vind and imag.
   part at Im_vind. Do integral over surface with index i=s in X-dir.
   Put coeffs into var with index co_ind.
   If var has zero imag. part set Im_vind=-1. */
int FSurf_CubSph_set_Ylm_coeffs(tBox *box, int s, int Re_vind, int Im_vind,
                                int lmax, int co_ind)
{
  int l,m, i,j,k, ijk, Ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  int nYs = (lmax*(lmax+1))/2 + lmax+1; /* number of Ylm's we use */
  int N1;                               /* X-dim of Re_Ylmp,Im_Ylmp arrays */
  int offset; /* offset used to write into var co_ind */
  double *Re_varp = box->v[Re_vind];
  double *Im_varp;
  double *Re_Ylmp;
  double *Im_Ylmp;
  double *Re_coeffp;
  double *Im_coeffp;
  double *co =  box->v[co_ind];

  //printf("FSurf_CubSph_set_Ylm_coeffs: Computing coeffs\n");
  //printf("lmax=%d\n", lmax);

  /* do we have imag. part in our var? */
  if(Im_vind>0) Im_varp = box->v[Im_vind];
  else          Im_varp = NULL;

  /* make room for all the Ylm's */
  N1 = n1; /* for now because spec_SurfaceIntegral only works then */
  Re_Ylmp = calloc(N1*n2*n3, sizeof(double));
  Im_Ylmp = calloc(N1*n2*n3, sizeof(double));
  Re_coeffp = calloc(N1*n2*n3, sizeof(double));
  Im_coeffp = calloc(N1*n2*n3, sizeof(double));
  if(Re_Ylmp || Im_Ylmp || Re_coeffp || Im_coeffp)
    errorexit("out of memory for Re_Ylmp, Im_Ylmp, ...");

  /* precompute the Ylm */
  FSurf_CubSph_set_Ylm(box, N1, Re_Ylmp, Im_Ylmp, lmax);
  if(N1!=n1) errorexit("implement case where N1!=n1");


  /* set integrands */
  for(k=0; k<n3; k++)
  for(j=0; j<n2; j++)
  {
    i=0;
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      double R,I, RYlm,IYlm;

      /* get Re and Im part of data at surface where i=s */
      ijk=Index(s,j,k);
      R = Re_varp[ijk];
      if(Im_vind<=0) I = 0.; /* imag. part is zero */
      else           I = Im_varp[ijk];
      
      /* get spherical harmonic Ylm */
      Ijk=Ind_n1n2(i,j,k, N1,n2);
      RYlm = Re_Ylmp[Ijk];
      IYlm = Im_Ylmp[Ijk];

      /* There is a choice of sign here: define the inner product by
         (f,g) = int f^* g
         and define
         psi_coeff = (Y, psi)  */
      Re_coeffp[Ijk] = RYlm * R + IYlm * I;
      Im_coeffp[Ijk] = RYlm * I - IYlm * R;
      i++;
    }
  }

  /* integrate over surfaces */
  spec_SurfaceIntegral(box,0, 1, Re_coeffp, Re_coeffp);
  spec_SurfaceIntegral(box,0, 1, Im_coeffp, Im_coeffp);

  /* Put coeffs into var with index co_ind */
  offset = ((box->nnodes/2)*s)/(n1-1);  /* offset for coeffs in co */
  ijk = offset;
  i = 0;
  for(l=0; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* set Re and Im part of coeff */
    co[ijk++] = Re_coeffp[i];
    co[ijk++] = Im_coeffp[i];
    i++;
  }

  free(Im_coeffp);
  free(Re_coeffp);
  free(Im_Ylmp);
  free(Re_Ylmp);
  return 0;
}


/* set var box->CI->iSurf and its derivs from FSurf_CubSph_sigma01_func */
int FSurf_CubSph_set_sigma01vars_from_sigma01_func(tBox *box, int si)
{
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *Yp = box->v[Ind("Y")];
  double *Zp = box->v[Ind("Z")];
  int i,j,k, ijk;
  double A,B, F;
  int isigma    = box->CI->iSurf[si];
  int isigma_dA = box->CI->idSurfdX[si][2];
  int isigma_dB = box->CI->idSurfdX[si][3];
  double *sigma = box->v[isigma];

  i = si*(n1-1);
  for(k=0; k<n3; k++)
  for(j=0; j<n2; j++)
  {
    /* get A,B at point ijk */
    ijk=Index(i,j,k);
    A = Yp[ijk];
    B = Zp[ijk];

    /* set sigma01 var */
    F = box->CI->FSurf[si](box,si, A,B);
    sigma[ijk] = F;
  }

  /* now set derivs */
  compute_CubedSphere_dsigma01(box, isigma, isigma_dA, isigma_dB);
  /* NOTE: Here we could also use spectral derivs of our Ylm expansions!!! */
}

/* initialize function FSurf_CubSph_sigma01_func and its coeffs in 
   isigma01_co of FSurf. Get coeffs from integrating over var in
   box->CI->iFS[si]. */
int FSurf_CubSph_init_from_CI_iFS(tBox *box, int si)
{
  int n1 = box->n1;
  int s  = si ? n1-1 : 0;  /* i-index of surface */
  int iFS = box->CI->iFS[si];
  int ret;

  /* set lmax we use */
  /* We need (lmax*(lmax+1))/2 + lmax+1  complex numbers at each point A,B
     to store the table.
     So when is (lmax*(lmax+1))/2 + lmax+1 = n1?
     set L = lmax ==> L^2/2 + 3L/2 + 1 = n1  <==> L^2 + 3 L + 2 - 2*n1 = 0
     so: 2L = -3 +- sqrt(9 - 4*(2 - 2*n1)) = -3 +- sqrt(8*n1 + 1) 
         L = (sqrt(8*n1 + 1) - 3)/2  */
  lmax = 0.5*(sqrt(8.*n1 + 1.) - 3.);

  /* save var indices */ 
  isigma01_co        = Ind("Coordinates_CubedSphere_sigma01_co");
  box->CI->FSurf[si] = FSurf_CubSph_sigma01_func;

  /* set coeffs co from values of sigma in FS */
  ret=FSurf_CubSph_set_Ylm_coeffs(box, s, iFS,-1, lmax, isigma01_co);

  /* set var box->CI->iSurf and its derivs */
  ret=FSurf_CubSph_set_sigma01vars_from_sigma01_func(box, si);

  return ret;
}