/* FSurf_CubedSpheres.c */
/* Wolfgang Tichy, Oct 2018 */
/* functions to to get sigma on surface of cubed sphere */

#include "sgrid.h"
#include "Coordinates.h"


/* global vars in this file */
int isigma01_co; /* index of Coordinates_CubedSphere_sigma01_co */
int lmax;        /* max l in Ylm expansion */

/* funcs in this file */
double sigma01_func(tBox *box, int si, double A, double B);


/**/
void init_CubedSphere_FSurf_from_iFS(tBox *box, int si)
{
  int N = box->nnodes;
  int offset = si ? N/2 : 0;
  int n1 = box->n1;
  int i  = si ? n1-1 : 0;
  int iFS = box->CI->iFS[si];
  int ijk,l,m;
  double *co, *FS;

  /* set lmax we use */
  /* We need (lmax*(lmax+1))/2 + lmax+1  double numbers at each point A,B
     to store the table.
     So when is (lmax*(lmax+1))/2 + lmax+1 = n1?
     set L = lmax ==> L^2/2 + 3L/2 + 1 = n1  <==> L^2 + 3 L + 2 - 2*n1 = 0
     so: 2L = -3 +- sqrt(9 - 4*(2 - 2*n1)) = -3 +- sqrt(8*n1 + 1) 
         L = (sqrt(8*n1 + 1) - 3)/2  */
  lmax = 0.5*(sqrt(8.*n1 + 1.) - 3.);

  /* save var indices */ 
  isigma01_co        = Ind("Coordinates_CubedSphere_sigma01_co");
  box->CI->FSurf[si] = sigma01_func;

  /* get pointers */
  co = box->v[isigma01_co];
  FS = box->v[iFS];

  /* set coeffs co from values of sigma in FS */
  ijk=offset;
  for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      co[ijk++] = 0.; // FIXME
    }
  co[0] = FS[0] // FIXME
}

/* return value of surface function sigma01 */
double sigma01_func(tBox *box, int si, double A, double B)
{
  int N = box->nnodes;
  int offset = si ? N/2 : 0;  /* offset for coeffs in co */
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
