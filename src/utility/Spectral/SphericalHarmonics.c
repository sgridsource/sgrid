/* SphericalHarmonics.c */
/* find Legendre Polys and Spherical Harmonics */
/* Wolfgang Tichy 10/2018 */

#include "sgrid.h"
#include "Spectral.h"

/* define PI */
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876

#define ind_lm(l,m)  ((l)*((l)+1))/2+(m)


/* use recursion recommended in numrec to get P_l^m for positive m */
double PLegendre(int l, int m, double x)
{
  double fact, Pkm, Pmm, Pmp1_m, somx2;
  int i, k;

  if(m < 0 || m > l || fabs(x) > 1.0) errorexit("Bad arguments in PLegendre");

  /* build P_m^m */
  Pmm = 1.0;
  if(m > 0)
  {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for(i=1; i<=m; i++)
    {
      Pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  /* are we done? */
  if(l == m) return Pmm;

  /* set P_{m+1}^m = x*(2*m+1) P_m^m */
  Pmp1_m = x*(2*m+1)*Pmm;
  /* is this all we need? */
  if(l == (m+1)) return Pmp1_m;

  /* use recursion: 
     P_l^m = ( x*(2*l-1)*P_{l-1}^m - (l+m-1)*P_{l-2}^m )/(l-m)   */
  for(k=m+2; k<=l; k++) /* k runs through all l we have */
  {
    Pkm    = (x*(2*k-1)*Pmp1_m - (k+m-1)*Pmm)/(k-m);
    Pmm    = Pmp1_m;
    Pmp1_m = Pkm;
  } /* now k=l */
  return Pkm;
}

/* alloc array to hold  lmax*(lmax+1)/2 + lmax+1  Plm */
double *alloc_Plm_Tab(int lmax)
{
  double *P = calloc( (lmax*(lmax+1))/2 + lmax+1, sizeof(double));
  return P;
}

/* Find part of P_l^m with l>=m and fill into P array.
   The array P will be filled with P[(l*(l+1))/2 + m] = P_l^m */
void set_PlmTab_for_l_GE_m(int lmax, int m, double x, double *P)
{
  double Pmm, Pmp1_m, Pkm;
  int k;

  if(m > lmax) errorexit("we need m <= lmax");

  /* get start value P_m^m */
  Pmm = PLegendre(m,m, x);
  P[ind_lm(m,m)] = Pmm;

  if(m == lmax) return;

  /* set P_{m+1}^m = x*(2*m+1) P_m^m */
  Pmp1_m = x*(2*m+1)*Pmm;
  P[ind_lm(m+1,m)] = Pmp1_m;

  /* use recursion: 
     P_l^m = ( x*(2*l-1)*P_{l-1}^m - (l+m-1)*P_{l-2}^m )/(l-m)   */
  for(k=m+2; k<=lmax; k++) /* k runs through all l we have */
  {
    Pkm    = (x*(2*k-1)*Pmp1_m - (k+m-1)*Pmm)/(k-m);
    Pmm    = Pmp1_m;
    Pmp1_m = Pkm;
    P[ind_lm(k,m)] = Pkm;
  } /* now k=l */
}

/* set Plm table P */
void set_PlmTab(int lmax, double x, double *Ptab)
{
  int m;
  for(m=0; m<=lmax; m++) set_PlmTab_for_l_GE_m(lmax, m, x, Ptab);
}

/* get P_l^m from table */
double Plm_from_Tab(int lmax, double *Ptab, int l, int m)
{
  return Ptab[ind_lm(l,m)];
}


/* set Ylm tables ReYtab and ImYtab */
void set_YlmTabs(int lmax, double th, double ph, double *ReYtab, double *ImYtab)
{
  int i, l, m;
  double fac, sfac;
  double oo4pi = 1./(4.*PI);

  /* put Plm into ReYtab for now */
  set_PlmTab(lmax, cos(th), ReYtab);

  for(l=0; l<=lmax; l++)
  {
    fac = (2*l+1) * oo4pi;
    for(m=0; m<=l; m++)
    {
      i = ind_lm(l,m);
      /* Y_l^m = sqrt(fac) P_l^m exp(I m phi), i.e.:
         Ytab[i] = sqrt(fac) * Ptab[i] * cexp(I * m*ph); */
      sfac  = ReYtab[i];  /* save Plm at index i */
      sfac *= sqrt(fac);
      ReYtab[i] = sfac * cos(m*ph); /* now Plm is no longer in ReYtab at i */
      ImYtab[i] = sfac * sin(m*ph);
      fac /= ((l+m+1)*(l-m)); 
    }
  }
}

/* get Y_l^m from table. In: ReYtab[],ImYtab[] => Out: ReYlm,ImYlm */
void Ylm_from_Tabs(int lmax, double *ReYtab, double *ImYtab, int l, int m,
                  double *ReYlm, double *ImYlm)
{
  int mm = abs(m);
  int i = ind_lm(l,mm);
  double pm, cc;

  /* take care of case with negative m */
  if(m<0)
  {
    pm = 1 - 2*(mm%2); /* (-1)^m */
    cc = -1.;          /* -1 for complex conj. */
    *ReYlm = pm*ReYtab[i];
    *ImYlm = pm*ImYtab[i]*cc;
  }
  else
  {
    *ReYlm = ReYtab[i];
    *ImYlm = ImYtab[i];
  }
}

/* some slow functions for testing */
double Plm_slowly(int l, int m, double x)
{
  double ret;
  double *Ptab = alloc_Plm_Tab(l);

  set_PlmTab(l, x, Ptab);
  ret = Plm_from_Tab(l, Ptab, l,m);

  free(Ptab);
  return ret;
}
double Re_Ylm_slowly(int l, int m, double theta, double phi)
{
  double R, I;
  double *ReYtab = alloc_Plm_Tab(l);
  double *ImYtab = alloc_Plm_Tab(l);

  set_YlmTabs(l, theta, phi, ReYtab, ImYtab);
  Ylm_from_Tabs(l, ReYtab, ImYtab, l,m, &R,&I);

  free(ImYtab);
  free(ReYtab);
  return R;
}
double Im_Ylm_slowly(int l, int m, double theta, double phi)
{
  double R, I;
  double *ReYtab = alloc_Plm_Tab(l);
  double *ImYtab = alloc_Plm_Tab(l);

  set_YlmTabs(l, theta, phi, ReYtab, ImYtab);
  Ylm_from_Tabs(l, ReYtab, ImYtab, l,m, &R,&I);

  free(ImYtab);
  free(ReYtab);
  return I;
}


/* from physics/ModeComputer/SpinWeightedSphericalHarmonics.c */
double Re_sYlm(int l, int m, int s, double theta, double phi);
double Im_sYlm(int l, int m, int s, double theta, double phi);

/* test: compare with Re_sYlm and Im_sYlm */
int test_SphericalHarmonics(void)
{
  int l,m, lmax=4;
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);

  for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      double theta, phi, R1,I1, R2,I2, dR,dI;
      for(theta = 0.11; theta<3.0; theta+=0.51)
        for(phi = 0.07; phi<6.0; phi+=0.53)
        {
          R1 = Re_Ylm_slowly(l,m, theta,phi);
          I1 = Im_Ylm_slowly(l,m, theta,phi);
          /*
          set_YlmTabs(l, theta, phi, ReYtab, ImYtab);
          Ylm_from_Tabs(l, ReYtab, ImYtab, l,m, &R1,&I1);
          */
          R2 = Re_sYlm(l,m,0, theta,phi);
          I2 = Im_sYlm(l,m,0, theta,phi);
          dR = R1 - R2;
          dI = I1 - I2;
          if(fabs(dR)>1e-10 || fabs(dI)>1e-10)
          {
            printf("l=%d m=%d  theta=%g phi=%g  dR=%g dI=%g\n",
                   l,m, theta,phi, dR,dI);
            printf("l=%d m=%d  theta=%g phi=%g  R1=%g I1=%g\n",
                   l,m, theta,phi, R1,I1);
            printf("l=%d m=%d  theta=%g phi=%g  R2=%g I2=%g\n",
                   l,m, theta,phi, R2,I2);
            /*
            printf("PLegendre=%g\n", PLegendre(l,m, cos(theta)));
            printf("Plm_slowly=%g\n", Plm_slowly(l,m, cos(theta)));
            */
            //return -1;
          }
        }
      printf("l=%d m=%d  final theta,phi  dR=%g dI=%g\n", l,m, dR,dI);
    }
  exit(0);

  free(ImYtab);
  free(ReYtab);
  return 0;
}


/******************************************************************/
/* routines for expansions in Ylm:  f = \sum_l \sum_m c_l^m Y_l^m */
/******************************************************************/

/* NOTE: Here coeffs are stored in this order in memory
   c_0^0, c_1^0,c_1^1, c_2^0,c_2^1,c_2^2, ...
   Each coeffs has a Re and Im part, i.e. in memory we have:
   Rc_0^0, Ic_0^0,  Rc_1^0, Ic_1^0, Rc_1^1, Ic_1^1, ... */
/* Also in some fields people use latitude instead of theta:
   latitude = \pi/2 - \theta
   \partial_{latitude} = -\partial_{\theta} */

/* use Sph. Hamonic coeffs c of a real function to compute coeffs cdphi of
   phi-deriv of function. We use only m>=0, since func is real. */
void SphHarm_dphi_forRealFunc(double *c, double *cdphi, int lmax)
{
  int i, l,m;
  double Rc,Ic, Rcd,Icd; /* Re and Im part of c and cdphi */

  for(i=0, l=0; l<=lmax; l++)
    for(m=0; m<=l; m++, i++) /* here we set only use  m>=0 */
    {
      Rc = c[i];
      Ic = c[i+1];
      Icd =  m*Rc;
      Rcd = -m*Ic;
      cdphi[i]   = Rcd;
      cdphi[i+1] = Icd;
    }
}

/* use Sph. Hamonic coeffs c of a real function to compute coeffs csdth of
   theta-deriv of function. We use only m>=0, since func is real. */
/* NOTE:
     \sin(\theta)\partial_{\theta} Y_l^m
       = l e_{l+1}^m  Y_{l+1}^m  -  (l+1) e_l^m  Y_{l-1}^m
   Thus
     (\sin(\theta)\partial_{\theta} f)_l^m
       = (l-1) e_l^m f_{l-1}^m  -  (l+2) e_{l+1}^m f_{l+1}^m
   where:
     e_0^0 := 0,   e_l^m := \sqrt{ (l^2 - m^2)/((2l+1)(2l-1)) }   */
void SphHarm_sin_theta_dtheta_forRealFunc(double *c, double *csdth, int lmax)
{
  int i, l,m, lp1, lm1, lp1_2, l_2, lp2;
  double elm; /* e_0^0 := 0, e_l^m := sqrt( (l^2 - m^2)/((2l+1)(2l-1)) ) */
  double elp1m; /* e_{l+1}^m */
  double Rclm1,Iclm1, Rclp1,Iclp1; /* Re, Im part of c_{l-1}^m and c_{l+1}^m */
  double Rcld,Icld; /* Re, Im part of csdth */

  for(i=0, l=0; l<=lmax; l++)
    for(m=0; m<=l; m++, i++) /* here we set only use  m>=0 */
    {
      l_2 = l*2;
      elm = l*l - m*m;
      elm = elm/( (l_2 + 1)*(l_2 - 1) );
      elm = sqrt(elm);

      lp1 = l+1;
      lp1_2 = lp1*2;
      elp1m = lp1*lp1 - m*m;
      elp1m = elp1m/( (lp1_2 + 1)*(lp1_2 - 1) );
      elp1m = sqrt(elp1m);

      Rclp1 = c[i+lp1_2];
      Iclp1 = c[i+lp1_2+1];
      Rclm1 = c[i-l_2];
      Iclm1 = c[i-l_2+1];

      lm1  = l - 1;
      lp2  = lp1 + 1;
      Rcld = lm1 * elm * Rclm1 - lp2 * elp1m * Rclp1;
      Icld = lm1 * elm * Iclm1 - lp2 * elp1m * Iclp1;
      csdth[i]   = Rcld;
      csdth[i+1] = Icld;
    }
}
