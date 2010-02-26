/* SpinWeightedSphericalHarmonics.c */
/* Wolfgang Tichy 2/2010 */

#include "sgrid.h"
#include "ModeComputer.h"

/*
compared to arXiv:0709.0093v2 [gr-qc] at http://arxiv.org/abs/0709.0093
and also arXiv:gr-qc/0610128, there are some differences here:
-We call usually call Wigner_d_function with arg s=2 not -2
hmmm???
*/

/* factorial */
double fact(double n)
{
  if(n<=1)
    return 1.0;
  else
  {
    n=n*fact(n-1);
    return n;
  }
}

/* Wigner d function, coded by Pablo Galaviz */
double Wigner_d_function_PG(int l, int m, int s, double theta)
{
    double dWig = 0;
    int C1 = s > m  ? 0 : m-s;
    int C2 = m < -s ? l+m : l-s;
    double costheta = cos(theta);
    double coshtheta = sqrt(0.5*(1.0 + costheta));
    double sinhtheta = sqrt(0.5*(1.0 - costheta));
    int c1 = C1 < C2 ? C1 : C2;
    int c2 = C1 < C2 ? C2 : C1;
    int t;

    if(c1%2==0)
    {
      for( t = c1; t <= c2; t+=2)
        dWig += pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );
      for( t = c1+1; t <= c2; t+=2)
        dWig -= pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );
    }
    else
    {
      for( t = c1; t <= c2; t+=2)
        dWig -= pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );
      for( t = c1+1; t <= c2; t+=2)
        dWig += pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );
    }

    return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
}

/* Wigner d function, coded by WT */
double Wigner_d_function_WT(int l, int m, int s, double theta)
{
  double Wigd = 0;
  int k1 = s > m  ? 0 : m-s;
  int k2 = m < -s ? l+m : l-s;
  double costhetao2 = cos(theta*0.5);
  double sinthetao2 = sin(theta*0.5);
  int k;

  if(k1%2==0)  /* k1 even */
  {
    for(k = k1; k <= k2; k+=2)
      Wigd += pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
    for(k = k1+1; k <= k2; k+=2)
      Wigd -= pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
  }
  else
  {
    for(k = k1; k <= k2; k+=2)
      Wigd -= pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
    for(k = k1+1; k <= k2; k+=2)
      Wigd += pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
  }
  return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * Wigd);
}

/* real part of spin-weighted spherical harmonic sYlm */
double Re_sYlm(int l, int m, int s, double theta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
  if(s%2 != 0) c=-c; /* multiply by (-1)^s */
  return c*Wigner_d_function_WT(l, m, s, theta) * cos(m*phi);
}

/* imaginary part of spin-weighted spherical harmonic sYlm */
double Im_sYlm(int l, int m, int s, double theta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
  if(s%2 != 0) c=-c; /* multiply by (-1)^s */
  return c*Wigner_d_function_WT(l, m, s, theta) * sin(m*phi);
}
