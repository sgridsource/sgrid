/* SpinWeightedSphericalHarmonics.c */
/* Wolfgang Tichy 2/2010 */

#include "sgrid.h"
#include "ModeComputer.h"


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

/* Wigner d function */
double Wigner_d_function(int l, int m, int s, double costheta)
{
    double dWig = 0;
    int C1 = s > m  ? 0 : m-s;
    int C2 = m < -s ? l+m : l-s;
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
// ??? is this correct
    return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
// should it be
//  return (sqrt( fact(l+m) * fact(l-m) /(fact(l+s) * fact(l-s)) ) * dWig); 
}


/* real part of spin-weighted spherical harmonic sYlm */
double Re_sYlm(int l, int m, int s, double costheta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
//if(m<0) m=-m;
  return c*Wigner_d_function(l, m, s, costheta) * cos(m*phi);
}

/* imaginary part of spin-weighted spherical harmonic sYlm */
double Im_sYlm(int l, int m, int s, double costheta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
//if(m<0) m=-m;

  return c*Wigner_d_function(l, m, s, costheta) * sin(m*phi);
}
