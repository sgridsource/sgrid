/* attenuation.c, Wolfgang Tichy 28.12.2003 */
/* useful attenuation functions */


#include <math.h>


/* Attenuation step which goes from 0 to 1 in [0,1]. 
   The Attenuation step lies at position p, and has slope s */
double Attenuation01(double x, double s, double p)
{
  if(x<=0.0) return 0.0;
  if(x>=1.0) return 1.0;
  
  return 0.5*( 1.0 + tanh(2.0*s*p*p*(1-p)*( -1/x + ((1-p)/p)/(1-x) ) ) );
}

