/* attenuation.c, Wolfgang Tichy 28.12.2003 */
/* useful attenuation functions */


#include <math.h>

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Tanh(x)    tanh((double) (x))
#define Tan(x)     tan((double) (x))
#define Cot(x)     (1.0/tan((double) (x)))
#define Sech(x)    (1.0/cosh((double) (x)))
#define Sec(x)     (1.0/cos((double) (x)))
#define Csc(x)     (1.0/sin((double) (x)))
#define Sqrt(x)    pow((double) (x), 0.5)
#define Abs(x)     fabs((double) (x))
#define Log(x)     log((double) (x)))

#define Pi 3.14159265358979323846264338328


/* Attenuation step which goes from 0 to 1 in [0,1]. 
   The Attenuation step lies at position p, and has slope s */
double Attenuation01(double x, double s, double p)
{
  if(x<=0.0) return 0.0;
  if(x>=1.0) return 1.0;

  return 0.5*( 1.0 + tanh(2.0*s*p*p*(1.0-p)*( -1.0/x + ((1.0-p)/p)/(1.0-x) ) ) );
}



/* my great Attenuation function */
/* good defaults are: q=1, s=2 */
double Att_0to1_x_x0_w_q_s(double x, double x0, double w, double q, double s)
{
  if(x>=x0+w)    return 1.0;
  else if(x<=x0) return 0.0;
  else return
  (1.0 + tanh((s/Pi)*( tan((x-x0)*Pi/(2.0*w)) -q*q/tan((x-x0)*Pi/(2.0*w)) )))/2.0;
}

/* first deriv of my great Attenuation function */
double dAtt_0to1_x_x0_w_q_s(double x, double x0, double w, double q, double s)
{
  if(x>=x0+w)    return 0.0;
  else if(x<=x0) return 0.0;
  else return
  (s*(Power(q,2)*Power(Csc((Pi*(x - x0))/(2.*w)),2) +
      Power(Sec((Pi*(x - x0))/(2.*w)),2))*
      Power(Sech((s*(-(Power(q,2)*Cot((Pi*(x - x0))/(2.*w))) +
                              Tan((Pi*(x - x0))/(2.*w))))/Pi),2))/(4.*w);
}

/* Inverted version of "my great Attenuation function */
double Att_1to0_x_x0_w_q_s(double x, double x0, double w, double q, double s)
{
  if(x>=x0+w)    return 0.0;
  else if(x<=x0) return 1.0;
  else return
  1.0-(1.0 + tanh((s/Pi)*( tan((x-x0)*Pi/(2.0*w)) -q*q/tan((x-x0)*Pi/(2.0*w)) )))/2.0;
}

/* Inverted version of "first deriv of my great Attenuation function */
double dAtt_1to0_x_x0_w_q_s(double x, double x0, double w, double q, double s)
{
  if(x>=x0+w)    return 0.0;
  else if(x<=x0) return 0.0;
  else return
  (-1.0)*(s*(Power(q,2)*Power(Csc((Pi*(x - x0))/(2.*w)),2) +
      Power(Sec((Pi*(x - x0))/(2.*w)),2))*
      Power(Sech((s*(-(Power(q,2)*Cot((Pi*(x - x0))/(2.*w))) +
                              Tan((Pi*(x - x0))/(2.*w))))/Pi),2))/(4.*w);
}
