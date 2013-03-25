/* ComplexFunctions.c */
/* Wolfgang Tichy 2007 */

#include "sgrid.h"
#include "Coordinates.h"

#define pow2(x)    ((x)*(x))

/* Arg function of z = x + iy : return value in (-PI,PI] */
double Arg(double x, double y)
{
  double arg;

  if(x==0.0)
  {
    if(y>0.0)  arg=PIh;
    else       arg=-PIh;
    if(y==0.0) arg=0.0;
  }
  else if(x<0)
  {
    if(y>=0.0)  arg=atan(y/x)+PI;
    else        arg=atan(y/x)-PI;
  }
  else
    arg=atan(y/x);

  return arg;
}

/* Arg function of z = x + iy : return value in [0,2PI) */
double Arg_plus(double x, double y)
{
  double arg;

  if(x==0.0)
  {
    if(y>0.0)  arg=PIh;
    else       arg=-PIh;
    if(y==0.0) arg=0.0;
  }
  else if(x<0)
  {
    if(y>=0.0)  arg=atan(y/x)+PI;
    else        arg=atan(y/x)-PI;
  }
  else
    arg=atan(y/x);

  if(arg<0.0) arg+=2.0*PI;
  return arg;
}

/* 1st derivs of Arg(x,y) and also Arg_plus(x,y) */
double dArgdx(double x, double y)
{
  return -y/(x*x + y*y);
}
double dArgdy(double x, double y)
{
  return x/(x*x + y*y);
}

/* 2nd derivs of Arg(x,y) and also Arg_plus(x,y) */
double ddArgdxdx(double x, double y)
{
  double r2 = x*x + y*y;
  return 2.0*x*y/(r2*r2);
}
double ddArgdxdy(double x, double y)
{
  double r2 = x*x + y*y;
  return 2.0*y*y/(r2*r2) - 1.0/r2;
}
double ddArgdydy(double x, double y)
{
  double r2 = x*x + y*y;
  return -2.0*x*y/(r2*r2);
}


/* Re, Im, Arg, and Abs of tanh(z)=tanh(x+iy) */
double Retanh(double x, double y)
{
  return sinh(2.0*x)/(cos(2.0*y) + cosh(2.0*x));
}
double Imtanh(double x, double y)
{
  return sin(2.0*y)/(cos(2.0*y) + cosh(2.0*x));
}
double Argtanh(double x, double y)
{
  /* return atan( sin(2.0*y)/sinh(2.0*x) ); */
  return Arg(sinh(2.0*x), sin(2.0*y));
}
double Abstanh(double x, double y)
{
  double coshx2 = cosh(2.0*x);
  double cosy2  = cos(2.0*y);
  /* return sqrt( (cosh(2.0*x)-cos(2.0*y))/(cos(2.0*y) + cosh(2.0*x)) ); */
  return sqrt( (coshx2-cosy2)/(coshx2+cosy2) );
}

/* D[Tanh[z],z] = 1/Cosh[z]^2 =: Sech[z]^2 */
/* Arg, and Abs of sech(z)=sech(x+iy) */
double Argsech(double x, double y)
{
  /* return atan( -sin(y)*sinh(x)/(cos(y)*cosh(x)) ); */
  return Arg(cos(y)*cosh(x), -sin(y)*sinh(x));
}
double Abssech(double x, double y)
{
  return sqrt( 2.0/(cos(2.0*y) + cosh(2.0*x)) );
}

/* Some derivs:
Abstanh[x_,y_] = Sqrt[ (Cosh[2 x]-Cos[2 y])/(Cos[2 y] + Cosh[2 x]) ];
Argtanh[x_,y_] = ArgWT[Sinh[2 x], Sin[2 y]];

d1Abstanh[x_,y_] = FullSimplify[D[Abstanh[x,y],x]];
d2Abstanh[x_,y_] = FullSimplify[D[Abstanh[x,y],y]];
d1Argtanh[x_,y_] = FullSimplify[D[Argtanh[x,y],x]];
d2Argtanh[x_,y_] = FullSimplify[D[Argtanh[x,y],y]];

d1d1Abstanh[x_,y_] = FullSimplify[D[d1Abstanh[x,y],x]];
d1d2Abstanh[x_,y_] = FullSimplify[D[d1Abstanh[x,y],y]];
d2d2Abstanh[x_,y_] = FullSimplify[D[d2Abstanh[x,y],y]];
d1d1Argtanh[x_,y_] = FullSimplify[D[d1Argtanh[x,y],x]];
d1d2Argtanh[x_,y_] = FullSimplify[D[d1Argtanh[x,y],y]];
d2d2Argtanh[x_,y_] = FullSimplify[D[d2Argtanh[x,y],y]];

replaceArgWTderivs = { Derivative[1,0][ArgWT][x_,y_] :> dArgdx[x,y],
                       Derivative[0,1][ArgWT][x_,y_] :> dArgdy[x,y],
                       Derivative[2,0][ArgWT][x_,y_] :> ddArgdxdx[x,y],
                       Derivative[1,1][ArgWT][x_,y_] :> ddArgdxdy[x,y],
                       Derivative[0,2][ArgWT][x_,y_] :> ddArgdydy[x,y] };
d1Argtanh[x_,y_] = d1Argtanh[x,y] /. replaceArgWTderivs;
d2Argtanh[x_,y_] = d2Argtanh[x,y] /. replaceArgWTderivs;
d1d1Argtanh[x_,y_] = d1d1Argtanh[x,y] /. replaceArgWTderivs;
d1d2Argtanh[x_,y_] = d1d2Argtanh[x,y] /. replaceArgWTderivs;
d2d2Argtanh[x_,y_] = d2d2Argtanh[x,y] /. replaceArgWTderivs;
*/
/* 1st and 2nd derivs of Argtanh(x,y) */
double dArgtanhdx(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  return 2.0*cosh(x2)*dArgdx(sinh(x2),sin(y2));
}
double dArgtanhdy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  return 2.0*cos(y2)*dArgdy(sinh(x2),sin(y2));
}
double ddArgtanhdxdx(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double sinh2x = sinh(x2);
  double sin2y = sin(y2);
  return 4.0*(pow2(cosh(x2))*ddArgdxdx(sinh2x,sin2y) + 
         dArgdx(sinh2x,sin2y)*sinh2x);
}
double ddArgtanhdxdy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  return 4.0*cos(y2)*cosh(x2)*ddArgdxdy(sinh(x2),sin(y2));
}
double ddArgtanhdydy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double sinh2x = sinh(x2);
  double sin2y = sin(y2);
  return 4.0*pow2(cos(y2))*ddArgdydy(sinh2x,sin2y) - 
         4.0*dArgdy(sinh2x,sin2y)*sin2y;
}

/* 1st and 2nd derivs of Abstanh(x,y) */
double dAbstanhdx(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double cosh2x = cosh(x2);
  double cos2y = cos(y2);
  return (2.0*cos2y*sinh(x2))/
         (pow2(cos2y + cosh2x)*
         sqrt(1.0 - (2.0*cos2y)/(cos2y + cosh2x)));
}
double dAbstanhdy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double cosh2x = cosh(x2);
  double cos2y = cos(y2);
  return (2.0*cosh2x*sin(y2))/
         (pow2(cos2y + cosh2x)*
         sqrt(1.0 - (2.0*cos2y)/(cos2y + cosh2x)));
}
double ddAbstanhdxdx(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double cosh2x = cosh(x2);
  double cos2y = cos(y2);
  return (cos2y*((-3 + 2.0*cos(y2*2.0))*cosh2x + cosh(x2*3.0) - 
         4*cos2y*pow2(sinh(x2))))/
         ((cos2y - cosh2x)*pow(cos2y + cosh2x,3)*
         sqrt(1.0 - (2.0*cos2y)/(cos2y + cosh2x)));
}
double ddAbstanhdxdy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double cosh2x = cosh(x2);
  double cos2y = cos(y2);
  return (-16.0*(pow2(cos2y) - cos2y*cosh2x + pow2(cosh2x))*
         sqrt(1.0 - (2.0*cos2y)/(cos2y + cosh2x))*sin(y2)*sinh(x2))/
         pow2(cos(y2*2.0) - cosh(x2*2.0));
}
double ddAbstanhdydy(double x, double y)
{
  double x2 = 2.0*x;
  double y2 = 2.0*y;
  double cosh2x = cosh(x2);
  double cos2y = cos(y2);
  return (4.0*cosh2x*sqrt(1.0 - (2.0*cos2y)/(cos2y + cosh2x))*
         (-5.0*cos2y + cos(y2*3.0) + 
         4.0*cosh2x*(cos2y*cosh2x + pow2(sin(y2)))))/
         pow2(cos(y2*2.0) - cosh(x2*2.0));
}
