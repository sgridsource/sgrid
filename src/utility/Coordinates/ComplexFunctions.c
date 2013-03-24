/* ComplexFunctions.c */
/* Wolfgang Tichy 2007 */

#include "sgrid.h"
#include "Coordinates.h"


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
