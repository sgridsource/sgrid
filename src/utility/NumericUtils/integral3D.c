/* Berechne 3D Integrale mit Hilfe eines 1D Integrators */
/* Wolfgang Tichy 3.2003 */
      
/* This routine works like quad3d.c .
   It uses the 1D Integrator passed into the function in int_meth.
   It does the integral of the supplied function fn(x,y,z) over the 3D region
   defined by [x_limit1, x_limit2] , 
   [y_limit1(x), y_limit2(x)] and [z_limit1(x,y), z_limit2(x,y)]:
      
   sx = accuracy integral tries to achieve in the outer function integral3d
   sy = accuracy integral tries to achieve in the inner function f1
   sz = accuracy integral tries to achieve in the inner function f2
   max,maxy,maxz = maximum number of interations the function integral will do
               
   Normally we should have sz<sy<sx
*/


/* static vars and functions, which are global in this file */
static double x_sav,y_sav;
static double sy_sav,sz_sav;
static int maxy_sav, maxz_sav;

static double (*integrator)(double (*f_int)(double), double a, double b,
                                            double eps, int max);
static double (*nrfunc)(double,double,double);

static double (*y_lim1func)(double);
static double (*y_lim2func)(double);
static double (*z_lim1func)(double,double);
static double (*z_lim2func)(double,double);



double integral3D(double (*int_meth)(double (*f_int)(double), 
                                     double a, double b,
                                     double eps, int max),
                  double (*func)(double, double, double), 
                  double x_limit1, double x_limit2,
                  double (*y_limit1)(double x), double (*y_limit2)(double x),
                  double (*z_limit1)(double x,double y), 
                  double (*z_limit2)(double x,double y),
                  double sx, double sy, double sz, 
                  int maxx, int maxy, int maxz)
{
	static double f1(double x);
        
        sy_sav=sy;
        sz_sav=sz;
        maxy_sav=maxy;
        maxz_sav=maxz;

        integrator=int_meth;
                        
	nrfunc=func;
	
	y_lim1func=y_limit1;
	y_lim2func=y_limit2;
	z_lim1func=z_limit1;
	z_lim2func=z_limit2;
	
	return integrator(f1, x_limit1, x_limit2, sx, maxx); 
}

static double f1(double x)
{
	static double f2(double y);

	x_sav=x;
	return integrator(f2, (*y_lim1func)(x), (*y_lim2func)(x), 
	                  sy_sav, maxy_sav);
}

static double f2(double y)
{
	static double f3(double z);

	y_sav=y;
	return integrator(f3, (*z_lim1func)(x_sav,y), (*z_lim2func)(x_sav,y),
	                  sz_sav, maxz_sav);
}

static double f3(double z)
{
	return (*nrfunc)(x_sav,y_sav,z);
}
