// WolfSource_Jan25-09c.cpp
// An error: early in SourceInKerrSchild() always ouputs the source as zero for r <=2.5M
// was corrected on 30Mar09

// WolfSource_Jan25-09b.cpp
// Generic instructions:
// To use this in a stand-alone fashion, just call standalone_main(),
// which you might want to edit.

// version: Jan25-09b included
//   #include <cstdlib>
//   #include <cstdlib>
// to make compilation easier for everyone.

// version: Jan25-09a corrected a logic error involving the initialization and
// defaults of the parameters of WolfWindow

// WolfSource_Jan25-09.cpp
// Add in a new window function suggested by Wolfgang

// RegSource_Dec15-08.cpp
// Add in function calls for psi^S
// Add in an integer parameter "smooth" that provides an unphysical
//   but smoother source function for testing convergence of numerical
//   techniques. See below.

// Add in an integer parameter "NoWindow which has the "window
//  function" set to 1 everywhere. This allows for exploring the
//  possibility of dealing with boundary conditions in different ways.

// RegSource_Oct1.cpp
//det  cd /cygdrive/c/work/code/scalarSource/RegSource
//det  cd /cygdrive/e/code/scalarSource/RegSource
//det  g++ -Wall -fno-rtti -o s  WolfSource_Jan25-09c.cpp

///////////////////////////////////////////////////////////////////////////////////
// WolfSource_Jan25-09.cpp
// computes the regularized source for a scalar charge in a circular orbit around
// Schwarzschild. This ought to be used to calculate the Scalar field self-force.
//
// Interface:
//
// Return the source-function in Kerr-Schild coordinates:
//   double SourceinKerrSchild(double tKS, double xKS, double yKS, double zKS,
//                             constants * param);
//
//
// Return the source-function in Schwarzschild coordinates:
//   double SourceInSchwarzschild(double t, double r, double theta, double phi,
//                              constants * param);
//
//
// Return Psi^S in Kerr-Schild coordinates:
//   double PsiSinKerrSchild(double tKS, double xKS, double yKS, double zKS,
//                           constants * param);
//
//
// Return Psi^S in Schwarzschild coordinates:
//   double PsiSinSchwarzschild(double t_s, double r_s, double th_s, double ph_s,
//                              constants* param);
//
//
// The original way to set parameters still works, but we now have some
// other preferred ways to set the parameters.
//   This is the (slightly modified) original way:
//   void set_parameters(double Mass, double Qcharge, double Rorbit,
//          int Nwindow, double width_window, constants * param, int smooth = 0);
//
//
// Note that we have a new parameter, an integer we call "smooth" that provides an
// unphysical, artificial source function that is guaranteed to be C^smooth.
// That is, if smooth=3, then the source will be C^3. This ought to be used to
// test numerical convergence of the numerical code when the step size is
// decreased (finite difference methods) or the number of basis functions is
// increased (spectral methods). By default, smooth=0 and no artificial smoothing
// is performed.
//
//
// Finally we have a new function that prints out the parameters in force at the time.
// The parameters go to "stdout" and each line is started with "##" which makes
// these lines ignorable by gnuplot.
//     void show_parameters(constants* p);
//
//
// Here is a caveat: We have never carefully debugged the consequences of
// setting M and Q to any values other than 1. If you want to try other
// values for M and Q it ought to work---probably! But if you are
// interested, let us know so that we can give it more attention.
//
// ******* New setup functions ********
//
// Here are completely new functions, introduced in this Jan25-09 version
//   to set parameters.
//
// New, preferred way to set parameters:  This uses default values for many of
// the parameters in the C++ style. If this is unfamiliar, let me know and I'll
// describe this in detail. I haven't used fortran in recent decades, but I
// believe that that the gnu compilers and linkers will handle this with no
// difficulty.
//
// Set the physical parameters of the system:
//   void set_orbit(constants* p, double R = 10, double M = 1, double Q = 1);
//
//   where R = radius of orbit, M = mass of point-charge, Q = charge of point-charge
//   and the default values are given in the argument list.
//
// WARNING: It is important to set the window parameters *after* the orbit is set.
//    This insures that the window ranges flank the particle orbit.
//
// Use the original window function and set its parameters:
//   void set_OrigWindow(constants* param, int Nwindow = 8, double width_window = 2,
//                       int smooth=0);
//
// Use the *New and Improved* window function and set its parameters:
//   void set_WolfWindow(constants* p,
//       double q1=1.2, double s1=1.9, double r1=0, double r2=0,
//       double q3=0, double s3=0, double r3=0, double r4=0,
//       int smooth=0);
//
// As you can see the new window function has many more parameters --- these are
//   described below. If any of the parameters are either set to 0, or are left
//   out of the list so that they pick up their default value, then the code makes
//   reasonable choices based upon the parameters given and upon the physical
//   constants R, M and Q,
//
//
// Note that the new interface sets "smooth" as a window parameter.
// By default, smooth=0 and no artificial smoothing is performed.
//
//
// Description of new WolfWindow:
// The window function has two independent parts, one part inside the
// particle, the other is outside. Between these two parts the window
// function is identically 1.
//
//
// In each of these parts there are two range parameters that determine
// where W = 0 and 1, and there are two more parameters that determine
// where W = 1/2 and what the derivative of W is at that same radius.
//
//
// There are four range parameters, r1 < r2 <= radial-orbit <= r3 < r4
//   and q1 and s1, associated with r1 to r2 and
//   and q3 and s3 associated with r3 to r4.
//   Most likely all of these parameters should be positive.
//
//
// For the range inside the particle orbit:
//   r1 is the inner radius of W where W=0 and all derivatives of W = 0,
//   and r2 is the outer radius where W=1 and all derivatives of W = 0.
//
//   q1 determines the value of r where the window function, W, is 1/2.
//   If q1=1 then the place where W=1/2 is near the middle of the range.
//
//   s1 is related to the radial derivative of W at r = q1.
//   If s1=1, then (r2-r1)*dW/dr ~ 1 at the location of W=1/2
//
//
// For the range outside the particle orbit, similar descriptions hold.
//   r3 is the inner radius where W=1 and r4 is where W=0.
//   q3 determines the value of r where the window function, W, is 1/2,
//   s3 is related to the radial derivative of W at r = q3,
//   All derivatives of W=0 at r3 and at r4.
//
//
// The window function is identically 1, with all derivatives = 0, between the
//   inner and outer ranges, that is for r2 <= r <= r3, which must include the
//   particle.
//   There appears to be no difficulty with choosing
//   r2 = r3 = location of particle.
//
//
// Currently, the defaults are chosen to make r2=r3=radius of the orbit, and
//   the inner range r1->r2 as wide as possible, but no bigger than about 20M
//   while keeping r1 a respectful distance away from the event horizon.
//
//
// The range outside the orbit defaults to make r3=orbital-radius,
//   and r4 is 20M larger than that. The effective source outside the orbit
//   appears visually very smooth and small --- presumably this is good.
//
//
// The default values for q3 and s3 are the same values as for q1 and s1.
//
//
// Here is a warning about the default range parameters: It appears that a
//   bigger range gives better results --- the evidence here is primarily
//   the desire for both Peter and Wolfgang to just set the window function
//   to zero at the start. (But this fouls boundary conditions!) But, in
//   running the code all of the grid points (I'm using finite difference
//   language) inside the window probably take a lot of time to calculate.
//   And the window is nonzero over a spherical shell --- so the number of
//   grid points inside will scale linearly with the radial extent r4-r1.
//   Now if this really makes the source smoother, then maybe it isn't
//   necessary to use so many angular grid points. This is the sort of thing
//   that should be experimented with at some point.
//
//
// If you wish to exercise the code and generate  data for a plot of Psi^S
//   or the source for one of the window functions, then just call
//             standalone_main()
//   which is the last function in the file.
//
/////////////////////////////////////////////////////////////////////////


#define SHOW(a)    " "<<#a<<" = "<<a<<" "


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
//  #include <iomanip>
//  #include <fstream>

// A reasonably arbitrary integer constant
#define CHECK 314159

struct constants
{
    int orbit_set;
    double M;
    double Q;
    double R;
    double Omega;
    double Omegad;

    int NoWindow;
    int window_set;
    int n;
    double width;
    int smooth;
    int wolf_set;
    double q1, q3;
    double s1, s3;
    double r1, r2, r3, r4;
};

// These are the prototypes for returning the source or Psi in Schw or K-S coordinates.

double SourceinKerrSchild(double tKS, double xKS, double yKS, double zKS, constants * param);
double SourceInSchwarzschild(double, double, double, double, constants *);
double PsiSinKerrSchild(double tKS, double xKS, double yKS, double zKS, constants * param);
double PsiSinSchwarzschild(double t_s, double r_s, double th_s, double ph_s, constants* param);


// These prototypes show the default vaues for the parameters.

void set_orbit(constants* p, double R = 10, double M = 1, double Q = 1);

void set_NoWindow(constants* p); // sets the window function to 1 everyhere

void set_OrigWindow(constants* p, int N = 8, double width = 2, int smooth=0);

void set_WolfWindow(constants* p,
                    double q1=1.2, double s1=1.9, double r1=0, double r2=0,
                    double q3=0, double s3=0, double r3=0, double r4=0,
                    int smooth=0);



// This function prints out the parameters in force at the time. They are
// preceded  by ## in the first two columns, which makes these lines ignorable by
// gnuplot.

void show_parameters(constants* p);



double powr(double r, int n)
{
    double temp=1.0;

    int i;
    for(i=0;i<n;i++){
        temp=temp*r;
    }

    return temp;
}

struct schw_elements
{
    double g[4][4];
    double detg;
    double Ddetg[4];
    double g_inv[4][4];
    double div_g_inv[4];
};

struct thz_elements
{
    double t,x,y,z;
    double Dx[4], Dy[4], Dz[4];
    double DDx[4][4], DDy[4][4], DDz[4][4];
};

struct source_elements
{
    double rho;
    double Psi;
    double DPsi[4];
    double DDPsi[4][4];

    double W;
    double DW[4];
    double DDW[4][4];
};

struct final_source
{
    double WPsi;
    double Del2Psi;
};


void set_orbit(constants* p, double R, double M, double Q){
    if(R<3.0*M){
        perror("Orbit too close to the black hole. Pick R>3M. \n");
        std::cout << SHOW(R) << SHOW(M) << SHOW(Q)  << std::endl;
        abort();
    }
    p->M = M;
    p->Q = Q;
    p->R = R*M;
    p->Omega = sqrt(M/(R*R*R));
    p->Omegad = p->Omega*sqrt(1.0-3.0*M/R);
    p->orbit_set = CHECK;
}

void set_OrigWindow(constants* p, int N, double width, int smooth){
    if (p->orbit_set==CHECK) {
      if (powr((p->R-2*p->M)/width, N) < 70) {
        perror("Original window extends into the black hole.\n"
               "increase ""N"" or decrease ""width"" \n");
        abort();
      }
    }
    p->n = N;
    p->width= width;
    p->wolf_set = 0;
    p->window_set = CHECK;
    p->smooth = smooth;
    p->NoWindow = 0;
}

void set_WolfWindow(constants* p,
                    double q1, double s1, double r1, double r2,
                    double q3, double s3, double r3, double r4,
                    int smooth) {
    if(p->orbit_set != CHECK) {
        perror("Please set the orbit with "
               "set_orbit(constants* p, double R = 10, double M = 1, double Q = 1)\n"
               "before setting the window parameters.");
        abort();
    }
    p->smooth = smooth;
    if (r1 == 0) r1 = p->R  - 20* p->M;
    if (r1 < 2.1* p->M ) r1 = 2.5* p->M;
    if (r2 <= r1 || r2 > p->R) r2 = p->R;

    if (r3 < p->R) r3 = p->R;
    if (r4 <= r3) r4 = r3 + 20* p->M;

    bool CheckRange = (2*(p->M) <= r1) && (r1 < r2)&& (r2 <= p->R)
                                       && (p->R <= r3 )&& (r3 < r4);
    if( !CheckRange){
        std::cout << SHOW(p->M) << SHOW(r1) << SHOW(r2) << SHOW(p->R)
                  << SHOW(r3)  << SHOW(r4) << std::endl;
        perror("Range problems. \n");
        abort();
    }

    (*p).wolf_set = CHECK;
    (*p).q1 = q1;    if (q3<=0) q3=q1; (*p).q3 = q3;
    (*p).s1 = s1;    if (s3<=0) s3=s1; (*p).s3 = s3;
    (*p).r1 = r1;    (*p).r3 = r3;
    (*p).r2 = r2;    (*p).r4 = r4;
    (*p).window_set = CHECK;
     p->NoWindow = 0;
}


void set_NoWindow(constants* p)
{ set_OrigWindow(p);
  p->NoWindow=CHECK;
}


void show_parameters(constants* p) {
  std::cout << "## " << SHOW(p->orbit_set) << SHOW(p->window_set)
            << SHOW(p->NoWindow) << SHOW(p->wolf_set) << std::endl;
  if(p->orbit_set==CHECK) std::cout << "## orbit is set: ";
  else std::cout << "## orbit is not set";
  std::cout << SHOW(p->M)
            << SHOW(p->Q)
            << SHOW(p->R)
            << SHOW(p->Omega)
            << SHOW(p->Omegad) << std::endl;
  if (p->window_set!=CHECK) {
      std::cout << "## window is not set." << std::endl;
  } else {
    if (p->wolf_set==CHECK) {
      std::cout << "## using wolfWindow: " << std::endl
                << "## "<< SHOW(p->q1) << SHOW(p->s1)
                << SHOW(p->r1) << SHOW(p->r2) << std::endl
                << "## "<< SHOW(p->q3) << SHOW(p->s3) << SHOW(p->r3) << SHOW(p->r4) << std::endl
                << "## "<< SHOW(p->smooth) << std::endl;
    } else if(p->NoWindow==CHECK) {
      std::cout <<  "## window function is set to unity everywhere."  << std::endl;
    } else {
      std::cout <<  "## using original Window: "  << std::endl
                << "## " << SHOW(p->n)
                << SHOW(p->width)
                << SHOW(p->smooth) << std::endl;
    }
  }
  std::cout << std::endl;
}


void set_parameters(double m, double q, double r, int N, double width_input,
                    constants * param, int smooth=0)
{
    if(r<6.0*m){
        perror("Orbit too close to the black hole. Pick R>6M. \n");
        abort();
    }
    (*param).orbit_set = CHECK;
    (*param).M = m;
    (*param).Q = q;
    (*param).R = r*m;
    (*param).Omega = sqrt(m/(r*r*r));
    (*param).Omegad = ((*param).Omega)*sqrt(1.0-3.0*m/r);

    (*param).n = N;
    (*param).width= width_input;
    (*param).wolf_set = 0;
    (*param).smooth = smooth;
    (*param).NoWindow = 0;
    (*param).window_set = CHECK;
}

//  set_wolf(M /*m*/, 1.0 /*q*/, 10*M /*radius of orbit*/,
//                   3.5 /*qw*/, 1.0 /*sw*/, 1. /*Wparticle*/,
//                 9.0 /*Win*/, 20 /*Wout*/, &Parameters);
void set_wolf(double m, double q, double r, double q1, double s1, double r1,  double r2,
              double q3, double s3, double r3,  double r4, constants *P)
{   if(r<=3.0*m){
        perror("Orbit too close to the black hole. Pick R>3M. \n");
        abort();
    }
    bool CheckRange = (2 <= r1) && (r1 <= r2)&& (r2 <= r)
                                       && (r <= r3 )&& (r3 <= r4);
    if( !CheckRange){
        std::cout << SHOW(2*m) << SHOW(r1) << SHOW(r2)
                  << SHOW(r) << SHOW(r3)  << SHOW(r4) << std::endl;
        perror("Range problems. \n");
        abort();
    }
    (*P).M = m;
    (*P).Q = q;
    (*P).R = r*m;
    (*P).Omega = sqrt(m/(r*r*r));
    (*P).Omegad = (*P).Omega*sqrt(1.0-3.0*m/r);
    (*P).orbit_set = CHECK;

    (*P).wolf_set = true;
    (*P).q1 = q1;      (*P).q3 = q3;
    (*P).s1 = s1;      (*P).s3 = s3;
    (*P).r1 = r1*m;    (*P).r3 = r3*m;
    (*P).r2 = r2*m;    (*P).r4 = r4*m;
    (*P).window_set = CHECK;
    (*P).NoWindow = 0;
}

void eval_schw_elements(double ts, double r, double theta, double phi,
                        constants * param, schw_elements * schw)
{
        double M=(*param).M;

        (*schw).detg=r*r*sin(theta);

        (*schw).Ddetg[0]=0.0;
        (*schw).Ddetg[1]=2.0*r*sin(theta);
        (*schw).Ddetg[2]=r*r*cos(theta);
        (*schw).Ddetg[3]=0.0;

        (*schw).div_g_inv[0]=0.0;
        (*schw).div_g_inv[1]=2.0*M/(r*r);
        (*schw).div_g_inv[2]=0.0;
        (*schw).div_g_inv[3]=0.0;

        int i;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                if(i!=j){
                    (*schw).g[i][j]=0.0;
                    (*schw).g_inv[i][j]=0.0;
                }
            }
        }

        (*schw).g[0][0]=-(1.0-2.0*M/r);
        (*schw).g[1][1]=1.0/(1.0-2.0*M/r);
        (*schw).g[2][2]=r*r;
        (*schw).g[3][3]=r*r*sin(theta)*sin(theta);

        (*schw).g_inv[0][0]=1.0/((*schw).g[0][0]);
        (*schw).g_inv[1][1]=1.0/((*schw).g[1][1]);
        (*schw).g_inv[2][2]=1.0/((*schw).g[2][2]);
        (*schw).g_inv[3][3]=1.0/((*schw).g[3][3]);
}

void eval_THZ(double ts, double r, double theta, double phi, constants * param, thz_elements * thz)
{
        static double tsLast=0, rLast=0, thetaLast=0, phiLast=0;
        if (fabs(tsLast-ts) + fabs(rLast-r) + fabs(thetaLast-theta)+ fabs(phiLast-phi) < 1.e-14){
            return;
        }
        else{
            tsLast=ts; rLast=r; thetaLast=theta; phiLast=phi;
        }

        double M=(*param).M;
        double R=(*param).R;
        double Omega=(*param).Omega;
        double Omegad=(*param).Omegad;

        double Kdelta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

        double cosTheta=cos(theta);
        double sinTheta=sin(theta);
        double cosPhi=cos(phi-Omega*ts);
        double sinPhi=sin(phi-Omega*ts);
        double cosTime=cos(Omegad*ts);
        double sinTime=sin(Omegad*ts);

        //evaluation of helping functions
        double X,Y,Z;                            //All of these are functions of (ts,r,theta, phi).
        double DX[4],DY[4],DZ[4];                //Capital letters stand for the helping functions.
        double DDX[4][4],DDY[4][4],DDZ[4][4];    //Small letters are the THZ coordinates.

        X=(r-R)/(sqrt(1.0-2.0*M/R));
        Y=R*sinTheta*sinPhi*sqrt((R-2.0*M)/(R-3.0*M));
        Z=R*cosTheta;

        DX[0]=0.0;
        DX[1]=1.0/(sqrt(1.0-2.0*M/R));
        DX[2]=0.0;
        DX[3]=0.0;

        double c=R*sqrt((R-2.0*M)/(R-3.0*M));

        DY[0]=-c*Omega*sinTheta*cosPhi;
        DY[1]=0.0;
        DY[2]=c*cosTheta*sinPhi;
        DY[3]=c*sinTheta*cosPhi;

        DZ[0]=0.0;
        DZ[1]=0.0;
        DZ[2]=-R*sinTheta;
        DZ[3]=0.0;

        int i;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                DDX[i][j]=0.0;
                DDY[i][j]=0.0;
                DDZ[i][j]=0.0;
            }
        }

        DDY[0][0]=-c*Omega*Omega*sinTheta*sinPhi;
        DDY[0][2]=-c*Omega*cosTheta*cosPhi;
        DDY[2][0]=DDY[0][2];
        DDY[0][3]=c*Omega*sinTheta*sinPhi;
        DDY[3][0]=DDY[0][3];
        DDY[2][2]=-c*sinTheta*sinPhi;
        DDY[2][3]=c*cosTheta*cosPhi;
        DDY[3][2]=DDY[2][3];
        DDY[3][3]=-c*sinTheta*sinPhi;

        DDZ[2][2]=-R*cosTheta;


        //evaluation of THZ coordinates (t,x,y,z) and their derivatives.

        (*thz).t=ts*sqrt(1.0-3.0*M/R)-(r*Omega*Y)/(sqrt(1-2.0*M/R))
            +((Omega*M*Y)/(R*R*sqrt(1-2.0*M/R)*(R-3.0*M)))
                 *(-0.5*X*X*(R-M)+M*Y*Y*(R-3.0*M)/(3.0*(R-2.0*M))+M*Z*Z)
            +(Omega*M*X*Y/(14.0*R*R*R*(R-2.0*M)*(R-3.0*M)))
            *(-X*X*(R*R-11.0*R*M+11.0*M*M)+Y*Y*(13.0*R*R-45.0*R*M+31.0*M*M)
                   +Z*Z*(13.0*R-5.0*M)*(R-2.0*M));

        double A1,B1,B2,B3,C1,D1,D2,D3;

        A1=M/(2.0*R*R*R*(R-3.0*M));
        B1=-(2.0*R-3.0*M);
        B2=(R-3.0*M);
        B3=(R-2.0*M);
        C1=M/(14.0*R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        D1=M*(13.0*R-19.0*M);
        D2=14.0*R*R-36.0*R*M+9.0*M*M;
        D3=(R-2.0*M)*(14.0*R-15.0*M);

        (*thz).z=r*cosTheta+A1*Z*(B1*X*X+B2*Y*Y+B3*Z*Z)+C1*X*Z*(D1*X*X+D2*Y*Y+D3*Z*Z);

        for(i=0;i<4;i++){

            (*thz).Dz[i]=cosTheta*Kdelta[1][i]-r*sinTheta*Kdelta[2][i]
                +A1*DZ[i]*(B1*X*X+B2*Y*Y+B3*Z*Z)
                +2.0*A1*Z*(B1*X*DX[i]+B2*Y*DY[i]+B3*Z*DZ[i])
                +C1*DX[i]*Z*(D1*X*X+D2*Y*Y+D3*Z*Z)
                +C1*X*DZ[i]*(D1*X*X+D2*Y*Y+D3*Z*Z)+2.0*C1*X*Z*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                (*thz).DDz[i][j]=-sinTheta*Kdelta[1][i]*Kdelta[2][j]
                        -sinTheta*Kdelta[1][j]*Kdelta[2][i]
                        -r*cosTheta*Kdelta[2][i]*Kdelta[2][j]
                        +A1*DDZ[i][j]*(B1*X*X+B2*Y*Y+B3*Z*Z)
                        +2.0*A1*DZ[i]*(B1*X*DX[j]+B2*Y*DY[j]+B3*Z*DZ[j])
                        +2.0*A1*DZ[j]*(B1*X*DX[i]+B2*Y*DY[i]+B3*Z*DZ[i])
                        +2.0*A1*Z*(B1*DX[i]*DX[j]+B1*X*DDX[i][j]+B2*DY[i]*DY[j]
                            +B2*Y*DDY[i][j]+B3*DZ[i]*DZ[j]+B3*Z*DDZ[i][j])
                        +C1*DDX[i][j]*Z*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +C1*DX[i]*DZ[j]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +2.0*C1*DX[i]*Z*(D1*X*DX[j]+D2*Y*DY[j]+D3*Z*DZ[j])
                        +C1*DX[j]*DZ[i]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +C1*X*DDZ[i][j]*(D1*X*X+D2*Y*Y+D3*Z*Z)
                        +2.0*C1*X*DZ[i]*(D1*X*DX[j]+D2*Y*DY[j]+D3*Z*DZ[j])
                        +2.0*C1*DX[j]*Z*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i])
                        +2.0*C1*X*DZ[j]*(D1*X*DX[i]+D2*Y*DY[i]+D3*Z*DZ[i])
                        +2.0*C1*X*Z*(D1*DX[i]*DX[j]+D1*X*DDX[i][j]+D2*DY[i]*DY[j]
                            +D2*Y*DDY[i][j]+D3*DZ[i]*DZ[j]+D3*Z*DDZ[i][j]);

            }
        }

               //This is for x_tilde and y_tilde, functions which are used to evaluate x and y.

        double xt,yt;
        double Dxt[4], Dyt[4];
        double DDxt[4][4], DDyt[4][4];

                                //y_tilde.

        double E1,E2,E3;

        A1=sqrt((R-2.0*M)/(R-3.0*M));
        B1=M/(2.0*R*R*R);
        C1=(R-3.0*M)/(R-2.0*M);
        D1=M/(14.0*R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        E1=2.0*M*(4.0*R-15.0*M);
        E2=14.0*R*R-69.0*M*R+89.0*M*M;
        E3=2.0*(R-2.0*M)*(7.0*R-24.0*M);

        yt=A1*r*sinTheta*sinPhi+B1*Y*(-2.0*X*X+C1*Y*Y+Z*Z)
            +D1*X*Y*(E1*X*X+E2*Y*Y+E3*Z*Z);

        for(i=0;i<4;i++){

            Dyt[i]=A1*sinTheta*sinPhi*Kdelta[1][i]
                  +A1*r*cosTheta*sinPhi*Kdelta[2][i]
                  +A1*r*sinTheta*cosPhi*Kdelta[3][i]
                  -A1*Omega*r*sinTheta*cosPhi*Kdelta[0][i]
                  +B1*DY[i]*(-2.0*X*X+C1*Y*Y+Z*Z)
                  +2.0*B1*Y*(-2.0*X*DX[i]+C1*Y*DY[i]+Z*DZ[i])
                  +D1*(Y*DX[i]+X*DY[i])*(E1*X*X+E2*Y*Y+E3*Z*Z)
                  +2.0*D1*X*Y*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDyt[i][j]=A1*cosTheta*sinPhi*Kdelta[1][i]*Kdelta[2][j]
                          +A1*sinTheta*cosPhi*Kdelta[1][i]*Kdelta[3][j]
                          -A1*Omega*sinTheta*cosPhi*Kdelta[1][i]*Kdelta[0][j]
                          +A1*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[1][j]
                          -A1*r*sinTheta*sinPhi*Kdelta[2][i]*Kdelta[2][j]
                          +A1*r*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[3][j]
                          -A1*Omega*r*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[0][j]
                          +A1*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[1][j]
                          +A1*r*cosTheta*cosPhi*Kdelta[3][i]*Kdelta[2][j]
                          -A1*r*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[3][j]
                          +A1*Omega*r*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[0][j]
                          -A1*Omega*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[1][j]
                          -A1*Omega*r*cosTheta*cosPhi*Kdelta[0][i]*Kdelta[2][j]
                          +A1*Omega*r*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[3][j]
                          -A1*Omega*Omega*r*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[0][j]
                          +B1*DDY[i][j]*(-2.0*X*X+C1*Y*Y+Z*Z)
                          +2.0*B1*DY[i]*(-2.0*X*DX[j]+C1*Y*DY[j]+Z*DZ[j])
                          +2.0*B1*DY[j]*(-2.0*X*DX[i]+C1*Y*DY[i]+Z*DZ[i])
                          +2.0*B1*Y*(-2.0*DX[j]*DX[i]-2.0*X*DDX[i][j]+C1*DY[i]*DY[j]+C1*Y*DDY[i][j]
                            +DZ[i]*DZ[j]+Z*DDZ[i][j])
                          +D1*DDX[i][j]*Y*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +D1*DX[i]*DY[j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*DX[i]*Y*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +D1*DX[j]*DY[i]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +D1*X*DDY[i][j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*X*DY[i]*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +2.0*D1*DX[j]*Y*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*DY[j]*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*Y*(E1*DX[i]*DX[j]+E1*X*DDX[i][j]+E2*DY[i]*DY[j]+E2*Y*DDY[i][j]
                            +E3*DZ[i]*DZ[j]+E3*Z*DDZ[i][j]);

            }
        }

                                                //x_tilde.

        double F1,G1,G2,G3,G4,G5,G6;

        A1=1.0/(sqrt(1.0-2.0*M/R));
        B1=M/(R*R*sqrt(1.0-2.0*M/R));
        C1=(R-3.0*M)/(R-2.0*M);
        D1=M/(2.0*R*R*R*(R-2.0*M)*(R-3.0*M));
        E1=-M*M;
        E2=(R-3.0*M)*(3.0*R-8.0*M);
        E3=3.0*(R-2.0*M)*(R-2.0*M);
        F1=M/(R*R*R*R*R*sqrt(1.0-2.0*M/R)*(R-3.0*M));
        G1=M*(R*R-R*M+3.0*M*M)/(8.0*(R-2.0*M));
        G2=(28.0*R*R-114.0*R*M+123.0*M*M)/28.0;
        G3=(14.0*R*R-48.0*R*M+33.0*M*M)/14.0;
        G4=M*(3.0*R*R*R-74.0*R*R*M+337.0*R*M*M-430.0*M*M*M)/(56.0*(R-2.0*M)*(R-2.0*M));
        G5=-M*M*(7.0*R-18.0*M)/(4.0*(R-2.0*M));
        G6=-M*(3.0*R+22.0*M)/56.0;

        xt=A1*(r*sinTheta*cosPhi-R)+B1*(-0.5*X*X+C1*Y*Y+Z*Z)
            +D1*X*(E1*X*X+E2*Y*Y+E3*Z*Z)
            +F1*(G1*X*X*X*X+G2*X*X*Y*Y+G3*X*X*Z*Z+G4*Y*Y*Y*Y+G5*Y*Y*Z*Z+G6*Z*Z*Z*Z);

        for(i=0;i<4;i++){

            Dxt[i]=A1*sinTheta*cosPhi*Kdelta[1][i]
                  +A1*r*cosTheta*cosPhi*Kdelta[2][i]
                  -A1*r*sinTheta*sinPhi*Kdelta[3][i]
                  +A1*Omega*r*sinTheta*sinPhi*Kdelta[0][i]
                  +B1*(-X*DX[i]+2.0*C1*Y*DY[i]+2.0*Z*DZ[i])
                  +D1*DX[i]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                  +2.0*D1*X*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                  +F1*(4.0*G1*X*X*X*DX[i]+2.0*G2*X*DX[i]*Y*Y+2.0*G2*X*X*Y*DY[i]
                    +2.0*G3*X*DX[i]*Z*Z+2.0*G3*X*X*Z*DZ[i]+4.0*G4*Y*Y*Y*DY[i]
                    +2.0*G5*Y*DY[i]*Z*Z+2.0*G5*Y*Y*Z*DZ[i]+4.0*G6*Z*Z*Z*DZ[i]);

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                DDxt[i][j]=A1*cosTheta*cosPhi*Kdelta[1][i]*Kdelta[2][j]
                          -A1*sinTheta*sinPhi*Kdelta[1][i]*Kdelta[3][j]
                          +A1*Omega*sinTheta*sinPhi*Kdelta[1][i]*Kdelta[0][j]
                          +A1*cosTheta*cosPhi*Kdelta[2][i]*Kdelta[1][j]
                          -A1*r*sinTheta*cosPhi*Kdelta[2][i]*Kdelta[2][j]
                          -A1*r*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[3][j]
                          +A1*Omega*r*cosTheta*sinPhi*Kdelta[2][i]*Kdelta[0][j]
                          -A1*sinTheta*sinPhi*Kdelta[3][i]*Kdelta[1][j]
                          -A1*r*cosTheta*sinPhi*Kdelta[3][i]*Kdelta[2][j]
                          -A1*r*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[3][j]
                          +A1*Omega*r*sinTheta*cosPhi*Kdelta[3][i]*Kdelta[0][j]
                          +A1*Omega*sinTheta*sinPhi*Kdelta[0][i]*Kdelta[1][j]
                          +A1*Omega*r*cosTheta*sinPhi*Kdelta[0][i]*Kdelta[2][j]
                          +A1*Omega*r*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[3][j]
                          -A1*Omega*Omega*r*sinTheta*cosPhi*Kdelta[0][i]*Kdelta[0][j]
                          +B1*(-DX[i]*DX[j]-X*DDX[i][j]+2.0*C1*DY[j]*DY[i]+2.0*C1*Y*DDY[i][j]
                               +2.0*DZ[i]*DZ[j]+2.0*Z*DDZ[i][j])
                          +D1*DDX[i][j]*(E1*X*X+E2*Y*Y+E3*Z*Z)
                          +2.0*D1*DX[i]*(E1*X*DX[j]+E2*Y*DY[j]+E3*Z*DZ[j])
                          +2.0*D1*DX[j]*(E1*X*DX[i]+E2*Y*DY[i]+E3*Z*DZ[i])
                          +2.0*D1*X*(E1*DX[i]*DX[j]+E1*X*DDX[i][j]+E2*DY[i]*DY[j]+E2*Y*DDY[i][j]
                                     +E3*DZ[i]*DZ[j]+E3*Z*DDZ[i][j])
                          +F1*(12.0*G1*X*X*DX[j]*DX[i]+4.0*G1*X*X*X*DDX[i][j]+2.0*G2*DX[j]*DX[i]*Y*Y
                              +2.0*G2*X*DDX[i][j]*Y*Y+4.0*G2*X*DX[i]*Y*DY[j]+4.0*G2*X*DX[j]*Y*DY[i]
                              +2.0*G2*X*X*DY[j]*DY[i]+2.0*G2*X*X*Y*DDY[i][j]+2.0*G3*DX[j]*DX[i]*Z*Z
                              +2.0*G3*X*DDX[i][j]*Z*Z+4.0*G3*X*DX[i]*Z*DZ[j]+4.0*G3*X*DX[j]*Z*DZ[i]
                              +2.0*G3*X*X*DZ[j]*DZ[i]+2.0*G3*X*X*Z*DDZ[i][j]+12.0*G4*Y*Y*DY[i]*DY[j]
                              +4.0*G4*Y*Y*Y*DDY[i][j]+2.0*G5*DY[j]*DY[i]*Z*Z+2.0*G5*Y*DDY[i][j]*Z*Z
                              +4.0*G5*Y*DY[i]*Z*DZ[j]+4.0*G5*Y*DY[j]*Z*DZ[i]+2.0*G5*Y*Y*DZ[j]*DZ[i]
                              +2.0*G5*Y*Y*Z*DDZ[i][j]+12.0*G6*Z*Z*DZ[j]*DZ[i]+4.0*G6*Z*Z*Z*DDZ[i][j]);
            }
        }

                                 //x(ts,r,theta,phi) and derivatives.

        (*thz).x=xt*cosTime-yt*sinTime;

        for(i=0;i<4;i++){

            (*thz).Dx[i]=Dxt[i]*cosTime-xt*Omegad*sinTime*Kdelta[0][i]
                -Dyt[i]*sinTime-yt*Omegad*cosTime*Kdelta[0][i];

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

            (*thz).DDx[i][j]=DDxt[i][j]*cosTime
                            -(Dxt[i]*Kdelta[0][j]+Dxt[j]*Kdelta[0][i])*Omegad*sinTime
                            -xt*Omegad*Omegad*cosTime*Kdelta[0][i]*Kdelta[0][j]
                            -DDyt[i][j]*sinTime
                            -(Dyt[i]*Kdelta[0][j]+Dyt[j]*Kdelta[0][i])*Omegad*cosTime
                            +yt*Omegad*Omegad*sinTime*Kdelta[0][i]*Kdelta[0][j];
            }
        }
                                            //y(ts,r,theta,phi) and derivatives.

        (*thz).y=xt*sinTime+yt*cosTime;

        for(i=0;i<4;i++){

            (*thz).Dy[i]=Dxt[i]*sinTime+xt*Omegad*cosTime*Kdelta[0][i]
                +Dyt[i]*cosTime-yt*Omegad*sinTime*Kdelta[0][i];

        }

        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){

                (*thz).DDy[i][j]=DDxt[i][j]*sinTime
                        +(Dxt[i]*Kdelta[0][j]+Dxt[j]*Kdelta[0][i])*Omegad*cosTime
                        -xt*Omegad*Omegad*sinTime*Kdelta[0][i]*Kdelta[0][j]
                        +DDyt[i][j]*cosTime
                        -(Dyt[i]*Kdelta[0][j]+Dyt[j]*Kdelta[0][i])*Omegad*sinTime
                        -yt*Omegad*Omegad*cosTime*Kdelta[0][i]*Kdelta[0][j];
            }
        }
}

// Schwarzschild coordinates
void eval_source_elements(double ts,double r, double theta, double phi,
                          constants * param, thz_elements * thz, source_elements * sourceW)
{
    //double M=(*param).M;
    double R=(*param).R;

    double Kdelta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    //rho, Psi and derivatives.

    double Q=(*param).Q;
    double x=(*thz).x;
    double y=(*thz).y;
    double z=(*thz).z;

    (*sourceW).rho=sqrt(x*x+y*y+z*z);
    (*sourceW).Psi=Q/((*sourceW).rho);

    int i;
    for(i=0;i<4;i++){
        (*sourceW).DPsi[i]=-(Q/(((*sourceW).rho)*((*sourceW).rho)*((*sourceW).rho)))
          * ( x*((*thz).Dx[i])+y*((*thz).Dy[i])+z*((*thz).Dz[i]) );
    }

    double rho = (*sourceW).rho;

    for(i=0;i<4;i++){
        int j;
        for(j=0;j<4;j++){
        (*sourceW).DDPsi[i][j]=(3.0*Q/(rho*rho*rho*rho*rho))*(x*((*thz).Dx[i])+y*((*thz).Dy[i])+z*((*thz).Dz[i]))
            *(x*((*thz).Dx[j])+y*((*thz).Dy[j])+z*((*thz).Dz[j]))-(Q/(rho*rho*rho))*(x*((*thz).DDx[i][j])+((*thz).Dx[i])*((*thz).Dx[j])+
            y*((*thz).DDy[i][j])+((*thz).Dy[i])*((*thz).Dy[j])+z*((*thz).DDz[i][j])+((*thz).Dz[i])*((*thz).Dz[j]));
        }
    }

    double dWdr=0;
    double dWdrdr=0;

    //window stuff
    if (param->NoWindow == CHECK) {
        (*sourceW).W = 1;
        for(i=0;i<4;i++){
            sourceW->DW[i]= 0;
            int j;
            for(j=0;j<4;j++){
                sourceW->DDW[i][j]= 0;
            }
        }
    } else if (! ((*param).wolf_set==CHECK) ) { // original window
        double W;
        double width=(*param).width;
        //double width = 1.0*M;
        double rl=(r-R)/width;
        int n = (*param).n;
        (*sourceW).W = W = exp(-powr(rl,n));
        for(i=0;i<4;i++){
            (*sourceW).DW[i]=-W*(n*powr(rl,n-1)*(Kdelta[1][i])/width) ;
        }
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                (*sourceW).DDW[i][j]=-((*sourceW).DW[j])*(n*powr(rl,n-1)*(Kdelta[1][i])/width)
                    -W*(n*(n-1)*powr(rl,n-2)*(Kdelta[1][i])*(Kdelta[1][j])/(width*width));
            }
        }
        return;
    } else {
        if((*param).wolf_set != CHECK) { //  wolf Window
          perror("Logic error in eval_source_elements()");
          exit(1);
        }
        double Pi=3.14159265358979323846264338;
        if (((*param).r1<r && r <(*param).r2) || ((*param).r3<r && r <(*param).r4)) {
//      within a window range
          double q2, s, w, Tan;
          int SIGN;
          if ((*param).r1<r && r <(*param).r2) { // within inner window range
            SIGN = 1;
            q2 = (*param).q1*(*param).q1;
            s = (*param).s1;
            w = (*param).r2 -(*param).r1;
            Tan = tan(Pi*(r-(*param).r1)/2./w);
          } else {  // within outer window range
            SIGN = -1;
//det corrected 1/27/2009    q2 = (*param).q3*(*param).q1;
            q2 = (*param).q3*(*param).q3;
            s = (*param).s3;
            w = (*param).r4 -(*param).r3;
            Tan = tan(Pi*(r-(*param).r4)/2./w);
          }
          double q4 = q2*q2;
          double w2 = w*w;
          double Tan2 = Tan*Tan;
          double Tan4 = Tan2*Tan2;
          double Tanh = tanh(s/Pi*(Tan-q2/Tan));
          double Tanh2 = Tanh*Tanh;
// SIGN=+1 for the inner window range, and SIGN=-1 for the outer window range.
          (*sourceW).W = 1./2.*(1+SIGN*Tanh);
// derivatives w.r.t. Schwarzschild r
          dWdr = SIGN*1/4.*(1-Tanh2)*s*(1+Tan2)/w*(1+q2/Tan2);
          double dWdrdr_check = SIGN*s/4.*(1-Tanh2)*(1+Tan2)
                    *(-s*(1+Tan2)*(3*Tan4+2*Tan2*(q2-Tan2)+q4)*Tanh
                        + Pi*Tan*(Tan4-q2))/Tan4/w2;
          dWdrdr = SIGN*s/4.*(1-Tanh2)*(1+Tan2)/Tan4/w2
               *( -s*(1+Tan2)*(Tan4+2*Tan2*q2+q4)*Tanh + Pi*Tan*(Tan4-q2) );
          assert(fabs(dWdrdr -dWdrdr_check) <= 1.e-14);
          // end inside range
        } else if (((*param).r2 <=r) && ( r <= (*param).r3)) { // between ranges
          (*sourceW).W =1;
        } else { // outside or inside the windows
          (*sourceW).W =0;
        }
        for(i=0;i<4;i++){
          (*sourceW).DW[i]=0;
        }
        for(i=0;i<4;i++){
          int j;
          for(j=0;j<4;j++){
            (*sourceW).DDW[i][j]=0;
          }
        }
        (*sourceW).DW[1] = dWdr;
        (*sourceW).DDW[1][1]=dWdrdr;
      }
//  printf("%f  %3.10e  %3.10e  %3.10e \n", r, (*sourceW).W, dWdr, dWdrdr);
}

double eval_source(double ts, double r, double theta, double phi,
                   constants * param, schw_elements * schw, thz_elements * thz,
                   source_elements * sourceW)
{
    double sum=0.0;
    double R=(*param).R;
    double M=(*param).M;
    double Pi=3.14159265358979323846264338;

    if(1==2 /*Det fabs(r-R)>3.5*M */ ){
        return 0.0;
    }
    else{

        eval_THZ(ts,r,theta,phi,param,thz);

        double t=(*thz).t;
        double x=(*thz).x;
        double y=(*thz).y;
        double z=(*thz).z;

// Keep the field point away from troublesome points too near the charge.
        if(sqrt(x*x+y*y+z*z)<0.005*M)
        {
            r=R+0.006*M;

            eval_THZ(ts,r,theta,phi,param,thz);
            t=(*thz).t;
            x=(*thz).x;
            y=(*thz).y;
            z=(*thz).z;
        }

        if(theta==0.0)
        {
            theta+=1.0e-12;

            eval_THZ(ts,r,theta,phi,param,thz);
            t=(*thz).t;
            x=(*thz).x;
            y=(*thz).y;
            z=(*thz).z;
        }

        if(fabs(theta-Pi)<1.0e-12)
        {
            theta-=1.0e-12;

            eval_THZ(ts,r,theta,phi,param,thz);
            t=(*thz).t;
            x=(*thz).x;
            y=(*thz).y;
            z=(*thz).z;
        }

        eval_schw_elements(ts,r,theta,phi,param,schw);
        eval_source_elements(ts,r,theta,phi,param,thz,sourceW);

        double Psi=(*sourceW).Psi;
        double W=(*sourceW).W;

        int i;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                sum+=((*schw).Ddetg[i])*((*schw).g_inv[i][j])*((*sourceW).DPsi[j]);
            }
        }
        sum=sum/((*schw).detg);
        for(i=0;i<4;i++){
            sum+=((*schw).div_g_inv[i])*((*sourceW).DPsi[i]);
        }
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                sum+=((*schw).g_inv[i][j])*((*sourceW).DDPsi[i][j]);
            }
        }

  //window modifications.
        sum=sum*W;
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                sum+=(Psi/((*schw).detg))*((*schw).Ddetg[i])*((*schw).g_inv[i][j])*((*sourceW).DW[j]);
            }
        }
        for(i=0;i<4;i++){
            sum+=Psi*((*schw).div_g_inv[i])*((*sourceW).DW[i]);
        }
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                sum+=2.0*((*schw).g_inv[i][j])*((*sourceW).DPsi[i])*((*sourceW).DW[j]);
            }
        }
        for(i=0;i<4;i++){
            int j;
            for(j=0;j<4;j++){
                sum+=Psi*((*schw).g_inv[i][j])*((*sourceW).DDW[i][j]);
            }
        }

        if((*param).smooth==0) return sum;
        else {
         double rho = (*sourceW).rho/param->M;
         return sum*powr(rho/(1+rho*rho),param->smooth);
        }
    }
}

double SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS, constants * param)
{
    schw_elements schw;
    thz_elements thz;
    source_elements sourceW;

    double M=(*param).M;

    double Pi = 3.14159265358979323846264338;
    double r = sqrt(xKS*xKS+yKS*yKS+zKS*zKS);

//    if(1==2 /*Det r<=2.5*M */){
// det corrected on 30Mar09:    if( r<=2.5*M ){
    if( r<=2.0*M ){
        return 0.0;
    }
    else{
        double theta;
        if(zKS != 0.0){
            if(zKS>0.0) theta = atan(sqrt(xKS*xKS+yKS*yKS)/zKS);
            else theta = Pi-atan(sqrt(xKS*xKS+yKS*yKS)/fabs(zKS));
        }
        else theta = Pi/2.0;

        double phi;  // 0 <= phi < 2*Pi
        if(xKS == 0.0){
            if(yKS>= 0.0) phi = Pi/2.0;
            else phi = -Pi/2.0;
        }
        else if(xKS > 0.0){
            if(yKS>=0.0) phi = atan(yKS/xKS);
            else phi = 2.0*Pi-atan(fabs(yKS)/xKS);
        }
        else{
            if(yKS>=0.0) phi = Pi-atan(yKS/fabs(xKS));
            else phi = Pi+atan(fabs(yKS/xKS));
        }

        double ts = tKS-2.0*M*log(r/(2.0*M)-1.0);
// std::cerr << SHOW(xKS) << SHOW(r*sin(theta)*cos(phi)) << SHOW(fabs(xKS-r*sin(theta)*cos(phi))) <<std::endl;
        assert(fabs(xKS-r*sin(theta)*cos(phi))<1.e-12);
        assert(fabs(yKS-r*sin(theta)*sin(phi))<1.e-12);
        assert(fabs(zKS-r*cos(theta))<1.e-12);

        //eval_source(ts,r,theta,phi,param,&schw,&thz,&sourceW,source);
        //return (*source).Del2Psi;

        return eval_source(ts,r,theta,phi,param,&schw,&thz,&sourceW);
   }
}


double PsiSinKerrSchild(double tKS, double xKS, double yKS, double zKS, constants * param)
{
    schw_elements schw;
    thz_elements thz;
    source_elements sourceW;

    double M=(*param).M;

    double Pi = 3.14159265358979323846264338;
    double r = sqrt(xKS*xKS+yKS*yKS+zKS*zKS);

    if(1==2 /*Det r<=2.5*M */){
        return 0.0;
    }
    else{
        double theta;
        if(zKS != 0.0){
            if(zKS>0.0) theta = atan(sqrt(xKS*xKS+yKS*yKS)/zKS);
            else theta = Pi-atan(sqrt(xKS*xKS+yKS*yKS)/fabs(zKS));
        }
        else theta = Pi/2.0;

        double phi;
        if(xKS == 0.0){
            if(yKS>= 0.0) phi = Pi/2.0;
            else phi = -Pi/2.0;
        }
        else if(xKS > 0.0){
            if(yKS>=0.0) phi = atan(yKS/xKS);
            else phi = 2.0*Pi-atan(fabs(yKS)/xKS);
        }
        else{
            if(yKS>=0.0) phi = Pi-atan(yKS/fabs(xKS));
            else phi = Pi+atan(fabs(yKS/xKS));
        }

        double ts = tKS-2.0*M*log(r/(2.0*M)-1.0);

        assert(fabs(xKS-r*sin(theta)*cos(phi))<1.e-12);
        assert(fabs(yKS-r*sin(theta)*sin(phi))<1.e-12);
        assert(fabs(zKS-r*cos(theta))<1.e-12);
        /* return */ eval_source(ts,r,theta,phi,param,&schw,&thz,&sourceW);
        return sourceW.Psi*sourceW.W;
   }
}



//Det stuff added below
double SourceInSchwarzschild(double t_s, double r_s, double th_s, double ph_s,
                             constants* param)
{
    schw_elements schw;
    thz_elements thz;
    source_elements sourceW;
    return eval_source(t_s,r_s,th_s,ph_s,param,&schw,&thz,&sourceW);
}


double PsiSinSchwarzschild(double t_s, double r_s, double th_s, double ph_s,
                           constants* param)
{
    schw_elements schw;
    thz_elements thz;
    source_elements sourceW;
    eval_source(t_s,r_s,th_s,ph_s,param,&schw,&thz,&sourceW);
//    return sourceW.Del2Psi;
    return sourceW.Psi*sourceW.W;
}



//  void set_orbit(constants* p, double R = 10, double M = 1, double Q = 1);
//
//  void set_OrigWindow(constants* p, int N = 8, double width = 2, int smooth=0);
//
//  void set_WolfWindow(constants* p,
//                      double q1=1.2, double s1=1.9, double r1=0, double r2=0,
//                      double q3=0, double s3=0, double r3=0, double r4=0,
//                      int smooth=0);
//
//  void show_parameters(constants* p);

//  void set_OrigWindow(constants* p, int N = 8, double width = 2, int smooth=0);
//
//  void set_WolfWindow(constants* p,
//                      double q1=1.2, double s1=1.9, double r1=0, double r2=0,
//                      double q3=0, double s3=0, double r3=0, double r4=0,
//                      int smooth=0);

// Set_NoWindow(constants* p);

int Det_main(int argc, char * argv[])
{
  constants PP;
  set_orbit(&PP);
  set_OrigWindow(&PP);
  show_parameters(&PP);

  set_WolfWindow(&PP);
  show_parameters(&PP);

exit(0);
  using namespace std;

// Read the command line
  int argn = 1;
// Command line arguments: q1, s1, r1, r2, q3, s3, r3, r4
  double Q1 = 1.2;
  double S1 = 1.9;
  double R1 = 2.0;
  double R2 = 10;
  double Q3 = 1.2;
  double S3 = 1.9;
  double R3 = 10.;
  double R4 = 20.0;
  if (argc > argn) {
    Q1 = atof(argv[argn]);
//    if ( Q1 <= 6) {
//      cerr << endl << " Radius " << RADIUS << " is out of bounds" << std::endl;
//      exit(1);
//    }
  }
  argn++;
  std::cout << "## " << Q1 << " ";

  if (argc > argn) {
    S1 = atof(argv[argn]);
  }
  argn++;
  std::cout << S1 << " ";

  if (argc > argn) {
    R1 = atof(argv[argn]);
    if ((R1 < 2) || (R1 > 20)) {
      std::cerr << endl << " R1 " << R1 << " is out of bounds" << std::endl;
      exit(1);
    }
  }
  argn++;
  std::cout << R1 << " ";

  if (argc > argn) {
    R2 = atof(argv[argn]);
    if ((R2 < 2) || (R2 > 20) || (R2 < R1)) {
      std::cerr << endl << " R2 " << R2 << "inner range is bad" << std::endl;
      exit(1);
    }
  }
  argn++;
  std::cout << R2 << " " << std::endl;

  if (argc > argn) {
    Q3 = atof(argv[argn]);
  }
  argn++;
  std::cout << "## " << Q3 << " ";

  if (argc > argn) {
    S3 = atof(argv[argn]);
  }
  argn++;
  std::cout << S3 << " ";

  if (argc > argn) {
    R3 = atof(argv[argn]);
    if ((R3 < 2) || (R3 > 20)) {
      std::cerr << endl << " R3 " << R3 << " is out of bounds" << std::endl;
      exit(1);
    }
  }
  argn++;
  std::cout << R3 << " ";

  if (argc > argn) {
    R4 = atof(argv[argn]);
    if ((R4 < 2) || (R4 < R3)) {
      std::cerr << endl << " R4 " << R4 << "outer range is bad" << std::endl;
      exit(1);
    }
  }
  argn++;
  std::cout << R4 << " " << std::endl;

// End of reading the command line

  double x,y,z;
  double M=1.0;

//  int n=8;             //exponent in window function.
//  double width=2.0*M; //width of window function.
  constants P;       //This structure contains all the initializing constants.


  set_orbit(&P, 10.0);  // void set_orbit(constants* p, double R = 10, double M = 1, double Q = 1);

  set_WolfWindow(&P, Q1, S1, R1, R2, Q3, S3, R3, R4);
  show_parameters(&P);

    for(x=0; x <20.1; x+=0.5) {
      for(y=-10; y <10.1; y+=0.5) {
        z = 0; // r*cos(th);
        if ((x*x+y*y+z*z) < 4.1*M*M) { // Inside the hole
          printf("%f %f 0.0 \n", x, y);
        } else {
          printf("%f %f %3.10e \n",
                 x, y,  SourceInKerrSchild(4.0*log(2.0)*M,x,y,z,&P));
        }
    }
    printf("\n");
  }
  return 0;
}

int standalone_main() {
// int main() {
    constants P;

     set_orbit(&P, 10);
     set_WolfWindow(&P, 0.6, 3.6, 2.1, 0, 1.2, 1.9);

//  > q1     = 0.6
//  > s1     = 3.6
//  > r1     = 2.1
//  > q3     = 1.2
//  > s3     = 1.9);
     show_parameters(&P);

    double x,y,z;
    for(x=0; x < 12. ; x+=0.3) {
      for(y=-10.0; y <10.1; y+=.3) {
        z = 0; // r*cos(th);
        if ((x*x+y*y+z*z) <= 4.0*(P.M)*(P.M)) { // Inside the hole
          printf("%f %f 0.0 \n", x, y);
        } else {

          printf("%f %f %3.10e \n",
                  x, y, SourceInKerrSchild(4.0*log(2.0)*(P.M),x,y,z,&P));

//          double PsiS = PsiSinKerrSchild(4.0*log(2.0)*P.M,x,y,z,&P);
//          if (fabs(PsiS) > 4) PsiS = 0;
//          printf("%f %f %3.10e \n", x, y, PsiS);

// //        printf("%f %f %3.10e \n",
// //               r, ph, SourceInSchwarzschild(0.0*M,r,th,ph*3.14159,&P));
      }
    }
    printf("\n");
  }
  return 0;
}
