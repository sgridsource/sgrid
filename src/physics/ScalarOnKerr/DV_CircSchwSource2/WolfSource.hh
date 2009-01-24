/* WolfSource.hh */
/* Wolfgang Tichy 1/2009 */


// structure def for constants 
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
double SourceInKerrSchild(double tKS, double xKS, double yKS, double zKS, constants * param);
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
