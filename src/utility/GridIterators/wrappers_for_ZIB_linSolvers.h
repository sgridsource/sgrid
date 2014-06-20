/* wrappers_for_ZIB_linSolvers.h */
/* Wolfgang Tichy 6/2014 */

/* some stuff from ZIB_linSolvers/itlin.h */
typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL ;
typedef enum {False=0, True=1} LOGICAL ;
typedef enum { CheckOnRestart=0, CheckEachIter=1 } TERM_CHECK ;
typedef enum { Absolute=0, Relative=1 } CONV_CHECK ;

struct ITLIN_OPT
{
   double tol, rho;
   int i_max, maxiter;
   TERM_CHECK termcheck;   /* GMRES only */
   CONV_CHECK convcheck;   /* PCG only   */
   LOGICAL rescale;        /* GBIT only  */
   PRINT_LEVEL errorlevel, monitorlevel, datalevel;
   FILE *errorfile, *monitorfile, *datafile,
        *iterfile, *resfile, *miscfile;
   double *scale;
};

struct ITLIN_INFO
{
   double precision, normdx, residuum;
   int iter, rcode, subcode, nomatvec, noprecon, noprecl, noprecr;
};

struct ITLIN_DATA
{
  double *res;
  double tau, t, normdx, residuum;
  enum { GMRES=0, GBIT=1, PCG=2 } codeid;
  enum {Initial=1,Intermediate=2,Solution=3,Final=4} mode;
};

struct ITLIN_IO
{
   FILE *errfile, *monfile, *datfile,
        *iterfile, *resfile, *miscfile;
   PRINT_LEVEL errlevel, monlevel, datlevel;
};

/* for funcs below */
typedef void MATVEC(int, double*, double*);
typedef void PRECON(int, double*, double*); 

/* headers for functions called in ZIB_linSolvers */
void ZIBgmres(int n, double *y, MATVEC *matvec,
              PRECON *preconr, PRECON *preconl, double *b,
              struct ITLIN_OPT *opt, struct ITLIN_INFO *info);
void ZIBgbit(int n, double *y, MATVEC *matvec, PRECON *precon, double *b,
             struct ITLIN_OPT *opt, struct ITLIN_INFO *info);
void ZIBpcg(int n, double *y, MATVEC *matvec, PRECON *precon, double *b,
            struct ITLIN_OPT *opt, struct ITLIN_INFO *info);
