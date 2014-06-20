/* wrappers_for_ZIB_linSolvers.c */
/* Wolfgang Tichy 6/2014 */

/* from http://elib.zib.de/pub/elib/codelib/NewtonLib/ */

#include "sgrid.h"
#include "GridIterators.h"
#include "wrappers_for_ZIB_linSolvers.h"


#define COMPILEZIBLINSOLVERS(str) errorexits("templates_%s_wrapper: "\
"to compile with ZIB_linSolvers use MyConfig with\n"\
"DFLAGS += -DZIBLINSOLVERS\n"\
"ZIBLINSOLVERSDIR = /home/wolf/Packages/ZIB_linSolvers\n"\
"SPECIALLIBS += -L$(ZIBLINSOLVERSDIR) -lZIB_linSolvers\n", (str))

#define COPY_ARRAY_INTO_VL copy_array_into_varlist
#define COPY_VL_INTO_ARRAY copy_varlist_into_array
// #define COPY_ARRAY_INTO_VL copy_array_into_varlist_LEVEL3
// #define COPY_VL_INTO_ARRAY copy_varlist_into_array_LEVEL3

#define MAX_NGLOBALS 1024

/* global var arrays in this file */
int iglobal_forZIB=-1; /* this is used as an index to pick the set of globals in the array */
void (*lop_forZIB[MAX_NGLOBALS])(tVarList *, tVarList *, tVarList *, tVarList *);
void (*precon_forZIB[MAX_NGLOBALS])(tVarList *, tVarList *, tVarList *, tVarList *);
tVarList *r_forZIB[MAX_NGLOBALS];
tVarList *x_forZIB[MAX_NGLOBALS];
tVarList *c1_forZIB[MAX_NGLOBALS];
tVarList *c2_forZIB[MAX_NGLOBALS];
long int dim_forZIB[MAX_NGLOBALS];
tSparseVector **Acol_forZIB[MAX_NGLOBALS];

/* extern globals */
extern double *DiagMinv_JacobiPrecon; /* from wrappers_for_JacobiPrecon.c */


/* print global vars in this file */
void print_globals_forZIB(void)
{
  printf("iglobal_forZIB = %d\n", iglobal_forZIB);
  printf("lop_forZIB[iglobal_forZIB] = %p\n", lop_forZIB[iglobal_forZIB]);
  printf("precon_forZIB[iglobal_forZIB] = %p\n", precon_forZIB[iglobal_forZIB]);
  printf("r_forZIB[iglobal_forZIB] = %p -> norm %g\n", r_forZIB[iglobal_forZIB], GridL2Norm(r_forZIB[iglobal_forZIB]));
  printf("x_forZIB[iglobal_forZIB] = %p -> norm %g\n", x_forZIB[iglobal_forZIB], GridL2Norm(x_forZIB[iglobal_forZIB]));
  printf("c1_forZIB[iglobal_forZIB] = %p -> norm %g\n", c1_forZIB[iglobal_forZIB], GridL2Norm(c1_forZIB[iglobal_forZIB]));
  printf("c2_forZIB[iglobal_forZIB] = %p -> norm %g\n", c2_forZIB[iglobal_forZIB], GridL2Norm(c2_forZIB[iglobal_forZIB]));
  printf("dim_forZIB[iglobal_forZIB] = %ld\n", dim_forZIB[iglobal_forZIB]);
}


/* y := A*x , where A is matrix defined by lop */
void ZIBmatvec(int n, double *x, double *y)
{
  /* compute r = A*x */
  COPY_ARRAY_INTO_VL(x, x_forZIB[iglobal_forZIB]);
  lop_forZIB[iglobal_forZIB](r_forZIB[iglobal_forZIB],
                             x_forZIB[iglobal_forZIB],
                             c1_forZIB[iglobal_forZIB],
                             c2_forZIB[iglobal_forZIB]);
  /* copy r into y */
  COPY_VL_INTO_ARRAY(r_forZIB[iglobal_forZIB], y);
}

/* ZIBmatvectrans: */
/* y := A'*x, where A' is transpose of A which is in 
   Acol_forZIB[iglobal_forZIB] */
void ZIBAcol_times_vec_trans(int n, double *x, double *y)
{
  int ncols = (int) dim_forZIB[iglobal_forZIB];

  /* compute y = A'*x */
  SparseMatrixLines_times_vector(Acol_forZIB[iglobal_forZIB], ncols, x, y);
}


/* Precon: solves M*x = b for x */
void ZIBpsolveLEFT(int n, double *b, double *x)
{
  /* speed up special Jacobi precon */
  if(precon_forZIB[iglobal_forZIB]==Jacobi_Preconditioner_from_DiagM)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_forZIB[iglobal_forZIB]; i++)
      x[i] = b[i]*DiagMinv_JacobiPrecon[i];
  }
  else if(precon_forZIB[iglobal_forZIB]==Preconditioner_I)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_forZIB[iglobal_forZIB]; i++)
      x[i] = b[i];
  }
  else /* generic case */
  {
    /* solve M*x = b for x, solution is x_forZIB[iglobal_forZIB],
       b is in r_forZIB[iglobal_forZIB] */
    COPY_ARRAY_INTO_VL(b, r_forZIB[iglobal_forZIB]);
    precon_forZIB[iglobal_forZIB](x_forZIB[iglobal_forZIB],
//    Preconditioner_I(x_forZIB[iglobal_forZIB],
                                  r_forZIB[iglobal_forZIB],
                                  c1_forZIB[iglobal_forZIB],
                                  c2_forZIB[iglobal_forZIB]);

    /* copy x_forZIB[iglobal_forZIB] into x */
    COPY_VL_INTO_ARRAY(x_forZIB[iglobal_forZIB], x);
  }
}

/* Precon: solves M'*x = b for x.  M' is transpose of procon matrix M */
void ZIBpsolveLEFTtrans(int n, double *b, double *x)
{
  /* speed up special Jacobi precon */
  if(precon_forZIB[iglobal_forZIB]==Jacobi_Preconditioner_from_DiagM)
  {
    int i;

    SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_forZIB[iglobal_forZIB]; i++)
      x[i] = b[i]*DiagMinv_JacobiPrecon[i];
  }
  else if(precon_forZIB[iglobal_forZIB]==Preconditioner_I)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_forZIB[iglobal_forZIB]; i++)
      x[i] = b[i];
  }
  else /* generic case */
    errorexit("psolveLEFTtrans works only for Jacobi_Preconditioner_from_DiagM");
}

/* Right Precon: solves M2*x = b for x */
void ZIBpsolveRIGHT(int n, double *b, double *x)
{
  /* psolveRIGHT is not implemented yet, so use  x = b */
  int dim = dim_forZIB[iglobal_forZIB];
  int i;
  for(i=0; i<dim; i++)  x[i] = b[i];
}

/* Right Precon: solves M2'*x = b for x.
   M2' is transpose of right procon matrix M2 */
void ZIBpsolveRIGHTtrans(int n, double *b, double *x)
{
  /* psolveRIGHTtrans is not implemented yet, so use  x = b */
  ZIBpsolveRIGHT(n, b, x);
}


/* call GMRES from ZIB_linSolvers */
int ZIB_gmres_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int leftprecon = Getv("GridIterators_GMRES_PreconSide", "left");
  int rightprecon = Getv("GridIterators_GMRES_PreconSide", "right");
  int N; /* dim of matrix */
  double *B;
  double *X;
  struct ITLIN_OPT *opt;
  struct ITLIN_INFO *info;
  int RESTRT, INFO, ITER;
  double RESID;
  double norm_b = GridL2Norm(b);
  int i;

  /* set N and RESTRT */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  if(Getv("GridIterators_GMRES_restart", "max"))  RESTRT = N;
  else	RESTRT = Geti("GridIterators_GMRES_restart");
  if(RESTRT>N) RESTRT = N;

  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE ZIBgmres: ");
  if(pr) printf("  ZIB_gmres_wrapper: itmax=%d tol=%.3e\n"
                "                     N=%d RESID=%.3e RESTRT=%d\n",
                itmax, tol, N, RESID, RESTRT);
  
  /* temporary storage */
  opt = malloc(sizeof(struct ITLIN_OPT));
  info= malloc(sizeof(struct ITLIN_INFO));
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  if(opt==NULL || info==NULL || B==NULL || X==NULL)
    errorexit("ZIB_gmres_wrapper: out of memory for X, B, opt, info");

  /* increase global var index to store globals in new place in array */
  iglobal_forZIB++;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);
  if(iglobal_forZIB>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in ZIBmatvec and psolveLEFT */
  lop_forZIB[iglobal_forZIB]	= lop;
  precon_forZIB[iglobal_forZIB]	= precon;
  r_forZIB[iglobal_forZIB]	= r;
  x_forZIB[iglobal_forZIB]	= x;
  c1_forZIB[iglobal_forZIB]	= c1;
  c2_forZIB[iglobal_forZIB]	= c2;
  dim_forZIB[iglobal_forZIB]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

  /* initialize solver options */
  opt->tol = RESID;
  opt->maxiter = itmax;
  opt->i_max = RESTRT;
  opt->termcheck = CheckEachIter;

  opt->errorlevel = Minimum;
  opt->monitorlevel = Minimum; // None;
  opt->datalevel = None;
  opt->errorfile = stdout;
  opt->monitorfile = stdout;
  opt->datafile = NULL;
  opt->iterfile = NULL;
  opt->resfile  = NULL;
  opt->miscfile = NULL;
  if(pr)
  {
    opt->errorlevel = Verbose;
    opt->monitorlevel = Verbose;
  }

#ifdef ZIBLINSOLVERS
  /* call ZIBgmres solver, and decide on which side we use precon */
  /* no precon is: ZIBgmres(N, X, &ZIBmatvec, NULL,NULL, B, opt,info); */
  if(leftprecon && rightprecon)
    ZIBgmres(N, X, &ZIBmatvec, &ZIBpsolveLEFT,&ZIBpsolveLEFT, B, opt,info);
  else if(leftprecon)
    ZIBgmres(N, X, &ZIBmatvec, NULL,&ZIBpsolveLEFT, B, opt,info);
  else if(rightprecon)
    ZIBgmres(N, X, &ZIBmatvec, &ZIBpsolveLEFT,NULL, B, opt,info);
  /* ZIBgmres computes P*x (where P is the precon). So to get x we need to
     apply P^{-1} to x */
  if(leftprecon)  ZIBpsolveLEFT(N, X, X); /* this works because X can be overwritten in place */
  if(rightprecon) ZIBpsolveLEFT(N, X, X); /* this works because X can be overwritten in place */
#else
  COMPILEZIBLINSOLVERS("ZIBgmres");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_forZIB--;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);

  /* get info from solver */
  INFO=info->rcode;
  ITER=info->iter;
  RESID=info->precision;

  /* read out x and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(info);
  free(opt);

  if(pr) printf("  ZIB_gmres_wrapper: ITER=%d RESID=%.3e INFO=%d\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER ZIBgmres: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  
  /* success! */
  return ITER;
}


/* call GBIT from ZIB_linSolvers */
int ZIB_gbit_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int N; /* dim of matrix */
  double *B;
  double *X;
  struct ITLIN_OPT *opt;
  struct ITLIN_INFO *info;
  int RESTRT, INFO, ITER;
  double RESID;
  double norm_b = GridL2Norm(b);
  int i;

  /* set N and RESTRT */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  if(Getv("GridIterators_GMRES_restart", "max"))  RESTRT = N;
  else	RESTRT = Geti("GridIterators_GMRES_restart");
  if(RESTRT>N) RESTRT = N;

  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE ZIBgbit: ");
  if(pr) printf("  ZIB_gbit_wrapper: itmax=%d tol=%.3e\n"
                "                     N=%d RESID=%.3e RESTRT=%d\n",
                itmax, tol, N, RESID, RESTRT);
  
  /* temporary storage */
  opt = malloc(sizeof(struct ITLIN_OPT));
  info= malloc(sizeof(struct ITLIN_INFO));
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  if(opt==NULL || info==NULL || B==NULL || X==NULL)
    errorexit("ZIB_gbit_wrapper: out of memory for X, B, opt, info");

  /* increase global var index to store globals in new place in array */
  iglobal_forZIB++;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);
  if(iglobal_forZIB>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in ZIBmatvec and psolveLEFT */
  lop_forZIB[iglobal_forZIB]	= lop;
  precon_forZIB[iglobal_forZIB]	= precon;
  r_forZIB[iglobal_forZIB]	= r;
  x_forZIB[iglobal_forZIB]	= x;
  c1_forZIB[iglobal_forZIB]	= c1;
  c2_forZIB[iglobal_forZIB]	= c2;
  dim_forZIB[iglobal_forZIB]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

  /* initialize solver options */
  opt->tol = RESID;
  opt->rho = 0.25;
  opt->maxiter = itmax;
  opt->i_max = RESTRT;
  opt->scale = NULL;
  opt->rescale = True;

  opt->errorlevel = Minimum;
  opt->monitorlevel = Minimum; // None;
  opt->datalevel = None;
  opt->errorfile = stdout;
  opt->monitorfile = stdout;
  opt->datafile = NULL;
  opt->iterfile = NULL;
  opt->resfile  = NULL;
  opt->miscfile = NULL;
  if(pr)
  {
    opt->errorlevel = Verbose;
    opt->monitorlevel = Verbose;
  }

#ifdef ZIBLINSOLVERS
  /* call ZIBgbit solver, and decide on which side we use precon */
  ZIBgbit(N, X, &ZIBmatvec, &ZIBpsolveLEFT, B, opt,info);
#else
  COMPILEZIBLINSOLVERS("ZIBgbit");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_forZIB--;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);

  /* get info from solver */
  INFO=info->rcode;
  ITER=info->iter;
  RESID=info->precision;

  /* read out x and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(info);
  free(opt);

  if(pr) printf("  ZIB_gbit_wrapper: ITER=%d RESID=%.3e INFO=%d\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER ZIBgbit: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  
  /* success! */
  return ITER;
}

/* call PCG from ZIB_linSolvers */
int ZIB_pcg_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int N; /* dim of matrix */
  double *B;
  double *X;
  struct ITLIN_OPT *opt;
  struct ITLIN_INFO *info;
  int RESTRT, INFO, ITER;
  double RESID;
  double norm_b = GridL2Norm(b);
  int i;

  /* set N and RESTRT */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  if(Getv("GridIterators_GMRES_restart", "max"))  RESTRT = N;
  else	RESTRT = Geti("GridIterators_GMRES_restart");
  if(RESTRT>N) RESTRT = N;

  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE ZIBpcg: ");
  if(pr) printf("  ZIB_pcg_wrapper: itmax=%d tol=%.3e\n"
                "                     N=%d RESID=%.3e RESTRT=%d\n",
                itmax, tol, N, RESID, RESTRT);
  
  /* temporary storage */
  opt = malloc(sizeof(struct ITLIN_OPT));
  info= malloc(sizeof(struct ITLIN_INFO));
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  if(opt==NULL || info==NULL || B==NULL || X==NULL)
    errorexit("ZIB_pcg_wrapper: out of memory for X, B, opt, info");

  /* increase global var index to store globals in new place in array */
  iglobal_forZIB++;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);
  if(iglobal_forZIB>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in ZIBmatvec and psolveLEFT */
  lop_forZIB[iglobal_forZIB]	= lop;
  precon_forZIB[iglobal_forZIB]	= precon;
  r_forZIB[iglobal_forZIB]	= r;
  x_forZIB[iglobal_forZIB]	= x;
  c1_forZIB[iglobal_forZIB]	= c1;
  c2_forZIB[iglobal_forZIB]	= c2;
  dim_forZIB[iglobal_forZIB]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

  /* initialize solver options */
  opt->tol = RESID;
  opt->maxiter = itmax;
  opt->i_max = RESTRT;
  opt->termcheck = CheckEachIter;

  opt->errorlevel = Minimum;
  opt->monitorlevel = Minimum; // None;
  opt->datalevel = None;
  opt->errorfile = stdout;
  opt->monitorfile = stdout;
  opt->datafile = NULL;
  opt->iterfile = NULL;
  opt->resfile  = NULL;
  opt->miscfile = NULL;
  if(pr)
  {
    opt->errorlevel = Verbose;
    opt->monitorlevel = Verbose;
  }

#ifdef ZIBLINSOLVERS
  /* call ZIBpcg solver, and decide on which side we use precon */
  /* no precon is: ZIBpcg(N, X, &ZIBmatvec, NULL,NULL, B, opt,info); */
  ZIBpcg(N, X, &ZIBmatvec, &ZIBpsolveLEFT, B, opt,info);
  /* ZIBpcg computes P*x (where P is the precon). So to get x we need to
     apply P^{-1} to x */
  ZIBpsolveLEFT(N, X, X); /* this works because X can be overwritten in place */
#else
  COMPILEZIBLINSOLVERS("ZIBpcg");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_forZIB--;
  if(pr) printf("  iglobal_forZIB=%d\n", iglobal_forZIB);

  /* get info from solver */
  INFO=info->rcode;
  ITER=info->iter;
  RESID=info->precision;

  /* read out x and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(info);
  free(opt);

  if(pr) printf("  ZIB_pcg_wrapper: ITER=%d RESID=%.3e INFO=%d\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER ZIBpcg: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  
  /* success! */
  return ITER;
}
