/* wrappers_for_templates.c */
/* Wolfgang Tichy 8/2007 */

/* from "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates/				*/

#include "sgrid.h"
#include "GridIterators.h"
#include "wrappers_for_templates.h"


#define COMPILETEMPLATES(str) errorexits("templates_%s_wrapper: "\
"to compile with templates use MyConfig with\n"\
"DFLAGS += -DTEMPLATES\n"\
"TEMPLATESDIR = /home/wolf/Packages/dctemplates_extBlasLapack\n"\
"SPECIALLIBS += -L$(TEMPLATESDIR) -L$(TEMPLATESDIR)/F2CLIBS \\\n"\
"-literatortemplates -lblas -llapack # -lI77 -lF77\n", (str))

#define COPY_ARRAY_INTO_VL copy_array_into_varlist
#define COPY_VL_INTO_ARRAY copy_varlist_into_array
// #define COPY_ARRAY_INTO_VL copy_array_into_varlist_LEVEL3
// #define COPY_VL_INTO_ARRAY copy_varlist_into_array_LEVEL3

#define MAX_NGLOBALS 1024

/* global var arrays in this file */
int iglobal_fortemplates=-1; /* this is used as an index to pick the set of globals in the array */
void (*lop_fortemplates[MAX_NGLOBALS])(tVarList *, tVarList *, tVarList *, tVarList *);
void (*precon_fortemplates[MAX_NGLOBALS])(tVarList *, tVarList *, tVarList *, tVarList *);
tVarList *r_fortemplates[MAX_NGLOBALS];
tVarList *x_fortemplates[MAX_NGLOBALS];
tVarList *c1_fortemplates[MAX_NGLOBALS];
tVarList *c2_fortemplates[MAX_NGLOBALS];
long int dim_fortemplates[MAX_NGLOBALS];
tSparseVector **Acol_fortemplates[MAX_NGLOBALS];

/* extern globals */
extern double *DiagMinv_JacobiPrecon; /* from wrappers_for_JacobiPrecon.c */


/* print global vars in this file */
void print_globals_fortemplates(void)
{
  printf("iglobal_fortemplates = %d\n", iglobal_fortemplates);
  printf("lop_fortemplates[iglobal_fortemplates] = %p\n", lop_fortemplates[iglobal_fortemplates]);
  printf("precon_fortemplates[iglobal_fortemplates] = %p\n", precon_fortemplates[iglobal_fortemplates]);
  printf("r_fortemplates[iglobal_fortemplates] = %p -> norm %g\n", r_fortemplates[iglobal_fortemplates], GridL2Norm(r_fortemplates[iglobal_fortemplates]));
  printf("x_fortemplates[iglobal_fortemplates] = %p -> norm %g\n", x_fortemplates[iglobal_fortemplates], GridL2Norm(x_fortemplates[iglobal_fortemplates]));
  printf("c1_fortemplates[iglobal_fortemplates] = %p -> norm %g\n", c1_fortemplates[iglobal_fortemplates], GridL2Norm(c1_fortemplates[iglobal_fortemplates]));
  printf("c2_fortemplates[iglobal_fortemplates] = %p -> norm %g\n", c2_fortemplates[iglobal_fortemplates], GridL2Norm(c2_fortemplates[iglobal_fortemplates]));
  printf("dim_fortemplates[iglobal_fortemplates] = %ld\n", dim_fortemplates[iglobal_fortemplates]);
}


/* y := alpha*A*x + beta*y, where A is matrix defined by lop */
int matvec(double *alpha, double *x, double *beta, double *y)
{
  /* compute r = A*x */
  COPY_ARRAY_INTO_VL(x, x_fortemplates[iglobal_fortemplates]);
  lop_fortemplates[iglobal_fortemplates](r_fortemplates[iglobal_fortemplates],
                                         x_fortemplates[iglobal_fortemplates],
                                         c1_fortemplates[iglobal_fortemplates],
                                         c2_fortemplates[iglobal_fortemplates]);

  /* set x_fortemplates[iglobal_fortemplates] = y */
  COPY_ARRAY_INTO_VL(y, x_fortemplates[iglobal_fortemplates]);

  /* set r=alpha*r + beta*x_fortemplates[iglobal_fortemplates] */
  vladd(r_fortemplates[iglobal_fortemplates], *alpha,
        r_fortemplates[iglobal_fortemplates], *beta,
        x_fortemplates[iglobal_fortemplates]);

  /* copy r into y */
  COPY_VL_INTO_ARRAY(r_fortemplates[iglobal_fortemplates], y);
  
  return 0;
}

/* matvectrans: */
/* y := alpha*A'*x + beta*y, where A' is transpose of A which is in 
   Acol_fortemplates[iglobal_fortemplates] */
int Acol_times_vec_trans(double *alpha, double *x, double *beta, double *y)
{
  int i;
  int ncols = (int) dim_fortemplates[iglobal_fortemplates];
  double *f = dmalloc(ncols); /* temp array */

  /* compute f = A'*x */
  SparseMatrixLines_times_vector(Acol_fortemplates[iglobal_fortemplates],
                                 ncols, x, f);
  /* y = alpha*f + beta*y */
  for(i=0; i<ncols; i++)
    y[i] = (*alpha) * f[i] + (*beta) * y[i];

  free(f);
  return 0;
}


/* Precon: solves M*x = b for x */
int psolveLEFT(double *x, double *b)
{
  /* speed up special Jacobi precon */
  if(precon_fortemplates[iglobal_fortemplates]==Jacobi_Preconditioner_from_DiagM)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_fortemplates[iglobal_fortemplates]; i++)
      x[i] = b[i]*DiagMinv_JacobiPrecon[i];
  }
  else if(precon_fortemplates[iglobal_fortemplates]==Preconditioner_I)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_fortemplates[iglobal_fortemplates]; i++)
      x[i] = b[i];
  }
  else /* generic case */
  {
    /* solve M*x = b for x, solution is x_fortemplates[iglobal_fortemplates],
       b is in r_fortemplates[iglobal_fortemplates] */
    COPY_ARRAY_INTO_VL(b, r_fortemplates[iglobal_fortemplates]);
    precon_fortemplates[iglobal_fortemplates](x_fortemplates[iglobal_fortemplates],
                        r_fortemplates[iglobal_fortemplates],
                        c1_fortemplates[iglobal_fortemplates],
                        c2_fortemplates[iglobal_fortemplates]);

    /* copy x_fortemplates[iglobal_fortemplates] into x */
    COPY_VL_INTO_ARRAY(x_fortemplates[iglobal_fortemplates], x);
  }
  return 0;
}

/* Precon: solves M'*x = b for x.  M' is transpose of procon matrix M */
int psolveLEFTtrans(double *x, double *b)
{
  /* speed up special Jacobi precon */
  if(precon_fortemplates[iglobal_fortemplates]==Jacobi_Preconditioner_from_DiagM)
  {
    int i;

    SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_fortemplates[iglobal_fortemplates]; i++)
      x[i] = b[i]*DiagMinv_JacobiPrecon[i];
  }
  else if(precon_fortemplates[iglobal_fortemplates]==Preconditioner_I)
  {
    int i;
    //SGRID_LEVEL4_Pragma(omp parallel for)
    for(i=0; i<dim_fortemplates[iglobal_fortemplates]; i++)
      x[i] = b[i];
  }
  else /* generic case */
    errorexit("psolveLEFTtrans works only for Jacobi_Preconditioner_from_DiagM");

  return 0;
}

/* Right Precon: solves M2*x = b for x */
int psolveRIGHT(double *x, double *b)
{
  /* psolveRIGHT is not implemented yet, so use  x = b */
  int dim = dim_fortemplates[iglobal_fortemplates];
  int i;
  for(i=0; i<dim; i++)  x[i] = b[i];
  return 0;
}

/* Right Precon: solves M2'*x = b for x.
   M2' is transpose of right procon matrix M2 */
int psolveRIGHTtrans(double *x, double *b)
{
  /* psolveRIGHTtrans is not implemented yet, so use  x = b */
  return psolveRIGHT(x, b);
}

/* use left/right precon */
int psolveQ(double *x, double *b, char *s, short *slen)
{
  if(s[0]=='R')  /* if s="RIGHT" */
    return psolveRIGHT(x, b);
  else /* if s="LEFT" */
    return psolveLEFT(x, b);
}

/* use left/right precon for transposed matrix */
int psolveQtrans(double *x, double *b, char *s, short *slen)
{
  if(s[0]=='R')  /* if s="RIGHT" */
    return psolveRIGHTtrans(x, b);
  else /* if s="LEFT" */
    return psolveLEFTtrans(x, b);
}

/* call GMRES from templates */
int templates_gmres_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  long int RESTRT;
  double *WORK;		long int LDW;
  double *H;		long int LDH;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  if(Getv("GridIterators_GMRES_restart", "max"))  RESTRT = N;
  else	RESTRT = Geti("GridIterators_GMRES_restart");
  if(RESTRT>N) RESTRT = N;
  LDW = (N) + 1; 
  LDH = (RESTRT+1) + 1;
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE gmres_: ");
  if(pr) printf("  templates_gmres_wrapper: itmax=%d tol=%.3e "
                "N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e  RESTRT=%ld LDH=%ld\n",
                itmax, tol, N, LDW, ITER, RESID, RESTRT, LDH);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*(RESTRT+4), sizeof(double));
  H = (double *) calloc(LDH*(RESTRT+2), sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL || H==NULL)
    errorexit("templates_gmres_wrapper: out of memory for X, B, WORK, H\n"
              "  consider reducing GridIterators_GMRES_restart");

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call gmres from templates */
  gmres_(&N, B, X, &RESTRT, WORK, &LDW, H, &LDH, &ITER, &RESID,
          matvec, psolveLEFT, &INFO);
#else
  COMPILETEMPLATES("gmres");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out x and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);
  free(H);

  if(pr) printf("  templates_gmres_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER gmres_: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}


/* call bicgstab_ from templates */
int templates_bicgstab_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  double *WORK;		long int LDW;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = (N) + 1; 
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE bicgstab_: ");
  if(pr) printf("  templates_bicgstab_wrapper: itmax=%d tol=%.3e "
                "N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e\n",
                itmax, tol, N, LDW, ITER, RESID);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*7, sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_bicgstab_wrapper: out of memory for X, B, WORK");

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call bicgstab from templates */
  bicgstab_(&N, B, X, WORK, &LDW, &ITER, &RESID, matvec, psolveLEFT, &INFO);
#else
  COMPILETEMPLATES("bicgstab");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out vlx and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  if(pr) printf("  templates_bicgstab_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER bicgstab_: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}


/* call CGS from templates */
int templates_cgs_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  double *WORK;		long int LDW;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = (N) + 1; 
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) printf("  templates_cgs_wrapper: itmax=%d tol=%.3e "
                "N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e\n",
                itmax, tol, N, LDW, ITER, RESID);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*7, sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_cgs_wrapper: out of memory for X, B, WORK");

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call cgs from templates */
  cgs_(&N, B, X, WORK, &LDW, &ITER, &RESID, matvec, psolveLEFT, &INFO);
#else
  COMPILETEMPLATES("cgs");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out vlx and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  if(pr) printf("  templates_cgs_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}


/* QMR needs the transpose of the matrix A appearing in A x = b */
/* so get A from lop: */
void SetMatrixColumns_And_InvDiag(
       tSparseVector **Acol, tSparseVector **AcolFD, double *DiagAinv, int ncols,
       void (*lop)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
       tVarList *r, tVarList *x, 
       tVarList *c1,tVarList *c2, int pr)
{
  tGrid *grid = r->grid;
  tGrid *grid_bak;
  int bi;
  int col, ent;
  tSparseVector **AcolD;
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");

  if(pr) printf("SetMatrixColumns_And_InvDiag\n");

  /* set Acol */                
  SetMatrixColumns_slowly(Acol, lop, r, x, c1, c2, pr);
  if(pr&&0) 
    for(col=0; col<ncols; col++) prSparseVector(Acol[col]);

  /* set another matrix of we want fd precon */
  if(use_fd)
  {
    /* save current grid in grid_bak and then convert grid to fin. diff. */
    grid_bak = make_empty_grid(grid->nvariables, 0);
    copy_grid_withoutvars(grid, grid_bak, 0);
    convert_grid_to_fd(grid);
    if(pr) printf("Using finite differencing to set matrix DiagAinv...\n");

    /* set AcolFD */                
    SetMatrixColumns_slowly(AcolFD, lop, r, x, c1, c2, pr);
    if(pr&&0) 
      for(col=0; col<ncols; col++) prSparseVector(AcolFD[col]);

    /* restore grid to spectral */
    copy_grid_withoutvars(grid_bak, grid, 0);
    free_grid(grid_bak);
    AcolD = AcolFD;
  }
  else /* use spectral matrix */
    AcolD = Acol;

  /* now set DiagAinv to diagonal of matrix AcolD */
  for(col=0; col<ncols; col++)
    for(ent = 0; ent < AcolD[col]->entries; ent++)
      if(AcolD[col]->pos[ent] == col) 
      {
        double DiagA = AcolD[col]->val[ent];
        if(DiagA==0.0) errorexit("DiagA is singular!!!");
        DiagAinv[col] = 1.0/DiagA;
        break;
      }
}

/* call QMR from templates */
int templates_qmr_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int nvars=b->n;
  int ncols;
  double *DiagMinv_JacobiPrecon_sav;
  tSparseVector **Acol;
  tSparseVector **AcolFD;
  double *DiagAinv;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  double *WORK;		long int LDW;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = (N) + 1; 
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) printf("  templates_qmr_wrapper: itmax=%d tol=%.3e "
                "N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e\n",
                itmax, tol, N, LDW, ITER, RESID);

  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*11, sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_qmr_wrapper: out of memory for X, B, WORK");

  /* QMR needs the transpose of the matrix A appearing in A x = b */
  /* so get A from lop: */
  ncols = (int) N;
  /* allocate Acol */
  Acol = AllocateSparseVectorArray(ncols);
  if(Acol) { if(pr) printf("allocated %d matrix columns for Acol\n", ncols); }
  else       errorexit("no memory for Acol");
  /* allocate AcolFD to hold matrix of fd version of lop */
  AcolFD = AllocateSparseVectorArray(ncols);
  if(AcolFD) { if(pr) printf("allocated %d matrix columns for AcolFD\n", ncols); }
  else       errorexit("no memory for AcolFD");
  /* allocate memory for diagonal of matrix in DiagAinv */
  DiagAinv = calloc(ncols, sizeof(*DiagAinv));

  /* set the matrix Acol, AcolFD, DiagAinv */
  SetMatrixColumns_And_InvDiag(Acol, AcolFD, DiagAinv, ncols, 
                               lop, r, x ,c1,c2, pr);

  /* free AcolFD, because for now we do not need it */
  FreeSparseVectorArray(AcolFD, ncols);
  AcolFD=NULL;
  /* save DiagMinv_JacobiPrecon and then set it to DiagAinv */
  DiagMinv_JacobiPrecon_sav = DiagMinv_JacobiPrecon;
  DiagMinv_JacobiPrecon = DiagAinv;

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;
  Acol_fortemplates[iglobal_fortemplates]	= Acol;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call qmr from templates */
  qmr_(&N, B, X, WORK, &LDW, &ITER, &RESID,
       matvec, Acol_times_vec_trans, psolveQ, psolveQtrans, &INFO);
#else
  COMPILETEMPLATES("qmr");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out vlx and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  /* restore DiagMinv_JacobiPrecon and free matrices */
  DiagMinv_JacobiPrecon = DiagMinv_JacobiPrecon_sav;
  FreeSparseVectorArray(Acol, ncols);
  FreeSparseVectorArray(AcolFD, ncols);
  free(DiagAinv);

  if(pr) printf("  templates_qmr_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}

/* call BiCG from templates */
int templates_bicg_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int nvars=b->n;
  int ncols;
  double *DiagMinv_JacobiPrecon_sav;
  tSparseVector **Acol;
  tSparseVector **AcolFD;
  double *DiagAinv;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  double *WORK;		long int LDW;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = (N) + 1; 
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) printf("  templates_bicg_wrapper: itmax=%d tol=%.3e "
                "N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e\n",
                itmax, tol, N, LDW, ITER, RESID);

  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(LDW*6, sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_bicg_wrapper: out of memory for X, B, WORK");

  /* QMR needs the transpose of the matrix A appearing in A x = b */
  /* so get A from lop: */
  ncols = (int) N;
  /* allocate Acol */
  Acol = AllocateSparseVectorArray(ncols);
  if(Acol) { if(pr) printf("allocated %d matrix columns for Acol\n", ncols); }
  else       errorexit("no memory for Acol");
  /* allocate AcolFD to hold matrix of fd version of lop */
  AcolFD = AllocateSparseVectorArray(ncols);
  if(AcolFD) { if(pr) printf("allocated %d matrix columns for AcolFD\n", ncols); }
  else       errorexit("no memory for AcolFD");
  /* allocate memory for diagonal of matrix in DiagAinv */
  DiagAinv = calloc(ncols, sizeof(*DiagAinv));

  /* set the matrix Acol, AcolFD, DiagAinv */
  SetMatrixColumns_And_InvDiag(Acol, AcolFD, DiagAinv, ncols, 
                               lop, r, x ,c1,c2, pr);

  /* free AcolFD, because for now we do not need it */
  FreeSparseVectorArray(AcolFD, ncols);
  AcolFD=NULL;
  /* save DiagMinv_JacobiPrecon and then set it to DiagAinv */
  DiagMinv_JacobiPrecon_sav = DiagMinv_JacobiPrecon;
  DiagMinv_JacobiPrecon = DiagAinv;

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;
  Acol_fortemplates[iglobal_fortemplates]	= Acol;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call bicg from templates */
  bicg_(&N, B, X, WORK, &LDW, &ITER, &RESID,
        matvec, Acol_times_vec_trans, psolveLEFT, psolveLEFTtrans, &INFO);
#else
  COMPILETEMPLATES("bicg");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out vlx and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  /* restore DiagMinv_JacobiPrecon and free matrices */
  DiagMinv_JacobiPrecon = DiagMinv_JacobiPrecon_sav;
  FreeSparseVectorArray(Acol, ncols);
  FreeSparseVectorArray(AcolFD, ncols);
  free(DiagAinv);

  if(pr) printf("  templates_bicg_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}

/* call SOR from templates */
int templates_sor_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int pr = Getv("GridIterators_verbose", "yes");
  int i,j;
  long int N; /* dim of matrix */
  double *B;
  double *X;
  double OMEGA = Getd("GridIterators_SOR_omega");
  double *WORK;		long int LDW;
  long int ITER;
  double RESID;
  long int INFO=-1;
  double norm_b = GridL2Norm(b);

  /* set long int vars */
  N = 0 ;
  forallboxes(grid,i)  N += grid->box[i]->nnodes;
  N = (b->n) * N; 
  LDW = N; 
  ITER = itmax;
  RESID = tol;
  /* do we scale RESID? */
  if(Getv("GridIterators_templates_RESID_mode", "tol/norm(b)") && norm_b>0.0)
    RESID = RESID / norm_b;

  if(pr) prTimeIn_s("Time BEFORE sor_: ");
  if(pr) printf("  templates_sor_wrapper: itmax=%d tol=%.3e N=%ld LDW=%ld\n"
                "  ITER=%ld RESID=%.3e\n",
                itmax, tol, N, LDW, ITER, RESID);
  
  /* temporary storage */
  B = (double *) calloc(N, sizeof(double));
  X = (double *) calloc(N, sizeof(double));
  WORK = (double *) calloc(N*(N+3), sizeof(double));
  if(B==NULL || X==NULL || WORK==NULL)
    errorexit("templates_sor_wrapper: out of memory for X, B, WORK\n");

  /* in sor_:  work_dim1 = *ldw;   omega = work[work_dim1 + 1] */
  WORK[LDW+1] = OMEGA;

  /* increase global var index to store globals in new place in array */
  iglobal_fortemplates++;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);
  if(iglobal_fortemplates>=MAX_NGLOBALS) errorexit("increase MAX_NGLOBALS");
  fflush(stdout);

  /* setup global vars and functions needed in matvec and psolveLEFT */
  lop_fortemplates[iglobal_fortemplates]	= lop;
  precon_fortemplates[iglobal_fortemplates]	= precon;
  r_fortemplates[iglobal_fortemplates]		= r;
  x_fortemplates[iglobal_fortemplates]		= x;
  c1_fortemplates[iglobal_fortemplates]		= c1;
  c2_fortemplates[iglobal_fortemplates]		= c2;
  dim_fortemplates[iglobal_fortemplates]	= N;

  /* setup local B and X */
  COPY_VL_INTO_ARRAY(b, B);
  COPY_VL_INTO_ARRAY(x, X);

#ifdef TEMPLATES
  /* call sor from templates */
  sor_(&N, B, X, WORK,&LDW, &ITER,&RESID, matvec, templates_backsolve, &INFO);
#else
  COMPILETEMPLATES("sor");
#endif
  /* decrease global var index since globals are no longer needed */
  iglobal_fortemplates--;
  if(pr) printf("  iglobal_fortemplates=%d\n", iglobal_fortemplates);

  /* read out x and normres */
  COPY_ARRAY_INTO_VL(X, x);
  *normres = RESID;

  /* free temporary storage */
  free(B);
  free(X);
  free(WORK);

  if(pr) printf("  templates_sor_wrapper: ITER=%ld RESID=%.3e INFO=%ld\n",
                ITER, RESID, INFO);
  fflush(stdout);

  if(pr) prTimeIn_s("Time AFTER sor_: ");

  /* iteration failed */
  if(INFO<0) return INFO;
  if(INFO>0) return -ITER;
  
  /* success! */
  return ITER;
}


/*************************************************************************/
/* use one of the wrappers (e.g. templates_gmres_wrapper) as precon for  */
/* one of the linear solvers in this file (e.g. templates_gmres_wrapper) */
/* this precon is called by psolveLEFT above                                 */
/*************************************************************************/
void templates_Preconditioner_for_templates_solver(tVarList *vlx,
                                                   tVarList *vlr,
                                                   tVarList *vlc1,
                                                   tVarList *vlc2)
{
  tGrid *grid = vlx->grid;
  tGrid *grid_bak;
  int pr = Getv("GridIterators_verbose", "yes");
  int use_fd = Getv("GridIterators_Preconditioner_type", "fd");
  int itmax = Getd("GridIterators_Preconditioner_itmax");
  double rtol = Getd("GridIterators_Preconditioner_reltol");
  double tol, res;
  int its;
  int (*templates_wrapper)
           (tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

  if(Getv("GridIterators_templates_as_Preconditioner", "GMRES"))
    templates_wrapper = templates_gmres_wrapper;
  else if(Getv("GridIterators_templates_as_Preconditioner", "BICGSTAB"))
    templates_wrapper = templates_bicgstab_wrapper;
  else if(Getv("GridIterators_templates_as_Preconditioner", "CGS"))
    templates_wrapper = templates_cgs_wrapper;
  else
    errorexit("unkown GridIterators_templates_as_Preconditioner");

  if(iglobal_fortemplates<0)
    errorexit("templates_Preconditioner_for_templates_solver "
              "works only if called from psolveLEFT in wrappers_for_templates.c");

  /* set tol */
  tol = GridL2Norm(vlr)*rtol;
  if(tol==0.0) tol=rtol;

  if(pr) printf("  templates_Preconditioner_for_templates_solver:\n");        
  if(use_fd)
  {
    /* save current grid in grid_bak and then convert grid to fin. diff. */
    grid_bak = make_empty_grid(grid->nvariables, 0);
    copy_grid_withoutvars(grid, grid_bak, 0);
    convert_grid_to_fd(grid);
    if(pr) printf("  Using finite differencing...\n");
  }
  fflush(stdout);
  /* Note: the 2nd and 3rd arg can be the same, because templates_*_wrapper
     saves the 2nd arg immediately in its B. Of course vlr is overwritten! */
  its = templates_wrapper(vlx, vlr, vlr, vlc1,vlc2, itmax, tol, &res,
                          lop_fortemplates[iglobal_fortemplates],
                          Preconditioner_I);
  if(use_fd)
  {
    /* restore grid to spectral */
    copy_grid_withoutvars(grid_bak, grid, 0);
    free_grid(grid_bak);
  }
  if(pr)
    printf("  templates_Preconditioner_for_templates_solver: its=%d\n", its);
  fflush(stdout);
}
