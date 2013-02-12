/* sgrid_LinAlg.h */
/* Wolfgang Tichy 8/2003 */



/* LongInt contains the long integer type used. For umfpack it must be 
   the same as SuiteSparse_long, which is a 64 bit integer on 64 
   bit linux, or any LP64 system */
#define LONGINT long


/* TYPES */

/* type def for Sparse Vectors */
typedef struct tSV
{
  int    entries; /* number of entries in list */
  int    *pos;    /* array containing positions */
  double *val;    /* val[i] contains value at position pos[i] */
} tSparseVector;

/* struct that contains all about a matrix that we need for umfpack_dl_solve */
typedef struct T_UMFPACK_A {
  LONGINT sys;       /* use sys=UMFPACK_A for a normal solve */
  LONGINT *Ap;       /* matrix */
  LONGINT *Ai;
  double  *Ax;
  void    *Numeric;  /* computed with umfpack_dl_numeric */
} tUMFPACK_A;

/* struct that contains all about a matrix, that we need for
   SuiteSparseQR_C_solve */
typedef struct T_SPQR_A {
  LONGINT sys;      /* default: SPQR_RX_EQUALS_B */
  LONGINT ordering; /* default: SPQR_ORDERING_DEFAULT */
  double  tol;      /* default: SPQR_DEFAULT_TOL */
  void *A;  /* covert to: cholmod_sparse *A; */
  void *QR; /* covert to: SuiteSparseQR_C_factorization *QR; */
  void *cc; /* covert to: cholmod_common *cc; */
} tSPQR_A;



/* FUNCTIONS */
      
/* SparseVector_utils.c */
tSparseVector *AllocateSparseVector(void);
void AddToSparseVector(tSparseVector *SV, int newpos, double newval);
void FreeSparseVector(tSparseVector *SV);
void prSparseVector(tSparseVector *SV);
tSparseVector **AllocateSparseVectorArray(int n);
void FreeSparseVectorArray(tSparseVector **A, int n);
double GetSparseVectorComponent(tSparseVector *SV, int comp);
void SparseMatrixLines_times_vector(tSparseVector **Aline, int nlines,
                                    double *x, double *f);
void vector_times_SparseMatrixLines(double *x, tSparseVector **Aline,
                                    int nlines, double *f);
void SparseMatrixLinesTranspose_times_vector(tSparseVector **Aline, int nlines,
                                             double *x, double *f);
 
/* ReadMatrixFrom_Fx.c */
void SetMatrixLines_slowly(tSparseVector **Aline,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixLines_forSortedVars_slowly(tSparseVector **Aline,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixColumns_slowly(tSparseVector **Acol,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixColumns_forSortedVars_slowly(tSparseVector **Acol,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixColumns_ForOneVarInOneBox_slowly(tSparseVector **Acol,
    int vlind, int b,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixColumns_ForOneVarInOneSubBox_slowly(tSparseVector **Acol,
    int vlind, int b, 
    int sbi, int sbj, int sbk,  int nsb1, int nsb2, int nsb3,
    void  (*Fx)(tVarList *Fdx, tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);

/* lapack_interface.c */
int lapack_dgesv(tSparseVector **Aline, tVarList *vlx, tVarList *vlb, int pr);

/* umfpack_di_interface.c */
void allocate_umfpack_di_matrix(int **Ap, int **Ai, double **Ax, int n, int nz);
void free_umfpack_di_matrix(int *Ap, int *Ai, double *Ax);
int set_umfpack_di_matrix_from_lines(int *Ap, int *Ai, double *Ax,
                                     tSparseVector **Aline, int nlines,
                                     double dropbelow, int pr);
int set_umfpack_di_matrix_from_columns(int *Ap, int *Ai, double *Ax,
                                       tSparseVector **Acol, int ncols,
                                       double dropbelow, int pr);
int umfpack_di_solve_fromAlines(tSparseVector **Aline, tVarList *vlx, tVarList *vlb,
                                double dropbelow, int pr);
int umfpack_di_solve_fromAcolumns(tSparseVector **Acol,
                                  tVarList *vlx, tVarList *vlb,
                                  double dropbelow, int pr);
int umfpack_di_solve_forSortedVars_fromAcolumns(tSparseVector **Acol,
      tVarList *vlx, tVarList *vlb,
      double dropbelow, int pr);
int umfpack_di_solve_from_Ap_Ai_Ax(int *Ap, int *Ai, double *Ax,
                                   tVarList *vlx, tVarList *vlb, int pr);

/* umfpack_dl_interface.c */
void allocate_umfpack_dl_matrix(LONGINT **Ap, LONGINT **Ai, double **Ax,
                                LONGINT n, LONGINT nz);
void free_umfpack_dl_matrix(LONGINT *Ap, LONGINT *Ai, double *Ax);
LONGINT set_umfpack_dl_matrix_from_lines(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                         tSparseVector **Aline, LONGINT nlines,
                                         double dropbelow, int pr);
LONGINT set_umfpack_dl_matrix_from_columns(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                           tSparseVector **Acol, LONGINT ncols,
                                           double dropbelow, int pr);
int umfpack_dl_solve_fromAcolumns(tSparseVector **Acol,
                                  tVarList *vlx, tVarList *vlb,
                                  double dropbelow, int pr);
int umfpack_dl_solve_from_Ap_Ai_Ax(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                   tVarList *vlx, tVarList *vlb, int pr);
int umfpack_dl_solve_from_Ap_Ai_Ax_x_b(LONGINT *Ap, LONGINT *Ai, double *Ax,
                                       double *x, double *b, LONGINT nrows,
                                       int pr);
int umfpack_dl_numeric_from_tUMFPACK_A(tUMFPACK_A *umfpackA,
                                       LONGINT nrows, int pr);
int umfpack_dl_solve_from_tUMFPACK_A_x_b(tUMFPACK_A umfpackA,
                                         double *x, double *b, int pr);

/* SuiteSparseQR_C_interface.c */
int SuiteSparseQR_solve_fromAcolumns(tSparseVector **Acol,
                                     tVarList *vlx, tVarList *vlb,
                                     double dropbelow, int pr);
int allocate_and_init_tSPQR_A_struct(tSPQR_A *SPQR, LONGINT nrows, 
                                     LONGINT nz, int pr);
int free_tSPQR_A_struct(tSPQR_A *SPQR_A);
int SuiteSparseQR_C_factorize_tSPQR_A(tSPQR_A *SPQR_A, int pr);
int SuiteSparseQR_solve_from_tSPQR_A_x_b(tSPQR_A SPQR_A,
                                         double *x, double *b, int pr);
