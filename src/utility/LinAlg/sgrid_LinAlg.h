/* sgrid_LinAlg.h */
/* Wolfgang Tichy 8/2003 */


/* TYPES */

/* type def for Sparse Vectors */
typedef struct tSV
{
  int    entries; /* number of entries in list */
  int    *pos;    /* array containing positions */
  double *val;    /* val[i] contains value at position pos[i] */
} tSparseVector;


/* FUNCTIONS */
      
/* SparseVector_utils.c */
tSparseVector *AllocateSparseVector(void);
void AddToSparseVector(tSparseVector *SV, int newpos, double newval);
void FreeSparseVector(tSparseVector *SV);
void prSparseVector(tSparseVector *SV);
double GetSparseVectorComponent(tSparseVector *SV, int comp);

/* ReadMatrixFrom_Fx.c */
void SetMatrixLines_slowly(tSparseVector **Aline,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);
void SetMatrixColumns_slowly(tSparseVector **Acol,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr);

/* lapack_interface.c */
int lapack_dgesv(tSparseVector **Aline, tVarList *vlx, tVarList *vlb, int pr);

/* umfpack_interface.c */
int umfpack_solve(tSparseVector **Aline, tVarList *vlx, tVarList *vlb,
                  double dropbelow, int pr);
int umfpack_solve_fromAcolumns(tSparseVector **Acol,
                               tVarList *vlx, tVarList *vlb,
                               double dropbelow, int pr);
