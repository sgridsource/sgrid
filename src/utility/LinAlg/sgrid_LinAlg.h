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
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2);
