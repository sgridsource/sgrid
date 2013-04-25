/* sgrid_main.h */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99 */

/* memory tracing */
// #define TRACEMEMORY

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876

/* Indices */
#define Index(i,j,k) ((i)+n1*((j)+n2*(k)))
#define Ind_n1n2(i,j,k,n1,n2) ((i)+(n1)*((j)+(n2)*(k)))
/* ijk = i + n1*j + n1*n2*k, thus:
   ijk/(n1*n2) = k 
   (ijk - n1*n2*k)/n1 = j
   (ijk - n1*n2*k - n1*j ) = i   */
#define kOfInd_n1n2(ijk,n1,n2)        ((ijk)/((n1)*(n2)))
#define jOfInd_n1n2_k(ijk,n1,n2,k)    (((ijk) - (n1)*(n2)*(k))/(n1))
#define iOfInd_n1n2_jk(ijk,n1,n2,j,k) ((ijk) - (n1)*(n2)*(k) - (n1)*(j))
#define forallijk(i,j,k) \
  for (k = 0; k < n3; k++) \
  for (j = 0; j < n2; j++) \
  for (i = 0; i < n1; i++)    


/* NOTE: In C99 these two have the same effect:
   #pragma omp parallel for
   _Pragma ( "omp parallel for" )
*/
/* To parallelize with OpenMP we need _Pragma ( "omp parallel for" ) 
   in many places. But for different applications we want to switch 
   them on or off depending on where they are.
   SGRID_LEVEL2_Pragma used for omp loops over a plane in a box (2d)
   SGRID_LEVEL3_Pragma used for omp loops over all points in a box (3d)
   SGRID_LEVEL4_Pragma used for omp loops over all boxes
   SGRID_LEVEL6_Pragma used for 6d omp loops (e.g. loop over box while interpolating onto each point)
   more can be defined easily.
   */
/* define SGRID_LEVEL2_Pragma macros that allow us to include
   certain pragmas only if certain things like LEVEL2_Pragmas are defined */
#ifdef LEVEL2_Pragmas
#define SGRID_LEVEL2_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL2_Pragma(x)
#endif

#ifdef LEVEL3_Pragmas
#define SGRID_LEVEL3_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL3_Pragma(x)
#endif

#ifdef LEVEL4_Pragmas
#define SGRID_LEVEL4_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL4_Pragma(x)
#endif

#ifdef LEVEL6_Pragmas
#define SGRID_LEVEL6_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6_Pragma(x)
#endif

#ifdef TOPLEVEL_Pragmas
#define SGRID_TOPLEVEL_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_TOPLEVEL_Pragma(x)
#endif


/* snap effect for grid coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))
#define dgreater(a,b) ((a)>(b)+dequaleps)

/* approx. <= and >= with dequaleps tolerance */
#define dlesseq(a,b) ( (a)<(b)+dequaleps )
#define dgreatereq(a,b) ( (a)>(b)-dequaleps )

#define signum(v) ((v) > 0.0 ? (1.0) : ((v) < 0.0 ? (-1.0) : (0.0)))


/* parameters.c */
/* parameter data base structure */
typedef struct {
  char *name;
  char *value;
  char *description;
  double numericalvalue; /* some pars are pure numbers, we cache them here */
  int booleanvalue; /* some pars are true/false, we cache them here as 1/0 */
} tParameter;
void makeparameter(char *name, char *value, char *description);
void AddPar(char *name, char *value, char *description);
void AddOrModifyPar(char *name, char *value, char *description);
void Sets(char *name, char *value);
void Seti(char *name, int i);
void Setd(char *name, double d);
char *Gets(char *name);
char *GetsLax(char *name);
int Geti(char *name);
double Getd(char *name);
int Getv(char *name, char *value);
int GetvLax(char *name, char *value);
char *NextEntry(char *list);
void Appends(char *name, char *value);
char *GetsInd(int i);
char *GetnameInd(int i);
tParameter *GetPointerTo_pbd();
int GetnParameters();
double GetCachedNumValByParIndex(int i);
int GetCachedBoolValByParIndex(int i);
int GetParIndex(char *name);
void print_pdb_i1_i2(tParameter *pdb, int i1, int i2, int pr_ind, int pr_cache);
void print_parameter_database();
void create_copy_of_pdb1_in_pdb2(tParameter *pdb1, int npdb1, int npdb1max,
                                 tParameter **pdb2);
tParameter *make_empty_pdb(int npdb1max);
void copy_pdb(tParameter *pdb1, int npdb1, tParameter *pdb2);
void free_pdb(tParameter *pdb1, int npdb1);

/* skeleton.c */
enum
{
  PRE_GRID,
  POST_GRID,
  PRE_INITIALDATA,
  INITIALDATA,
  POST_INITIALDATA,
  PRE_EVOLVE,
  EVOLVE,
  PRE_POST_EVOLVE,
  POST_EVOLVE,
  ANALYZE, 
  OUTPUT,
  POST_OUTPUT,
  POST_FINALIZE_GRID,
  NFUNCTIONS
};
void AddFun(int step, int (*f)(tGrid *), char *name);
void RunFun(int step, tGrid *grid);

/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(char *tensorindices, int *nilist, char **ilist, int *sym);

/* utilities.c */
void errorexit(char *file, int line, char *s);
void errorexits(char *file, int line, char *s, char *t);
void errorexiti(char *file, int line, char *s, int i);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))
#define errorexiti(s,i) errorexiti(__FILE__, __LINE__, (s), (i))

void yo(void);
void Yo(double x);
void prdivider(int n);
double getTimeIn_s();
void prTimeIn_s(char *comment);
double min2(double x, double y);
double min3(double x, double y, double z);
double max2(double x, double y);
double max3(double x, double y, double z);
double min_in_1d_array(double *f, int n, int *imin);
double max_in_1d_array(double *f, int n, int *imax);
double min2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imin);
double max2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imax);
double min3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imin);
double max3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imax);
int system2(char *s1, char *s2);
int system3(char *s1, char *s2, char *s3);
double *dmalloc(int n);
int *imalloc(int n);
char *cmalloc(int n);
void *pmalloc(int n);

/* variables.c */
#define PtrFromInd(box,i)     ((box)->v[(i)])
#define vlldataptr(vl,box,i)  ((box)->v[(vl)->index[i]])
extern int globalnvariables;
int Ind(char *name);
int IndLax(char *name);
void AddVar(char *name, char *indices, char *description);
void AddConstantVar(char *name, char *tensorindices, char *description);
tVarList *AddDuplicate(tVarList *vl, char *postfix);
tVarList *AddDuplicateEnable(tVarList *vl, char *postfix);

char *VarName(int i);
int VarNComponents(int i);
int VarComponent(int i);
int IndComponent0(int i);
char *VarNameComponent0(char *name);
char *VarTensorIndices(int i);
void VarNameSetBoundaryInfo(char *name, 
			    double farlimit, double falloff, double propspeed);
double VarFallOff(int i);
double VarFarLimit(int i);
double VarPropSpeed(int i);
int VarSymmetry(int i, int dir);
void VarNameSetConstantFlag(char *name);
int VarConstantFlag(int i);

void prvarlist(tVarList *v);
void prvarlist_inbox(tBox *box, tVarList *v);
tVarList *vlalloc(tGrid *grid);
void vlenable(tVarList *v);
void vlenablegrid(tGrid *grid, tVarList *v);
void vldisable(tVarList *v);
void vlfree(tVarList *u);
void vlpushone(tVarList *v, int vi);
void vlpush(tVarList *v, int vi);
void vlpushvl(tVarList *v, tVarList *u);
void vldropone(tVarList *v, int vi);
void vldrop(tVarList *v, int vi);
void vldropn(tVarList *v, int n);
tVarList *vlduplicate(tVarList *v);
void vlsetconstant(tVarList *u, const double c);
void vlcopy(tVarList *v, tVarList *u);
void vlcopygrid(tGrid *grid, tVarList *v, tVarList *u);
void varcopy(tGrid *grid, int iv, int iu);
void vlswap(tVarList *v, tVarList *u);
void varswap(tGrid *grid, int iv, int iu);
void vlaverage(tVarList *r, tVarList *a, tVarList *b);
void vlsubtract(tVarList *r, tVarList *a, tVarList *b);
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList*b);
void varadd(tGrid *grid, int ir, double ca, int ia, double cb, int ib);
void vladdto(tVarList *r, const double ca, tVarList *a);

tVarList *VLPtrEnable1(tGrid *grid, char *varname);
void VLDisableFree(tVarList *vl);
