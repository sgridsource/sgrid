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

#define forinnerijk(i,j,k, n1,n2,n3) \
  for(k = 1; k < (n3-1); k++) \
  for(j = 1; j < (n2-1); j++) \
  for(i = 1; i < (n1-1); i++)


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

#if defined(LEVEL6_Pragmas) || defined(TOPLEVEL_Pragmas)
#define SGRID_LEVEL6orTOP_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6orTOP_Pragma(x)
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
tParameter *GetPointerTo_pbd(void);
int GetnParameters(void);
double GetCachedNumValByParIndex(int i);
int GetCachedBoolValByParIndex(int i);
int GetParIndex(char *name);
int Set_pdb_iStart_AtPar(char *name);
void print_pdb_i1_i2(tParameter *pdb, int i1, int i2, int pr_ind, int pr_cache);
void print_parameter_database(void);
void create_copy_of_pdb1_in_pdb2(tParameter *pdb1, int npdb1, int npdb1max,
                                 tParameter **pdb2);
tParameter *make_empty_pdb(int npdb1max);
void copy_pdb(tParameter *pdb1, int npdb1, tParameter *pdb2);
void free_pdb(tParameter *pdb1, int npdb1);
void free_global_parameter_database_contents(void);

/* skeleton.c */
enum
{
  POST_PARAMETERS,
  PRE_GRID,
  POST_GRID,
  PRE_COORDINATES,
  COORDINATES,
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
  PRE_FINALIZE_GRID,
  POST_FINALIZE_GRID,
  LAST,
  NFUNCTIONS
};
void AddFun(int step, int (*f)(tGrid *), char *name);
void RunFun(int step, tGrid *grid);
void PrintFuncs(void);

/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(char *tensorindices, int *nilist, char **ilist, int *sym);

/* utilities.c */
void SGRID_errorexit(char *file, int line, char *s);
void SGRID_errorexits(char *file, int line, char *s, char *t);
void SGRID_errorexiti(char *file, int line, char *s, int i);
#define errorexit(s)    SGRID_errorexit(__FILE__, __LINE__, (s))
#define errorexits(s,t) SGRID_errorexits(__FILE__, __LINE__, (s), (t))
#define errorexiti(s,i) SGRID_errorexiti(__FILE__, __LINE__, (s), (i))
#define PRF     printf("%s", __func__)
#define PRFs(s) printf("%s%s", __func__, s)

void yo(void);
void Yo(double x);
void prdivider(int n);
void initTimeIn_s(void);
double getTimeIn_s(void);
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
int finit(double x);
int get_par_from_str(const char *str, char *name, const char *delim,
                     char *value, int n);
int copy_file_into_dir(char *fname, char *dir);
int system2(char *s1, char *s2);
int system3(char *s1, char *s2, char *s3);
int system_emu(const char *command);
int construct_argv(char *str, char ***argv);
double *dmalloc(int n);
int *imalloc(int n);
char *cmalloc(int n);
void *pmalloc(int n);

/* variables.c */
#define PtrFromInd(box,i)     ((box)->v[(i)])
#define vlldataptr(vl,box,i)  ((box)->v[(vl)->index[i]])
extern int globalnvariables;
int IndLax(char *name);
int Ind(char *name);
int Set_vdb_iStart_AtPar(char *name);
void AddVar(char *name, char *indices, char *description);
void AddConstantVar(char *name, char *tensorindices, char *description);
void AddVarToGrid(tGrid *grid, char *name, char *tensorindices,
                  char *description);
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

/* endianIO.c */
size_t SGRID_fwrite_swapbytes(const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t SGRID_fread_swapbytes(void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t SGRID_fwrite_little(const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t SGRID_fread_little(void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t SGRID_fwrite_big(const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t SGRID_fread_big(void *ptr, size_t size, size_t nmemb, FILE *fp);
int SGRID_return_BYTE_ORDER_LITTLE(void);
int SGRID_print_endian_info(tGrid *grid);
size_t fwrite_double_little(const double *buf, size_t nmemb, FILE *fp);
size_t fread_double_little(double *buf, size_t nmemb, FILE *fp);

/* sgrid_MPI.c */
int sgrid_MPI_Init(int *pargc, char ***pargv);
int sgrid_MPI_Finalize(void);
int sgrid_MPI_rank(void);
int sgrid_MPI_size(void);  
int sgrid_MPI_barrier(void);
