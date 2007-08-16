/* sgrid_main.h */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99 */

/* memory tracing */
// #define TRACEMEMORY

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.1415926535897932
#define PIh 1.5707963267948966
#define PIq 0.78539816339744831
#define Index(i,j,k) ((i)+n1*((j)+n2*(k)))
#define forallijk(i,j,k) \
  for (k = 0; k < n3; k++) \
  for (j = 0; j < n2; j++) \
  for (i = 0; i < n1; i++)    

/* snap effect for grid coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))


/* parameters.c */
void AddPar(char *name, char *value, char *description);
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

/* skeleton.c */
enum
{
  PRE_GRID,
  PRE_INITIALDATA,
  INITIALDATA,
  POST_INITIALDATA,
  PRE_EVOLVE,
  EVOLVE,
  POST_EVOLVE,
  ANALYZE, 
  OUTPUT,
  POST_OUTPUT,
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

void yo(), yo1(), yo2(), yo3(), yob();
void prdivider(int n);
double min2(double x, double y);
double min3(double x, double y, double z);
double max2(double x, double y);
double max3(double x, double y, double z);
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
tVarList *vlduplicate(tVarList *v);
void vlsetconstant(tVarList *u, const double c);
void vlcopy(tVarList *v, tVarList *u);
void vlcopybox(tBox *level, tVarList *v, tVarList *u);
void vlaverage(tVarList *r, tVarList *a, tVarList *b);
void vlsubtract(tVarList *r, tVarList *a, tVarList *b);
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList*b);

tVarList *VLPtrEnable1(tGrid *grid, char *varname);
void VLDisableFree(tVarList *vl);
