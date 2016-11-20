/* sgrid_MemoryMan.h */
/* Wolfgang Tichy, April 2005 */


/*************************************************************************/
/* basic structures */

/* a node */
typedef struct tNODE {
  int i;	/* index in box list */
  int type;	/* type of node */
} tNode;

/* a box */
typedef struct tBOX {
  struct tGRID *grid;	/* pointer to parent grid */
  int b;		/* index of this box */
  int n1;		/* number of points in x-direction */
  int n2;		/* number of points in y-direction */
  int n3;		/* number of points in z-direction */
  int nnodes;		/* number of nodes */
  tNode *node;		/* list of nodes */
  double **v;		/* list of data pointers to variables */
  double *D1;		/* differentiation matrix for direction 1 */
  double *D2;		/* differentiation matrix for direction 2 */
  double *D3;		/* differentiation matrix for direction 3 */
  double *DD1;		/* diff. matrix for two derivs in direction 1 */
  double *DD2;		/* diff. matrix for two derivs in direction 2 */
  double *DD3;		/* diff. matrix for two derivs in direction 3 */
  double (*x_of_X[4])(void *aux, int ind, double X, double Y, double Z);	    /* func to compute x from X */
  double (*dX_dx[4][4])(void *aux, int ind, double X, double Y, double Z);	    /* dX_l_{spec}/dx_m_{cart} */
  void (*Sing_d_dx[4])(void *aux, void *v, void *v1, void *v2, void *v3);  /* func to compute d/dx_m_{cart} at singular points */
  double (*ddX_dxdx[4][4][4])(void *aux, int ind, double X, double Y, double Z); /* d^2X_l_{spec}/(dx_m_{cart} dx_n_{cart}) */
  double *F1;		/* filter matrix for direction 1 */
  double *F2;		/* filter matrix for direction 2 */
  double *F3;		/* filter matrix for direction 3 */
  double (*basis1)(void *aux, double a, double b, int k, int n1, double X); /* basisfunction k in direction 1 */
  double (*basis2)(void *aux, double a, double b, int k, int n2, double X); /* basisfunction k in direction 2 */
  double (*basis3)(void *aux, double a, double b, int k, int n3, double X); /* basisfunction k in direction 3 */
  double *Mcoeffs1;	/* matrix to obtain coeffs in direction 1 */
  double *Mcoeffs2;	/* matrix to obtain coeffs in direction 2 */
  double *Mcoeffs3;	/* matrix to obtain coeffs in direction 3 */
  double *Meval1;	/* matrix to evaluate in direction 1 */
  double *Meval2;	/* matrix to evaluate in direction 2 */
  double *Meval3;	/* matrix to evaluate in direction 3 */
  double *Int1;		/* matrix to integrate in direction 1 */
  double *Int2;		/* matrix to integrate in direction 2 */
  double *Int3;		/* matrix to integrate in direction 3 */
  void (*get_coeffs1)(double *,double *, int); /* get coeffs in dir. 1 */
  void (*get_coeffs2)(double *,double *, int); /* get coeffs in dir. 2 */
  void (*get_coeffs3)(double *,double *, int); /* get coeffs in dir. 3 */
  void (*coeffs_of_deriv1)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_deriv2)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_deriv3)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_2ndderiv1)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_2ndderiv2)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_2ndderiv3)(void *aux, double, double, double *,double *, int);
  void (*coeffs_of_int1)(double, double, double *,double *, int);
  void (*coeffs_of_int2)(double, double, double *,double *, int);
  void (*coeffs_of_int3)(double, double, double *,double *, int);
  void (*eval_onPoints1)(double *,double *, int);
  void (*eval_onPoints2)(double *,double *, int);
  void (*eval_onPoints3)(double *,double *, int);
  void (*filter_coeffs1)(double *, int, int);
  void (*filter_coeffs2)(double *, int, int);
  void (*filter_coeffs3)(double *, int, int);
  int TransformType1;
  int TransformType2;
  int TransformType3;
  double bbox[6];	/* global bounding box */
  int ibbox[6];	/* global bounding box in index range */
} tBox;

/* several boxes make up a numerical grid */
typedef struct tGRID {
  tBox **box;	/* list of pointers to boxes */
  int nboxes;	/* number of boxes */
  int nvariables;	/* number of variables */
  int iteration;	/* number of iterations */
  double time;	/* physical time */
  double dt;	/* time step */
} tGrid;

/* Note: To include info about which boxes touch or overlap tGrid
could also contain info about transitions between boxes:
  tTransition **trans;  // list of pointers to transitions
  int ntrans;           // number of transitions
where
typedef struct tFACE {
  int b;  // index of box on which face is located
  int f;  // face index. Runs from 0 to 5 like bbox (for each box)
} tFace;
typedef struct tTRANSITION {
  tFace f0, f1;  // the 2 faces that coincide
                 // if f1.f<0, f0 coincides only with interior of box f1.b
} tTransition;
*/



/**************************************************************************/

/* variable lists 
   there should be a better place for it */
typedef struct tVARLIST {
  struct tGRID *grid;  /* struct tBOX *box;  is not needed in here */
  double time;
  int n;
  int *index;
  void *vlPars; /* A pointer that is usually NULL, but can point to some
                   object that contains special extra pars or info. This
                   pointer is not touched by the funcs in variables.c (such
                   as vlduplicate, vlcopy, vlfree, ...). So the user has 
                   to manage it: e.g. free it, before calling vlfree. */
} tVarList;


/**************************************************************************/
/* loops */

#include "sgrid_MemoryMan_loops.h"


/**************************************************************************/
/* functions */

/* grid.c */
tGrid *make_grid(int pr);
int set_BoxStructures_fromPars(tGrid *g, int pr);
tGrid *make_empty_grid(int nvariables, int pr);
int copy_grid_withoutvars(tGrid *g_old, tGrid *g_new, int pr);
int copy_grid(tGrid *g_old, tGrid *g_new, int pr);
int point_grid_tosamevars(tGrid *g_old, tGrid *g_new, int pr);
int set_gridvars_toNULL(tGrid *g_new, int pr);


/* print.c */
void printgrid(tGrid *g);
void printbox(tBox *l);
void printnode(tNode *n);
void printbbox(tBox *box, double *bbox, int *ibbox);
void printboundary(tBox *box);
void printvar_inbox(tBox *box, char *name);
void printvar(tGrid *grid, char *name);
void printVarList(tVarList *vl);
void printmatrix(double *M, int n);


/* storage.c */
tGrid *alloc_grid(int nboxes, int nvariables);
void make_box(tGrid *g, tBox *box, int l);
void free_box(tBox *box);
void free_grid_only(tGrid *g);
void free_grid(tGrid *g);
void enablevarcomp_inbox(tBox *box, int i);
void disablevarcomp_inbox(tBox *box, int i);
void enablevar_inbox(tBox *box, int i);
void disablevar_inbox(tBox *box, int i);
void enablevarcomp(tGrid *grid, int i);
void disablevarcomp(tGrid *grid, int i);
void enablevar(tGrid *grid, int i);
void disablevar(tGrid *grid, int i);
void enablevarlist(tVarList *v);
void disablevarlist(tVarList *v); 
void enablesamevars(tGrid *grid, tGrid *newgrid);
//double *PtrEnable(tBox *box, char *name);
//double *Ptr(tBox *box, char *name);

void realloc_gridvariables(tGrid *g, int nvariables);

/* points.c */
tPointList *AllocatePointList(tGrid *grid);
void AddToPointList(tPointList *PL, int boxindex, int newpoint);
void FreePointList(tPointList *PL);
void prPointList(tPointList *PL);
