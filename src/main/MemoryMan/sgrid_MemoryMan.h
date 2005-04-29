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
  double (*x_of_X[4])(void *aux, double X, double Y, double Z);	    /* func to compute x from X */
  double (*dX_dx[4][4])(void *aux, double x, double y, double z);	    /* dX_l_{spec}/dx_m_{cart} */
  void (*Sing_d_dx[4])(void *aux, void *vx, void *v1, void *v2, void *v3);  /* func to compute d/dx_m_{cart} at singular points */
  // double (*dX_dxdx[4][4][4])(void *aux, double x, double y, double z); /* d^2X_l_{spec}/(dx_m_{cart} dx_n_{cart}) */
  double bbox[6];	/* global bounding box */
  int ibbox[6];	/* global bounding box in index range */
} tBox;

/* several boxs make up a numerical grid */
typedef struct tGRID {
  tBox **box;	/* list of pointers to boxs */
  int nboxes;	/* number of boxes */
  int nvariables;	/* number of variables */
  int iteration;	/* number of iterations */
  double time;	/* physical time */
  double dt;	/* time step */
} tGrid;




/**************************************************************************/

/* variable lists 
   there should be a better place for it */
typedef struct tVARLIST {
  struct tGRID *grid;  /* struct tBOX *box;  is not needed in here */
  int n;
  int *index;
} tVarList;


/**************************************************************************/
/* loops */

#include "sgrid_MemoryMan_loops.h"


/**************************************************************************/
/* functions */

/* grid.c */
tGrid *make_grid(int pr);


/* print.c */
void printgrid(tGrid *g);
void printbox(tBox *l);
void printnode(tNode *n);
void printbbox(tBox *box, double *bbox, int *ibbox);
void printboundary(tBox *box);
void printvar_inbox(tBox *box, char *name);
void printvar(tGrid *grid, char *name);
void printVarList(tVarList *vl);


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

