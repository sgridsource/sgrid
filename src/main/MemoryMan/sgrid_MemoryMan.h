/* sgrid_MemoryMan.h */
/* Wolfgang Tichy, April 2005 */

/* number of attributes a box can have */
#define NATTRIBS 10


/*************************************************************************/
/* auxilliary types and structures we need to define early */

/* struct that contains info about a coordinate singularity */
/* We store chars here:
   {'0', '.', 'i', 'n'} means {0, nonzero, +/-inf, nan }  */
typedef struct tSINGINFO {
  char x[4];        /* info about Cart. coords */
  char dx_dX[4][4]; /* and their derivs */
  int f[6];         /* face sing. is on (runs from 0 to 5 like bbox) */
} tSingInfo;

/* incomplete struct declaration so we can use struct tBOX * below */
struct tBOX;

/* struct that contains info that helps with coord trafo */
typedef struct tCOORDINFO {
  int iFS[6]; /* Index of var that defines a surface, e.g. FS[0]=sigma of
                 AnsorgNS. Both the vars iSurf[] and FSurf[] are computed from
                 iFS[] by functions in Coordinates. The user should only ever
                 set the var with index iFS[] */
  int iSurf[6];       /* index of surface at grid points, computed from iFS[] */
  int idSurfdX[6][4]; /* index of derivs at grid points, dSurfdX[0][1]=dFS[0]/dX */
  double s[6];        /* some values, e.g. value of surface var in case it is const */
  double xc[4];       /* xc[1..3] = (x,y,z) of coord center for this box */
  int dom;  /* domain index, e.g. 0-5 to indicate cubed sphere wedge */
  int type; /* coordinate type, e.g. outerCubedSphere */
  int useF; /* if useF=1 we use the funcs below to get values between grid
               points and to initialize vals inside iSurf and idSurfdX */
  double (*FSurf[6])(struct tBOX *box, int f, double C1, double C2); /* 6 funcs that return surface val, e.g. FSurf[0]=sigma */
  double (*dFSurfdX[6][4])(struct tBOX *box, int f, double C1, double C2); /* funcs that return derivs of FSurf, dFSurfdX[0][1]=dSurf[0]/dX */
} tCoordInfo;

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
  double (*x_of_X[4])(void *aux, int ind, double X, double Y, double Z);    /* func to compute x from X. Note: xyz_Of_XYZ would be better */
  int (*XYZ_Of_xyz)(struct tBOX *box, int ind, double x, double y, double z, double *X, double *Y, double *Z);  /* func to compute X,Y,Z from x,y,z */
  double (*dX_dx[4][4])(void *aux, int ind, double X, double Y, double Z);    /* dX_l_{spec}/dx_m_{cart} */
  void (*Sing_d_dx[4])(void *aux, void *v, void *v1, void *v2, void *v3);  /* func to compute d/dx_m_{cart} at singular points */
  double (*ddX_dxdx[4][4][4])(void *aux, int ind, double X, double Y, double Z); /* d^2X_l_{spec}/(dx_m_{cart} dx_n_{cart}) */
  double (*dx_dX[4][4])(void *aux, int ind, double X, double Y, double Z);	    /* dX_l_{cart}/dx_m_{spec} */
  int (*isSing)(void *aux, double X, double Y, double Z, int update, tSingInfo *si); /* indicates if there is a sing. at X,Y,Z */
  tCoordInfo CI[1];     /* info about coords, access e.g. as: box->CI->xc[1] */
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
  int TransformType1; /* If we use FFT and which FFT we use */
  int TransformType2; /* all 3 from set_TransformType_flags_inbox(tBox *box) */
  int TransformType3;
  int periodic[4];      /* periodic[1]=1 means direction1 is periodic */
  double bbox[6];	/* global bounding box */
  int ibbox[6];	        /* global bounding box in index range */
  int inactive;         /* if inactive=1 we skip this box in forallboxes, forallVarsBoxesAndSubboxes_defIndices, ... */
  int Attrib[NATTRIBS]; /* Each box can have attributes. They are usually set
                           and used inside certain modules/Projects. */
  struct tBFACE **bface;  // list of pointers to tBface
  int nbfaces;            // number of Bfaces on this box
} tBox;

/* several boxes make up a numerical grid */
typedef struct tGRID {
  tBox **box;	/* list of pointers to boxes */
  int nboxes;	/* number of boxes */
  int ninactive; /*  number of boxes with inactive=1 */ 
         /* sometimes e.g. Newton may need extra eqns that are not on any grid
            point. We can collect them in an extra box that we set inactive
            in certain calculations. Or Newton may need to skip some boxes. */
  int nvariables;	/* number of variables */
  int iteration;	/* number of iterations */
  double time;	/* physical time */
  double dt;	/* time step */
} tGrid;

/* point lists */
typedef struct tPL
{
  tGrid  *grid;  /* grid to which points belong */
  int *npoints;  /* npoints[b] = number of points in list in box b */
  int **point;   /* point[b] = array containing indices of all the points 
                    in the list in box b */
  int *blist;    /* list of boxes in which we have points */
  int nblist;    /* number of boxes in blist */
} tPointList;

/* Note: To include info about which boxes touch or overlap, tBox
         also contains info about BCs for boxes:
  struct tBFACE **bface;  // list of pointers to tBface
  int nbfaces;            // number of Bfaces on this box
where tBface is a part of a box face that touches or overlaps at most one
other box. We use the same BC on all of tBface. */
typedef struct tBFACE {
  tGrid *grid; // grid in which our boxface is
  int b;       // index of box on which Bface is located
  int f;       // face, runs from 0 to 5 like bbox (for each box)
  int fi;      // face index: 0 <= fi < nbfaces
  tPointList *fpts; // list of points on/near face, access w. forPointList_inbox
   // The normal vector is n^i_{a}=dx^i/dX^a, e.g. X^1=const face has n^i_{1}
   // dx^i/dX^a can be obtained from dX^a/dx^i using dXdx_from_dxdX
  int ob;          // ind. of other box that touches or overlaps, -1 if none
  int ofi;         // face index of other box (if they are touching)
  int oXi,oYi,oZi; // ind of vars in this box that contain coords in other box
  unsigned touch          : 1; // 1 if tBface only touches other face
  unsigned sameX          : 1; // 1 if X-coord of points in neighboring faces is same
  unsigned sameY          : 1; // 1 if Y-ccord of points in neighboring faces is same
  unsigned sameZ          : 1; // 1 if Z-ccord of points in neighboring faces is same
  unsigned same_fpts      : 1; // 1 if point i in fpts is same as point i in other fpts
  unsigned fpts_off_face  : 1; // 1 if points in fpts are not in face (e.g. ChebZeros)
  unsigned setnormalderiv : 1; // 1 if we set normal derivs of field and not field itself
  unsigned innerbound     : 1; // 1 if bface is inner boundary (e.g. horizon)
  unsigned outerbound     : 1; // 1 if bface is outer grid boundary (e.g. infinity)
} tBface;
/* The flags sameX/Y/Z refer to the X,Y,Z-coords of the neigboring box ob,
   not the coords of box b the bface is on! E.g. if the points on the
   touching face of box ob have the same Y and Z coords (of box ob) as the
   points in our bface from box b we set sameY=sameZ=1. This means once we
   find the X,Y,Z-coords (of box ob) of any point in our bface in box b,
   there is a coordinate line (with Y=const1, Z=const2) on which this point
   is located. Then we can later do 1D interpolation in box ob along its
   X-direction in some BC. */
/* NOTE: ob, ofi, oXi,oYi,oZi are all set to -1 when a bface is allocated with
   add_empty_bface. 
   So ofi=-1 means we do not know the other face index yet.
   In set_ofi_in_all_bfaces of utility/Coordinates/Coordinates_set_bfaces.c
   we set ofi<-1 if box b OVERLAPS with several faces of box ob. This makes
   sense, since for overlapping boxes a bface is not in the same place as
   another bface, so that BCs involve 3d interpolation. In this case we set
   ofi = -( (1<<(ofi1+1)) | (1<<(ofi2+1)) | ... )
   so that (assuming int ofi has enough space) each other face index is in
   one bit of ofi. We also make sure that then ofi<-1 and even.
   In set_ofi_in_all_bfaces we also set
   ofi = -( (1<<(ofi1+1)) | (1<<(ofi2+1)) | ... )
   if a bface of box b TOUCHES several bfaces from other boxes. BUT this
   should be avoided!!! When touching, bfaces should be split so that each
   bface only touches one other bface. We can, after all, have as many
   bfaces as we want on each boxface. */
/* ALSO: We need a function that sets oXi,oYi,oZi and sets the corresponding
   vars. Probably we just need the 3 vars oX,oY,oZ in each box. We can then
   loop over all bfaces of a box and set oX,oY,oZ to whatever they need
   to be. If edges or corners belong to several bfaces, the last bface wins.
   If the BC routine later loops in the same order over the bfaces, all
   should be consistent and work out. So use forallbfaces everywhere!!! */


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
/* use lists with int entries */
#define TYP int
#include "list_templates.h"
#undef TYP

/**************************************************************************/
/* functions */

/* grid.c */
tGrid *make_grid(int pr);
int set_BoxStructures_fromPars(tGrid *g, int pr);
tGrid *make_empty_grid(int nvariables, int pr);
int copy_grid_withoutvars(tGrid *g_old, tGrid *g_new, int pr);
int copy_grid(tGrid *g_old, tGrid *g_new, int pr);
int copy_gridvar(int vi, tGrid *g_old, tGrid *g_new);
int point_grid_tosamevars(tGrid *g_old, tGrid *g_new, int pr);
int set_gridvars_toNULL(tGrid *g_new, int pr);
void deactivate_boxes(tGrid *grid, intList *ilist);
void activate_boxes(tGrid *grid, intList *ilist);
void deactivate_allboxes(tGrid *grid);
void activate_allboxes(tGrid *grid);
int IndexOf_ith_ActiveBox(tGrid *grid, int i);

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
void printbface(tBface *bface);
void printbfaces(tBox *box);
void printCI(tBox *box);


/* storage.c */
tGrid *alloc_grid(int nboxes, int nvariables);
void make_box(tGrid *g, tBox *box, int l);
int add_empty_bface(tBox *box, int f);
int add_point_to_bface_inbox(tBox *box, int fi, int ijk, int f);
tBface *duplicate_bface_without_fpts_for_grid(tBface *bface0, tGrid *grid);
tBface *duplicate_bface_for_grid(tBface *bface0, tGrid *grid);
void free_bface(tBface *bface);
void free_all_bfaces(tBox *box);
int remove_bface(tBox *box, int fi);
int remove_bfaces_with_NULL_fpts(tBox *box);
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

/* boxlists.c */
void pr_boxlist(intList *bl);
int bladd_ifAttrib(tGrid *grid, int iAttr, int AttrVal, intList *bl);
int bladd_ifnotAttrib(tGrid *grid, int iAttr, int AttrVal, intList *bl);
int bl2_ifAttrib(tGrid *grid, intList *bl1, int iAttr, int AttrVal, intList *bl2);
int bl2_ifnotAttrib(tGrid *grid, intList *bl1, int iAttr, int AttrVal, intList *bl2);
int bladd_neighbors(tBox *box, intList *bl);
int bladd_neighbors_of_neighbors(tBox *box, intList *bl);

/* points.c */
tPointList *AllocatePointList(tGrid *grid);
void AddToPointList(tPointList *PL, int boxindex, int newpoint);
void DropFromPointList(tPointList *PL, int boxindex, int point);
void DropListFromPointList(tPointList *PL, int boxindex, intList *drop);
void FreePointList(tPointList *PL);
void prPointList(tPointList *PL);
void prPointList_inbox(tPointList *PL, int b);
void prPointList_ijk_inbox(tPointList *PL, int b);
tPointList *DuplicatePointList_for_grid(tPointList *PL0, tGrid *grid);
tPointList *DuplicatePointList(tPointList *PL0);
