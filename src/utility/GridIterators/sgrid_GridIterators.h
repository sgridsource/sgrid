/* sgrid_GridIterators.h */
/* Wolfgang Tichy 8/2008 */


/* macro that loops over vars, active boxes and subboxes */
#define forallVarsBoxesAndSubboxes_defIndices(vlv, blocki, vi,bi, sbi,sbj,sbk, nsb1,nsb2,nsb3) \
  for(blocki=0; (blocki)<(nsb1)*(nsb2)*(nsb3)*(activeboxes((vlv)->grid))*((vlv)->n); (blocki)++) { \
   int acu_, vi, abi_,bi, sbk, sbj, sbi; \
   vi   = (blocki)/((nsb1)*(nsb2)*(nsb3)*(activeboxes((vlv)->grid))); \
   acu_ = (blocki) - (nsb1)*(nsb2)*(nsb3)*(activeboxes((vlv)->grid)) * vi; \
   abi_ = acu_/((nsb1)*(nsb2)*(nsb3)); \
   bi   = IndexOf_ith_ActiveBox((vlv)->grid, abi_); \
   acu_ = acu_ - (nsb1)*(nsb2)*(nsb3) * abi_; \
   sbk  = acu_/((nsb1)*(nsb2)); \
   acu_ = acu_ - (nsb1)*(nsb2) * sbk; \
   sbj  = acu_/(nsb1); \
   sbi  = acu_ - (nsb1) * sbj;
/* to end the above marco use this: */
#define End_forallVarsBoxesAndSubboxes_defIndices }

/* gives index we need to loop over in a subbox, e.g.: i1 <= i < i2 */
#define IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3) \
  i1 = ((sbi)*n1)/(nsb1); \
  i2 = (((sbi)+1)*n1)/(nsb1); \
  j1 = ((sbj)*n2)/(nsb2); \
  j2 = (((sbj)+1)*n2)/(nsb2); \
  k1 = ((sbk)*n3)/(nsb3); \
  k2 = (((sbk)+1)*n3)/(nsb3); \


/* structure we use to store results of a Newton step */
typedef struct tNEWTONRESULTS {
  double lambda;
} tNewtonResults;

/* structure we use to store pointers to varlists and funcs */
typedef struct tNEWTONARGPOINTERS {
  int (**linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *));
} tNewtonArgPointers;

/* since use it often we define the type tLinSolver. This is a function
   pointer to one of our linear solvers */
typedef 
int (*tLinSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                  tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                  int lin_itmax, double lin_tol, double *lin_normres,
                  void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
                  void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *));
/* FIXME: use tLinSolver in more places in this file and also in Newton.c
          and elsewhere in the code */

/* return a funtion pointer to a linear solver based on its name */
tLinSolver get_linSolver_by_name(char *str);

/* dot product and L2-Norm over entire grid or box for varlists */
double GridDotProduct(tVarList *vlu, tVarList *vlv);
double GridL2Norm(tVarList *vlu);
double BoxL2Norm(int b, tVarList *vlu);
double varBoxL2Norm(tBox *box, int iu);


/* functions for grid iterations */

int Newton(
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlu, tVarList *vlFu, tVarList *vlc1, tVarList *vlc2,
  int itmax, double tol, double *normres, int pr,
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  void (*linPrecon)(tVarList *Hinv_v, tVarList *v, tVarList *, tVarList *),
  tVarList *vldu, tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int linSolv_itmax, double linSolv_tolFac, double linSolv_tol);

void Preconditioner_I(tVarList *vlJdu, tVarList *vldu,
                      tVarList *vlduDerivs, tVarList *vlu);
int Newton_linesrch(
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_c1, tVarList *vl_c2),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_d1, tVarList *vl_d2),
  tVarList *vlu, tVarList *vlFu, tVarList *vlc1, tVarList *vlc2,
  int itmax, double tol, double *normres, int pr,
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                   tVarList *vl_res, tVarList *vl_d1, tVarList *vl_d2,
                   int lin_itmax, double lin_tol, double *lin_normres,
	           void (*J_du)(tVarList *, tVarList *, tVarList *, tVarList *), 
	           void (*lin_precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  void (*Precon)(tVarList *Hinv_v, tVarList *v, tVarList *, tVarList *),
  tVarList *vldu, tVarList *vlres, tVarList *vld1, tVarList *vld2,
  int linSolv_itmax, double linSolv_tolFac, double linSolv_tol);


/* Below comes a long list of possible linear solvers. Most of them are
   called "wrapper" because they call external solvers */
/* BEGIN list of linear solvers */
int bicgstab(tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	     int imax, double tol, double *res,
	     void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* wrappers for funcs from templates */
int templates_gmres_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicgstab_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_cgs_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_qmr_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicg_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_sor_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
void templates_Preconditioner_for_templates_solver(tVarList *vlx,
                                                   tVarList *vlr,
                                                   tVarList *vlc1,
                                                   tVarList *vlc2);

/* wrappers for LAPACK */
int LAPACK_dgesv_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* wrappers for UMFPACK */
int bicgstab_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_gmres_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicgstab_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_cgs_wrapper_with_fd_UMFPACK_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int UMFPACK_solve_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int UMFPACK_solve_forSortedVars_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* wrappers for Jacobi Precon */
int bicgstab_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_gmres_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicgstab_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_cgs_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_qmr_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicg_wrapper_with_Jacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int bicgstab_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_gmres_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicgstab_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_cgs_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_qmr_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicg_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int ZIB_gmres_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int ZIB_gbit_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int ZIB_pcg_wrapper_with_BlockJacobi_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* from WTsolver.c */
int WTsolver(tVarList *vlx, tVarList *vlb, 
             tVarList *r, tVarList *c1,tVarList *c2,
	     int itmax, double tol, double *normres,
	     void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* from SOR.c */
int SOR_Iterator(tVarList *vlx, tVarList *vlb, 
                 tVarList *r, tVarList *c1,tVarList *c2,
    	         int itmax, double tol, double *normres,
	         void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	         void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int bicgstab_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_gmres_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_bicgstab_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int templates_cgs_wrapper_with_SOR_precon(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* wrappers for SuiteSparseQR */
int SuiteSparseQR_solve_wrapper(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* wrappers for ZIB_linSolvers */
int ZIB_gmres_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int ZIB_gbit_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
int ZIB_pcg_wrapper(
            tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
/* END list of linear solvers */
