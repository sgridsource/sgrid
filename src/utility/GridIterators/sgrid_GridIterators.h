/* sgrid_GridIterators.h */
/* Wolfgang Tichy 8/2008 */

int bicgstab(tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
	     int imax, double tol, double *res,
	     void (*lop)(tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *));

int Newton(
  void  (*Fu)(tVarList *vl_Fu,  tVarList *vl_u,  tVarList *vl_aux),
  void (*Jdu)(tVarList *vl_Jdu, tVarList *vl_du, tVarList *vl_aux2),
  tVarList *vlu, tVarList *vlFu, tVarList *vlaux, 
  int itmax, double tol,double *normres, int pr,
  int (*linSolver)(tVarList *vl_du, tVarList *vl_Fu, 
                    tVarList *vl_res, tVarList *vl_aux2,
                    int lin_itmax, double lin_tol, double *lin_normres,
	            void (*J_du)(tVarList *, tVarList *, tVarList *), 
	            void (*lin_precon)(tVarList *, tVarList *, tVarList *)),
  void (*linPrecon)(tVarList *Hinv_v, tVarList *v, tVarList *),
  tVarList *vldu, tVarList *vlres, tVarList *vlaux2,
  int linSolv_itmax, double linSolv_tolFac );
