/* Poisson.h */
/* Wolfgang Tichy, April 2005 */


int Poisson_startup(tGrid *grid);
int Poisson_analyze(tGrid *grid);
int Poisson_solve(tGrid *grid);

void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_Poisson(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);
void Precon_I(tVarList *vlJdu, tVarList *vldu,
              tVarList *vlduDerivs, tVarList *vlu);


int LinSolve_withLAPACK(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

int LinSolve_withUMFPACK(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
