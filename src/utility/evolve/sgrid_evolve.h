/* sgrid_evolve.h */
/* Wolfgang Tichy, April 2005 */

void evolve_rhsregister(void (*f)(tVarList *, tVarList *, double, tVarList *));
void evolve_vlregister(tVarList *u);
void evolve_vlretrieve(tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_pp);
void addDissipation(tVarList *vlu_n, tVarList *vlu_p, double dt, 
                    tVarList *vlu_c);
