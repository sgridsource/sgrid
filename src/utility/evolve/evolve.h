/* evolve.h */
/* Wolfgang Tichy, April 2005 */


extern tVarList *u_c, *u_p, *u_q, *u_r;

int evolve(tGrid *grid);
int evolve_test_startup(tGrid *grid);
int evolve_test_analyze(tGrid *grid);

void evolve_icn(tGrid *grid);
void evolve_euler(tGrid *grid);
void evolve_rk(tGrid *grid);

void (*evolve_rhs)
  (tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void (*evolve_algebraicConditions)(tVarList *unew, tVarList *upre);
