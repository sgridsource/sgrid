/* evolve_test.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 10/03 */

/* test time integrators by integrating ODE 
   del_t u = -u,  u(0) = 1    =>   u(t) = exp(-t) 
*/

#include "sgrid.h"
#include "evolve.h"




/* evolution step for testing */
void evolve_test_evolve(tVarList *vlunew, 
			tVarList *vlupre, double c, tVarList *vlucur)
{
  tGrid *grid = vlucur->grid;
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];

    double *ucur = vlldataptr(vlucur, box, 0);
    double *upre = vlldataptr(vlupre, box, 0);
    double *unew = vlldataptr(vlunew, box, 0);
    int i;
  
    if (c == 0) 
      forallpoints(box, i) unew[i] = -ucur[i];
    else
      forallpoints(box, i) unew[i] = upre[i] - c * ucur[i];

printvar(grid, "evolve_test_u");
//printvar(grid, "evolve_test_u_p");
//printvar(grid, "evolve_test_u_q");
//printvar(grid, "evolve_test_u_c");
printvar(grid, "evolve_test_error");

  }
}




/* initialize test */
int evolve_test_startup(tGrid *grid) 
{
  tVarList *vlu = VLPtrEnable1(grid, "evolve_test_u");
  int b;

  /* I need to enable evolve_test_error  */
  enablevar(grid, Ind("evolve_test_error"));

  for (b = 0; b < grid->nboxes; b++)
  {  
    tBox *box = grid->box[b];
    int i;
    double *u = box->v[Ind("evolve_test_u")];
//double *u;
//printgrid(grid);
//printbox(box);
//printVarList(vlu);
//printf("Ind(\"evolve_test_u\")=%d\n", Ind("evolve_test_u")); yo();

    printf("Testing time integrators by integrating ODE:\n");

    forallpoints(box, i) u[i] = 1;

    evolve_vlregister(vlu);
    evolve_rhsregister(evolve_test_evolve);
  }
  return 0;
}




/* compute normalized absolute error in ANALYSIS */
int evolve_test_analyze(tGrid *grid)
{
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *u = box->v[Ind("evolve_test_u")];
    double *e = box->v[Ind("evolve_test_error")];
    double c;
    int i;

    /* computer error */
    forallpoints(box, i) 
      e[i] = fabs(u[i]/exp(-grid->time) - 1);

    /* scale error according to some expected order of convergence */
    c = pow(grid->dt, Getd("evolution_method_order"));
    forallpoints(box, i)
      e[i] /= c;
  }
  return 0;
}
