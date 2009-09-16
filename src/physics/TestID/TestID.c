/* TestID.c */
/* Wolfgang Tichy 09/2009 */

#include "sgrid.h"
#include "TestID.h"



/* compute TestID data */
int TestID(tGrid *grid)
{
  /* allocate memory for psi and its derivs, and for Var List psiandderivs */
  psiandderivs->grid = grid;
  enablevarlist(psiandderivs);

  /* allocate memory for ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));

  printf("TestID: Computing initial data for testing:\n");
  printf("  TestID_type = %s\n", Gets("TestID_type"));

  if(Getv("TestID_type", "1dGaugeWave"))
    TestID_1dGaugeWave(grid, Ind("x"), Ind("gxx"), Ind("Kxx"), 
                       Ind("psi"), Ind("dpsiopsix"), Ind("ddpsiopsixx"), 
                       Ind("alpha"), Ind("betax") );
  else if(Getv("TestID_type", "1dLinearWave"))
    TestID_1dLinearWave(grid, Ind("x"), Ind("gxx"), Ind("Kxx"), 
                        Ind("psi"), Ind("dpsiopsix"), Ind("ddpsiopsixx"), 
                        Ind("alpha"), Ind("betax") );
  else
    errorexit("unknown TestID_type");
    
  /* set initial TrK for (needed) for some gauges in BSSN */
  set_K_initial(grid);

  /* set initial lapse and shift? */
  if (Getv("TestID_initial_lapse", "one"))
  {
    int bi;
    for(bi = 0; bi < grid->nboxes; bi++)
    {
      tBox *box = grid->box[bi];
      int i;
      double *alpha = box->v[Ind("alpha")];

      forallpoints(box, i)  alpha[i] = 1;
      printf("  set lapse to one\n");
    }
  }
  else 
    printf("  left lapse unchanged\n");

  if (Getv("TestID_initial_shift", "zero"))
  {
    int bi;

    enablevar(grid, Ind("betadotx"));
    for(bi = 0; bi < grid->nboxes; bi++)
    {
      tBox *box = grid->box[bi];
      int i;
      double *betax = box->v[Ind("betax")];
      double *betay = box->v[Ind("betay")];
      double *betaz = box->v[Ind("betaz")];
      double *betadotx = box->v[Ind("betadotx")];
      double *betadoty = box->v[Ind("betadoty")];
      double *betadotz = box->v[Ind("betadotz")];

      forallpoints(box, i)
      {
        betax[i] = betay[i] = betaz[i] = 0;
        betadotx[i] = betadoty[i] = betadotz[i] = 0;
      }
    }
    printf("  set shift to zero\n");
  }
  else
    printf("  left shift unchanged\n");

  return 0;
}
