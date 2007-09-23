/* rk.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 10/03 */

/* First version in admc_rk.c for Cactus, Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/99.
   Now generalized with the more general setup found in Mathematica
   http://documents.wolfram.com/v5/Built-inFunctions/NumericalComputation/
     EquationSolving/AdvancedDocumentation/NDSolve.html
   See also the description for Matlab, 
   "Behind and Beyond the MATLAB ODE Suite"
   okumedia.cc.osaka-kyoiku.ac.jp/~ashino/pdf/2651.pdf
   www.site.uottawa.ca/~remi/odesuite.ps

   Issues:
   - which of the "modern" Runge-Kutta features matter for method of lines?
     time step control? implicity?
   - currently not memory optimized, let's experiment first
*/

#include "sgrid.h"
#include "evolve.h"


/* pointers to temporary variable lists */
#define NSTAGES 6
tVarList *u_k[NSTAGES];




/* various Butcher tables that provide the coefficients for Runge-Kutta */
typedef struct {
  int nstages;
  double a[NSTAGES][NSTAGES];
  double b[NSTAGES];
  double c[NSTAGES];
} tButcher;

/* RK1: Euler */
tButcher rk1 = {
  1,
  {{0}},
  {1}, 
  {1}
};

/* RK2: modified Euler (?) */
tButcher rk2 = {
  2,
  {{0}, {1}},
  {0.5, 0.5}, 
  {1, 1}
};

/* RK2: midpoint method */
tButcher rk2a = {
  2,
  {{0}, {0.5}},
  {0, 1  }, 
  {0, 0.5}
};


/* RK3: classic Heun, version 1, Mathematica's 3(2) */
tButcher rk3 = {
  3,
  {{0}, {0.5}, {-1, 2}},
  {1./6, 2./3, 1./6}, 
  {0,    1./2, 1   }
};

/* RK3: classic Heun, version 2 */
tButcher rk3a = {
  3,
  {{0}, {1./3}, {0, 2./3}},
  {1./4, 0   , 3./4}, 
  {0,    1./3, 2./3}
};

/* RK3: Matlab's ode23, Ralston */
tButcher rk3b = {
  3,
  {{0}, {1./2}, {0, 3./4}},
  {2./9, 1./3, 4./9}, 
  {0,    1./2, 3./4}
};


/* RK4: classic Runge-Kutta */
tButcher rk4 = {
  4,
  {{0}, {0.5}, {0, 0.5}, {0, 0, 1}},
  {1./6, 1./3, 1./3, 1./6}, 
  {0,    1./2, 1./2, 1   }
};

/* RK4: Mathematica's 4(3) */
tButcher rk4a = {
  4,
  {{0}, {2./5}, {-3./20, 3./4}, {19./44, -15./44, 10./11}},
  {11./72, 25./72, 25./72, 11./72}, 
  {2./5, 3./5, 1, 1}
};


/* RK5: Dormand-Prince, ode45 in Matlab */
tButcher rk5 = {
  6,
  {{0}, 
   {1./5}, 
   {3./40, 9./40}, 
   {44./45, -56./15, 32./9},
   {19372./6561, -25360./2187, 64448./6561, -212./729},
   {9017./3168, -355./33, 46732./5247, 49./176, -5103./18656}
  },
  {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84}, 
  {1./5, 3./10, 4./5, 8./9, 1, 1}
};





/* for debugging: print Butcher table */
void prbutchertable(tButcher *b)
{
  int i, j;

  printf("Runge-Kutta: Butcher table for %s\n", Gets("evolution_method_rk")); 
  for (i = 0; i < b->nstages; i++) {
    printf("%6.3f  ", b->c[i]);
    for (j = 0; j < b->nstages; j++) 
      printf("%6.3f ", b->a[i][j]);
    printf("\n");
  }
  printf("        ");
  for (j = 0; j < b->nstages; j++)
    printf("%6.3f ", b->b[j]);
  printf("\n");
}





/* generic explicit Runge-Kutta */
void evolve_rk_generic(tGrid *grid, tButcher *rk) 
{
  double dt = grid->dt;
  int i, j;

  /* initialize: u_p = u_c */
  vlcopy(u_p, u_c);

  /* iterate over stages */
  for (i = 0; i < rk->nstages; i++) {
    if (i) vlcopy(u_c, u_p);

    /* u_c = u_p + Sum aij * kj
       note that vladd does ignores terms with zero coefficient */
    for (j = 0; j < i; j++) 
      vladd(u_c, 1.0, u_c, dt*rk->a[i][j], u_k[j]); 

    /* set algebraic BCs for u_c, and any other stuff needed */
    if(evolve_algebraicConditions) evolve_algebraicConditions(u_c, u_p);

    /* u_ki = F(u_c) 
       step size 0.0 means right-hand side is returned in u_ki
       u_p is not needed but passed in as dummy
    */
    evolve_rhs(u_k[i], u_p, 0.0, u_c);

    /* synchronize */
    //sgridpi_vlsynchronize(u_k[i]);
  }

  /* compute final time level, u_c = u_p + Sum bj * kj */
  vlcopy(u_c, u_p);
  for (j = 0; j < rk->nstages; j++)
    vladd(u_c, 1.0, u_c, dt*rk->b[j], u_k[j]);

  /* set algebraic BCs for u_c, and any other stuff needed */
  if(evolve_algebraicConditions) evolve_algebraicConditions(u_c, u_p);
}




/* explicit Runge-Kutta */
void evolve_rk(tGrid *grid) 
{
  static int firstcall = 1;
  tButcher *rk;
  char s[3] = "_r";
  int i;

  /* pick method */
  if      (Getv("evolution_method_rk", "rk5"))      rk = &rk5;
  else if (Getv("evolution_method_rk", "rk4a"))     rk = &rk4a;
  else if (Getv("evolution_method_rk", "rk4"))      rk = &rk4;
  else if (Getv("evolution_method_rk", "rk3b"))     rk = &rk3b;
  else if (Getv("evolution_method_rk", "rk3a"))     rk = &rk3a;
  else if (Getv("evolution_method_rk", "rk3"))      rk = &rk3;
  else if (Getv("evolution_method_rk", "rk2a"))     rk = &rk2a;
  else if (Getv("evolution_method_rk", "rk2"))      rk = &rk2;
  else if (Getv("evolution_method_rk", "rk1"))      rk = &rk1;
  else if (Getv("evolution_method_rk", "midpoint")) rk = &rk2a;
  else if (Getv("evolution_method_rk", "euler"))    rk = &rk1;
  else
    errorexits("unknown evolution_method_rk %s", s);

  /* add auxiliary fields to data base, no storage yet 
     incrementing second letter in s gives _r, _s, _t, ...
  */
  if (firstcall) {
    firstcall = 0;
    for (i = 0; i < rk->nstages; i++, s[1]++) 
      u_k[i] = AddDuplicate(u_c, s);

    if (0) prbutchertable(rk);
  }

  // strange bug? rk4 shadows rk4a
  // printf("%s, getv %d\n", Gets("evolution_method_rk"), Getv("evolution_method_rk", "rk3"));

  /* store current grid and time=1 in u_k[i] */
  for (i = 0; i < rk->nstages; i++)
    if (u_k[i]) { u_k[i]->grid = grid;  u_k[i]->time = 1.0; }

  /* turn on memory (if varlist is non null) */
  for (i = 0; i < rk->nstages; i++)
    enablevarlist(u_k[i]);

  /* evolve */
  evolve_rk_generic(grid, rk);

  /* turn off memory */
  for (i = 0; i < rk->nstages; i++)
    disablevarlist(u_k[i]);
}
