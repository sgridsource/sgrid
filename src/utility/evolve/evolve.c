/* evolve.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/02 */

/* provide the glue for the different evolution schemes */

#include "sgrid.h"
#include "evolve.h"



/* variables and their copies for evolution 
   initialized to 0 for clarity since we rely on that
*/
tVarList *u_c = 0, *u_p = 0, *u_q = 0, *u_r = 0;
tVarList *u_pp = 0;


/* pointer to evolution routine */
void (*evolve_rhs)
  (tVarList *unew, tVarList *upre, double c, tVarList *ucur) = 0;

/* Pointer to routine that can set algebraic Boundary Conditions
   after each sub time step (e.g. in rk).
   Note that any time deriv BCs such as d_t u = something, should be
   set in evolve_rhs itself. Here we set only stuff like
   unew = upre, or u1new = -u2new  <-- at new time */
void (*evolve_algebraicConditions)
  (tVarList *unew, tVarList *upre) = 0;



/* store variables locally */ 
void evolve_vlregister(tVarList *u) 
{
  u_c = u;
}




/* retrieve variables */
void evolve_vlretrieve(tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_pp)
{
  if (vlu_c)  *vlu_c  = u_c;
  if (vlu_p)  *vlu_p  = u_p;
  if (vlu_pp) *vlu_pp = u_pp;
}




/* store evolution routine locally */
void evolve_rhsregister(void (*f)(tVarList *, tVarList *, double, tVarList *)) 
{
  evolve_rhs = f;
}

/* set pointer to evolve_algebraicConditions */
void evolve_algebraicConditionsregister(void (*f)(tVarList *, tVarList *))
{
  evolve_algebraicConditions = f;
}


/* evolve wrapper for function skeleton */
int evolve(tGrid *grid) 
{
  int i, j;

  /* make sure we have an evolution routine and a list of variables */
  if (!evolve_rhs) {
    return 0;
    // let's assume someone else provides evolution
    // errorexit("evolve: no evolution routine");
  }
  if (!u_c) 
    errorexit("evolve: no list of variables");

  /* store grid for this call to evolve */
  u_c->grid = grid;

  /* add new variables before first time step */
  if (grid->iteration == 0) {
    printf("Adding variables for evolution:\n");

    /* add additional variables */
    if (Getv("evolution_method", "euler")) {
      u_p = AddDuplicate(u_c, "_p");
    }
    if (Getv("evolution_method", "icn")) {
      u_p = AddDuplicate(u_c, "_p");
      u_q = AddDuplicate(u_c, "_q");
    }
    if (Getv("evolution_method", "rk")) {
      u_p = AddDuplicate(u_c, "_p");
    }
    //if (Geti("amr_lmax") > 0) {
    //  if (!u_p) AddDuplicate(u_c, "_p");
    //  u_pp = AddDuplicate(u_c, "_pp");
    //}
  }

  /* store current grid in existing variable lists */
  if (u_p) u_p->grid = grid;
  if (u_q) u_q->grid = grid;
  if (u_r) u_r->grid = grid;
  // if (u_pp) u_pp->grid = grid;

  /* turn on memory (if varlist is non null) */
  enablevarlist(u_p);
  enablevarlist(u_q);
  enablevarlist(u_r);
  // enablevarlist(u_pp);

  /* choose between different method of line integrators */
  if (Getv("evolution_method", "icn"))
    evolve_icn(grid);
  else if (Getv("evolution_method", "euler"))
    evolve_euler(grid);
  else if (Getv("evolution_method", "rk"))
    evolve_rk(grid);
  else
    errorexit("evolution_method not known, use icn");

  /* compute difference between new and old data, store in u_change */
  for (j = 0; j < u_c->n; j++)
  {
    if (Getv("evolve_compute_change", VarName(u_c->index[j])))
    {
      tVarList *var, *change;
      int b;

      /* make temporary varlist containing only the current variable */
      var = vlalloc(grid);
      vlpush( var, IndComponent0(u_c->index[j]) );
      enablevarlist(var);

      /* duplicate var, to create var_change */
      change = AddDuplicateEnable(var, "_change");

      for (b = 0; b < grid->nboxes; b++)
      {
        tBox *box = grid->box[b];
        double *old = box->v[u_p->index[j]];
        double *new = box->v[u_c->index[j]];
        double *ch = box->v[change->index[VarComponent(u_c->index[j])]];

        /* loop over all points and write into var_change */
        forallpoints(box, i)
        {
          ch[i] = new[i] - old[i];
        }
        /* free temp. var lists */
        vlfree(change);
        vlfree(var);
      }
    }
  }

  /* turn off memory (if varlist is non null) */
  if (Getv("evolve_persist", "no")) {
    disablevarlist(u_p);
    disablevarlist(u_q);
    disablevarlist(u_r);
  }

  return 0;
}
