/* ADMvars.c */
/* Bernd Bruegmann 11/02 */

#include "sgrid.h"
#include "ADMvars.h"




/* compute ADM constraints at ANALYZE time */
int computeADMconstraints(tGrid *grid)
{
  int NANcheck = Getv("ADMvars_NANcheck", "yes");
  tVarList *vl, *wl;
  
  vl = vlalloc(grid);
  vlpush(vl, Ind("ham"));
  vlpush(vl, Ind("momx"));
  vlpush(vl, Ind("trK"));
  /* add normalzed constraints to list if ADMvars_normalizedConstraints=yes */
  if(Getv("ADMvars_normalizedConstraints", "yes"))
  {
    vlpush(vl, Ind("normham"));
    vlpush(vl, Ind("normmomx"));
  }
  else /* if ADMvars_normalizedConstraints=no add constraints again to list */
  {
    vlpush(vl, Ind("ham"));
    vlpush(vl, Ind("momx"));
  }

  if (timeforoutput(grid, vl)) {

    /*    printf("computing ADM constraints\n"); */

    enablevarlist(vl);
    wl = vlalloc(grid);
    vlpush(wl, Ind("gxx"));
    vlpush(wl, Ind("Kxx"));
    vlpushvl(wl, vl);

    /* compute the adm constraints */
    ADMconstraints(wl);

    /* set boundaries */
    //set_boundary_symmetry(level, vl);
    //set_boundary_periodic(level, vl);
    //set_boundary_refinement(level, vl);

    /* sunchronize */
    //bampi_vlsynchronize(vl);
    
    /* clean up */
    vlfree(wl);
    
    /* Check for NANs and INFs */
    if (NANcheck)
      if (CheckIfFinite(grid, "ham") != 0) {
	printf("time %f, iteration %d:\n",
	       grid->time, grid->iteration);
	errorexit("CheckIfFinite: NAN/INF error: ham is not finite!");
      }

    /* output */
    /* this unfinished, but we want something like:
    write_varlist(grid, vl);
    disablevarlist(vl);
    */
  } 
    
  vlfree(vl);
  return 0;
}


int allocateADMvars(tGrid *grid)
{
  
  /* enable g, K and their derivs */
  if(Getv("ADMvars_memory_for_g_K", "yes"))
  {
    enablevar(grid, Ind("gxx"));
    enablevar(grid, Ind("Kxx"));
  }
  
  /* enable psiandderivs */
  if(Getv("ADMvars_memory_for_psiandderivs", "yes"))
  {
    psiandderivs->grid = grid;
    enablevarlist(psiandderivs);
  }

  /* enable derivs of g and K */
  if(Getv("ADMvars_memory_for_dg_ddg_dK", "yes"))
  {
    enablevar(grid, Ind("ADMvars_dgxxx"));
    enablevar(grid, Ind("ADMvars_ddgxxxx"));
    enablevar(grid, Ind("ADMvars_dKxxx"));
  }

  return 0;
}
