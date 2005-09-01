/* ADMvars.c */
/* Bernd Bruegmann 11/02 */

#include "sgrid.h"
#include "ADMvars.h"




/* compute ADM constraints at ANALYZE time */
int computeADMconstraints(tGrid *grid)
{
  int NANcheck = Getv("ADMvars_NANcheck", "yes");
  tVarList *vl, *wl;

  /* compute the constraints */
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

  if (timeforoutput(grid, vl))
  {
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

  /* compute ADM energy */
  vl = vlalloc(grid);
  vlpush(vl, Ind("E_ADM"));

  if (timeforoutput(grid, vl))
  {
    int bi;

    enablevarlist(vl);
    wl = vlalloc(grid);
    vlpushvl(wl, vl);
    vlpush(wl, Ind("gxx"));
    vlpush(wl, Ind("x"));
    vlpush(wl, Ind("y"));
    vlpush(wl, Ind("z"));
    vlpush(wl, Ind("psi"));
    vlpush(wl, Ind("dpsiopsix"));

    /* get the integrand for the ADM energy in var E_ADM */
    ADMenergy_spheric_intergrand(wl);

    /* now integrate integrand in var E_ADM to get the real E_ADM */
    forallboxes(grid, bi)
    {
      tBox *box = grid->box[bi];
      char str[1000];
  
      snprintf(str, 999, "box%d_Coordinates", bi);
      if( Getv(str, "SphericalDF") )
      {
        double *E_ADM = box->v[Ind("E_ADM")];
        spec_sphericalDF2dIntegral(box, E_ADM, E_ADM);
      }
      else
        errorexits("ADMvars.c: I don't know how to do surface integrals in "
                   "%s coordinates", str);
    }
    vlfree(wl);
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

  /* enable temp1, temp2, temp3 */
  if(Getv("ADMvars_memory_for_temp123", "yes"))
  {
    enablevar(grid, Ind("temp1"));
    enablevar(grid, Ind("temp2"));
    enablevar(grid, Ind("temp3"));
  }

  return 0;
}
