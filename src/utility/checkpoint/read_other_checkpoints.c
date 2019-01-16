/* read_other_checkpoints.c */
/* Wolfgang Tichy 11/2011 */

/* read checkpoint files from a different grid */

#include "sgrid.h"
#include "checkpoint.h"



/* global vars from parameter data base */
extern tParameter *pdb;
extern int npdb, npdbmax;


/* swap two tParameter pointers
   If we have: tParameter *pdb1; tParameter *pdb2;
   call like:  swap_pdb_pointers(&pdb1, &pdb2); */
void swap_tParameter_pointers(tParameter **p1, tParameter **p2)
{
  tParameter *pdb0;

  /* swap pointers p1 and p2 */
   pdb0 = *p1;
   *p1  = *p2;  
   *p2  = pdb0;
}

/* Read a list of vars and pars from a specific checkpoint file.
   The Vars are interpolated, in case the grids are different. */
void checkpoint_interpolate_Vars_get_Pars(char *filename, 
                                          tVarList *varlist, char *parlist1)
{
  tGrid *grid = varlist->grid;
  tGrid *grid2;
  tParameter *pdb2;
  char *par;
  char *val;
  int tempind = Ind("temp1");
  char *parlist;
  int b, j;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

//printVarList(varlist);         

  /* read the file into grid2, pdb2 */
  printf("Setting parameters and variables in pdb2 and grid2 by "
         "reading from file:\n%s\n", filename);
  read_checkpoint_into_grid2_pdb2(grid, filename, &grid2, &pdb2);

  /* we can only interpolate if the box boundaries are in the same place
     => read them from pdb2 as well */
  /* make list of pars that need to be taken from pdb2: */
  /* add parlist1 to parlist */
  parlist = cmalloc(strlen(parlist1) + (grid->nboxes)*8*14);
  strcpy(parlist, parlist1);

  /* add box?_min?, box?_max?, box?_CI_s, box?_CI_xc to parlist */
  forallboxes(grid, b)
  {
    char str[1000];
    snprintf(str, 999, " box%d_min1", b);
    strncat(parlist, str, 999);
    snprintf(str, 999, " box%d_max1", b);
    strncat(parlist, str, 999);

    snprintf(str, 999, " box%d_min2", b);
    strncat(parlist, str, 999);
    snprintf(str, 999, " box%d_max2", b);
    strncat(parlist, str, 999);

    snprintf(str, 999, " box%d_min3", b);
    strncat(parlist, str, 999);
    snprintf(str, 999, " box%d_max3", b);
    strncat(parlist, str, 999);

    snprintf(str, 999, " box%d_CI_s", b);
    strncat(parlist, str, 999);
    snprintf(str, 999, " box%d_CI_xc", b);
    strncat(parlist, str, 999);
  }

  /* loop over pars in parlist */
  printf("Setting the following parameters in pdb to values from pdb2:\n");
  NextEntry(0);
  while( (par=NextEntry(parlist)) != NULL)
  {
    /* swap parameter databases, so that we can use e.g. Gets on pdb2 */
    swap_tParameter_pointers(&pdb, &pdb2);
    val = Gets(par); /* get par value from pdb2 */
    
    swap_tParameter_pointers(&pdb, &pdb2);
    printf(" %s = %s\n", par, val);
    Sets(par, val);  /* set par value in pdb */
    fflush(stdout);
  }
  printf(" Note: Here we added box?_min?, box?_max?, box?_CI_s, box?_CI_xc "
         " because we\n"
         " can only interpolate if the box boundaries are in the same place.\n");

  /* we can only interpolate if the box boundaries are in the same place
     on both grids => Adjust grid1. */
  set_BoxStructures_fromPars(grid, 0);

  /* reset x,y,z, dXdx and such */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* interpolate varlist from grid2 to grid */
  printf("Interpolating from grid2 to grid:\n");
  for(j = 0; j < varlist->n; j++)
  {
    printf(" %s ... ", VarName(varlist->index[j]));
    fflush(stdout);
    spec_interpolate_Var_from_grid2_to_grid1(grid, grid2, 
                                             varlist->index[j], tempind);
    printf("done.\n");
    fflush(stdout);
  }

  /* free grid2 and pdb2, which were read from file */
  free_grid(grid2);
  free_pdb(pdb2,npdb);
  free(parlist);
}

/* read checkpoint file into grid2 and pdb2 */
/* NOTE: grid2 and pdb2 have to be freed by the caller.
   If we have: tGrid *grid2;  tParameter *pdb2; 
   call like:  read_checkpoint_into_grid2_pdb2(grid, filename, &grid2, &pdb2); */
void read_checkpoint_into_grid2_pdb2(tGrid *grid, char *filename,
                                     tGrid **grid2, tParameter **pdb2)
{
  tGrid *g2;
  tGrid *g_dummy;
  tParameter *p2;
  FILE *fp;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  fp = fopen(filename, "rb");
  if(fp==NULL)
  {
    printf("read_checkpoint_into_grid2_pdb2: cannot open %s\n", filename);
    errorexit("read_checkpoint_into_grid2_pdb2: checkpoint file does not exist!");
  }

  /* Create a copy of parameter database pdb in p2.
     Note: pdb is a global var. */
  create_copy_of_pdb1_in_pdb2(pdb,npdb,npdbmax, &p2);

  /* swap pointers pdb and p2 */
  swap_tParameter_pointers(&pdb, &p2);

  /* Create an empty dummy grid */
  g_dummy = make_empty_grid(globalnvariables, 0);
  
  /* get new pars And Iterations */
  checkpoint_read_ParsAndIterations_local(g_dummy, fp);
 
  /* make grid g2 */
  g2 = make_grid(0);

  /* set g2->time, g2->iteration */
  g2->time      = g_dummy->time;
  g2->iteration = g_dummy->iteration;

  /* free dummy grid */
  free_grid(g_dummy);

  /* enable all vars of grid also on grid g2 */
  enablesamevars(grid, g2);

  /* initialize coords on g2 */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(g2);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* read Vars from checkpoint into g2 */
  checkpoint_read_Vars_local(g2, fp);

  /* re-swap pointers pdb and p2 */
  swap_tParameter_pointers(&pdb, &p2);

  /* set grid2 and pdb2 */
  *grid2 = g2;
  *pdb2  = p2;

  /* clean up */
  fclose(fp);
}
