/* checkpoint.c */
/* Wolfgang Tichy 11/2008 */

/* write and read checkpoint files 
*/

#include "sgrid.h"
#include "checkpoint.h"



/* fake the mpi calls this module makes */
int sgrid_mpi_size(void){ return 1; }
int sgrid_mpi_rank(void){ return 0; }
int sgrid_mpi_allreduce_sum_int(int n){ return n; }
void sgrid_mpi_barrier(void){;}
    

/* main entry into checkpointing */
int checkpoint(tGrid *grid)
{
  int di;
  double dt, dT;
  int restart_it;
  int pr = 0;

  /* check whether we are on */
  if(Getv("checkpoint", "no")) return 0;
  if(pr) printf("checkpoint called at top grid iteration %d\n",
		 grid->iteration);

  /* iteration at which we will do a restart if there are checkpoint files */
  restart_it = Geti("checkpoint_restart_it");

  /* if we are in the restart phase */
  if (Getv("checkpoint", "restart"))
  {
    /* here we could do restarts from initial data */ 
    if (grid->iteration == 0)
      ;
    
    /* clobber data at restart_it, usually after one full evolution step */
    if (grid->iteration == restart_it || restart_it<0)
    {
      printf("checkpoint restart: read files\n");
      checkpoint_read(grid);
      checkpoint_copy_output();
      Sets("checkpoint", "yes");

      /* run for that much longer */
      dT = Getd("checkpoint_DeltaT");
      if (dT > 0) {
	Setd("finaltime", grid->time + dT);
	Seti("iterations", ((int) ((grid->time + dT)/grid->dt + 0.5)));
	printf("new time is %.3f, adding %.3f to the clock, run until %.3f\n",
	       grid->time, dT, grid->time + dT);
	printf("embarking on iterations %d to %d\n",
	       grid->iteration + 1, Geti("iterations")); 
      }
      return 1; /* return 1 after a real restart */
    }

    /* return since we are still starting up and 
       don't want new checkpoint files            */
    return 0;
  }

  /* if we get here, we are not in the restart phase */

  /* determine how often we want to checkpoint */
  dt = Getd("checkpoint_dt");
  if (dt > 0)
    di = dt/grid->dt + dequaleps;
  else
    di = Geti("checkpoint_di");
  if(pr) printf("di = %d\n", di);
  
  /* check whether it is time to write a checkpoint */
  if (grid->iteration % di) {
    if(pr) printf("not time to write checkpoint, return\n");
    return 0;
  }

  /* we do not write checkpoints after initial data, if restart_it>0 */
  if(grid->iteration == 0 && restart_it>0)  return 0;

  /* write checkpoint */
  checkpoint_write(grid);

  return 0;
}




/* return pointer to filename, has to be freed by caller 
   this is actually quite tricky since machine dependent
   - assume for now that each processor can write its own output!
   - it may be more efficient to collect all data onto processor 0
     and then write one file
*/
char *checkpoint_filename(char *suffix)
{
  char *dir = cmalloc(strlen(Gets("outdir")) + strlen(suffix) + 1);
  char *filename = cmalloc(1000);
  char formatstring[100];

  /* name of directory */
  sprintf(dir, "%s%s", Gets("outdir"), suffix);
    
  /* variable length format for different number of processors
     (same as for stdout.01 etc.)
  */
  snprintf(formatstring, 100, "%%s/checkpoint.%%0%dd",
	   (int) log10(sgrid_mpi_size())+1);
  snprintf(filename, 1000, formatstring, dir, sgrid_mpi_rank());

  if (1) printf("  checkpoint filename = %s\n", filename);

  free(dir);
  return filename;
}




/* open checkpoint file for this processor 
   return 0 if not
   return file pointer to opened file
   this file has to be closed by caller
*/
FILE *checkpoint_openfiles(char *suffix)
{
  char *filename = checkpoint_filename(suffix);
  FILE *fp;
  int nfp;

  fp = fopen(filename, "rb");
  nfp = fp ? 1 : 0;
  nfp = sgrid_mpi_allreduce_sum_int(nfp);

  if (nfp == 0) {
    if (1) printf(
      "  there are no checkpoint files to read yet\n");
    return 0;
  }

  if (nfp < sgrid_mpi_size()) {
    if (fp) fclose(fp);
    errorexit("there are some, but not all of the expected checkpoint files");
  }

  free(filename);
  return fp;
}




/* check whether there is a checkpoint file for each processor */
int checkpoint_checkforfiles(char *suffix)
{
  FILE *fp = checkpoint_openfiles(suffix);

  if (fp) {
    fclose(fp);
    if (1) printf("  found checkpoint files, going into restart mode\n");
    return 1;
  }

  if (1) printf("  no checkpoint files, no restart\n");
  return 0;
}

/* check whether we have a complete set of checkpoint files in
   the "outdir_previous" directory */
int checkpoint_checkifrestart(tGrid* grid)
{
  printf("checkpoint_checkifrestart: checking if we can restart...\n");
  /* special for checkpointing:
     check whether we have a complete set of checkpoint files in
     the "outdir_previous" directory */
  if (checkpoint_checkforfiles("_previous")) 
  {
    /* go into restart mode */
    Sets("checkpoint", "restart");
  }
  return 0;
}


/* copy output from previous run into current location so that all
   additional output can continue/append where we left off */
void checkpoint_copy_output()
{
  char *suffix = "_previous";
  char *current = Gets("outdir");
  char *previous = cmalloc(strlen(current) + strlen(suffix) + 1);

  /* name of previous directory */
  sprintf(previous, "%s%s", current, suffix);

  /* the regular output is all done through processor 0 */
  if (processor0) {
    
    /* remove current output directory
       it only contains iteration 0 and perhaps 1 of the restarted run 
    */
    system2("rm -rf", current);

    /* move or copy previous output */
    if( Getv("checkpoint_previous", "no") &&
       !GetvLax("sgrid_options", "--keep_previous") )
      system3("mv", previous, current);
    else
      system3("cp -pr", previous, current);
  }
  /* wait here to make sure proc0 has copied all there is */
  sgrid_mpi_barrier();
         
  /* redirect stdout and stderr for multiprocessor jobs */
  if (!processor0)
  {
    char s[1000], f[100];

    sprintf(f, "%%s/stdout.%%0%dd", (int) log10(sgrid_mpi_size())+1);
    sprintf(s, f, current, sgrid_mpi_rank());
    freopen(s, "a", stdout);
    freopen(s, "a", stderr);
  }
  
  /* clean up */
  free(previous);
}




/* create list of variables that are really needed for checkpoint
   caller has to free storage
   - needs to be tested, made more flexible! */
tVarList *checkpoint_varlist(tGrid *grid) 
{
  tVarList *u, *up, *upp;
  tVarList *vl = vlalloc(grid);
  int i;
  
  /* get all variables that are maintained for evolution
     note that the ADM variables are derived quantities e.g. for BSSN */
  evolve_vlretrieve(&u, &up, &upp);
  vlpushvl(vl, u);
  vlpushvl(vl, up);
  vlpushvl(vl, upp);

  /* add all variables with storage (in box0) contained in
     checkpoint_additional_variables */
  for (i = 0; i < grid->nvariables; i++)
  {
    tBox *box=grid->box[0];
    if (box->v[i] != 0)
      if( Getv("checkpoint_additional_variables", VarName(i)) )
      {
        vldrop(vl, i); /* <--drop it in case it is there already */
        vlpush(vl, i); /* <--now add it again */
      }
  }

  /* there are other variables that we may need/want
     - x, y, z: available from startup
     - psi and derivs: should be set during initial data (!)
     - K_initial: ?
     some of these could be added to the par checkpoint_additional_variables
     by other modules...
  */

  /* print variables that are considered to be essential */
  if (0) prvarlist(vl);

  return vl;
}
