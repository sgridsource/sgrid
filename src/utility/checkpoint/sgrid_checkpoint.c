/* sgrid_checkpoint.c */
/* Wolfgang Tichy 11/2008 */

#include "sgrid.h"
#include "checkpoint.h"



int sgrid_checkpoint(void)
{
  printf("Adding checkpoint\n");
  AddPar("checkpoint", "no", "whether to checkpoint [no,yes]");
  if (Getv("checkpoint", "no")) return 0;

  AddPar("checkpoint_di", "1", "how often to checkpoint");
  AddPar("checkpoint_dt", "0", "how often to checkpoint");

  AddPar("checkpoint_DeltaT", "0",
	 "run for that much beyond previous checkpoint");

  AddPar("checkpoint_variables", "auto",
	 "which variables to checkpoint [auto,all]");
  /* additional variables which also need to be saved,
     besides the ones implied by auto */
  AddPar("checkpoint_additional_variables", "",
	 "additional variables we need to save [any varnames]");
  AddPar("checkpoint_previous", "no",
	 "whether to keep previous directory [no]");

  /* special for checkpointing:
     check whether we have a complete set of checkpoint files in 
     the "outdir_previous" directory */
  if (checkpoint_checkforfiles("_previous")) {
    /* go into restart mode */
    Sets("checkpoint", "restart");
  }
  return 0;
}
