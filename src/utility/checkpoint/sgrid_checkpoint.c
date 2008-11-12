/* sgrid_checkpoint.c */
/* Wolfgang Tichy 11/2008 */

#include "sgrid.h"
#include "checkpoint.h"



int sgrid_checkpoint(void)
{
  printf("Adding checkpoint\n");

  /* parameters */           
  AddPar("checkpoint", "no", "whether to checkpoint [no,yes]");
  if (Getv("checkpoint", "no")) return 0;

  /* functions */
  AddFun(PRE_GRID, checkpoint_checkifrestart, "check if we can restart");

  /* rest of parameters */
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

  return 0;
}
