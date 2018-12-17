/* bam_checkpoint.h */
/* Wolfgang Tichy 11/2008 */

#include "wolfio.h"

/* checkpoint.c */
char *checkpoint_filename(char *suffix, char *flag);
int checkpoint(tGrid *grid);
int parameterio_write_current_pars(tGrid *g);
int parameterio_update_pars(tGrid *g);
void checkpoint_interpolate_Vars_get_Pars(char *filename, 
                                          tVarList *varlist, char *parlist);
