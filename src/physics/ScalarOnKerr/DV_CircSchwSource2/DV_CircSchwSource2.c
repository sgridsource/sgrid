/* DV_CircSchwSource2.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "DV_CircSchwSource2.h"


/* initialize DV_CircSchwSource2 */
int DV_CircSchwSource2_startup(tGrid *grid)
{
  double x,y;
  double M=1.0;
  int i,j;
  printf("Initializing DV_CircSchwSource2:\n");

  DV_set_parameters(M, 1.0, 10.0*M);
  DV_show_parameters();

  return 0;
}
