/* sgrid_DV_CircSchwSource.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "DV_CircSchwSource.h"




int sgrid_DV_CircSchwSource(void) 
{
  if (!Getv("physics", "ScalarOnKerr")) return 0;
  printf("Adding DV_CircSchwSource\n");

  /* functions */
  AddFun(POST_INITIALDATA, DV_CircSchwSource_startup, 
	 "initialize DV_CircSchwSource system from adm initial data");

  /* variables */
  //AddVar("DV_CircSchwSource_psi", "", "scalar");

  /* parameters */
  AddPar("DV_CircSchwSource_useWindow", "no", "if we use Windowfunc. [no,yes]");
  AddPar("DV_CircSchwSource_DVWindow_n",     "8", "ian_n in Windowfunc.");
  AddPar("DV_CircSchwSource_DVWindow_width", "2", "ian_width in Windowfunc.");
         
  return 0;
}
