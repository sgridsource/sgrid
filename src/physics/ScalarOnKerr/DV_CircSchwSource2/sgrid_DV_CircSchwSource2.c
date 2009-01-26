/* sgrid_DV_CircSchwSource2.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "DV_CircSchwSource2.h"




int sgrid_DV_CircSchwSource2(void) 
{
  if (!Getv("physics", "ScalarOnKerr")) return 0;
  printf("Adding DV_CircSchwSource2\n");

  /* functions */
  AddFun(POST_INITIALDATA, DV_CircSchwSource2_startup, 
	 "initialize DV_CircSchwSource2");

  /* variables */
  //AddVar("DV_CircSchwSource_psi", "", "scalar");

  /* parameters */
  AddPar("DV_CircSchwSource_useWindow", "yes", "if we use Windowfunc. [no,yes]");
  AddPar("DV_CircSchwSource_Window_type", "orig",
         "window type we use [orig,wolf,no]");
  AddPar("DV_CircSchwSource_DVWindow_n",     "8", "ian_n in set_OrigWindow");
  AddPar("DV_CircSchwSource_DVWindow_width", "2", "ian_width set_OrigWindow");
  AddPar("DV_CircSchwSource_Window_smooth", "0",
         "if smooth=3 source will be artificially smoothed to C^3");
  AddPar("DV_CircSchwSource_Window_q1", "1.2", "q1 determines the value of r where the window function, W=1/2");
  AddPar("DV_CircSchwSource_Window_s1", "1.9", "slope of W");
  AddPar("DV_CircSchwSource_Window_r1", "0", "r1 is the inner radius where W=0");
  AddPar("DV_CircSchwSource_Window_r2", "0", "r2 is the outer radius where W=1");
  AddPar("DV_CircSchwSource_Window_q3", "0", "q3 determines the value of r where the window function, W=1/2");
  AddPar("DV_CircSchwSource_Window_s3", "0", "slope of W");
  AddPar("DV_CircSchwSource_Window_r3", "0", "r3 is the inner radius where W=1");
  AddPar("DV_CircSchwSource_Window_r4", "0", "r4 is where W=0");
         
  return 0;
}
