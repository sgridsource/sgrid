/* sgrid_s2kit.c */
/* Wolfgang Tich 8/2007 */

#include "sgrid.h"
#include "s2kit.h"


int sgrid_s2kit() 
{
  /* if (!Getv("physics", "s2kit")) return; */
  printf("Adding s2kit\n");

  /* functions */

  /* variables */

  /* parameters */

  /* test */
  if (Getv("physics", "s2kit_test"))
  {
    AddFun(INITIALDATA, s2kit_test_Naive_YlmFilter, "do test");
    //AddVar("temp1", "", "temp var");
    AddVar("s2kit_test_var1", "", "test var1");
    AddVar("s2kit_test_var2", "", "test var2");
    AddVar("s2kit_test_fil1", "", "after filtering");
    AddVar("s2kit_test_fil2", "", "after filtering");
    AddVar("s2kit_test_dif1", "", "difference after filtering");
    AddVar("s2kit_test_dif2", "", "difference after filtering");
    // AddPar("s2kit_method_order", "0", "expected order of convergence");
  }  
  return 0;
}
