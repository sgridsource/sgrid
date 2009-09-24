/* sgrid_TestID.c */
/* Wolfgang Tichy, 09/2009 */

#include "sgrid.h"
#include "TestID.h"


int sgrid_TestID() 
{
  if (!Getv("physics", "TestID")) return 0;
  printf("Adding TestID\n");

  /* functions: */
  AddFun(INITIALDATA, TestID, "TestID: set initial data for a single BH");

  /* Parameters: */
  AddPar("TestID_type", "1dGaugeWave",
         "type of ini. data [1dGaugeWave,1dLinearWave,1dRandomNoise]");
  AddPar("TestID_amplitude", "1", "amplitude A of test data");
  AddPar("TestID_initial_lapse", "donothing",
         "initial lapse [donothing,one]");
  AddPar("TestID_initial_shift", "donothing",
         "initial shift [donothing,zero]");

  return 0;
}
