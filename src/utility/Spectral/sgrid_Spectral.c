/* sgrid_Spectral.c */
/* Wolfgang Tichy 2003 */

#include "sgrid.h"
#include "Spectral.h"


int sgrid_Spectral(void) 
{
  printf("Adding Spectral\n");

  AddPar("Spectral_second_deriv_order", "123",
         "order in which mixed second derivs are taken "
         "[123,132,213,231,312,321]");

  return 0;
}
