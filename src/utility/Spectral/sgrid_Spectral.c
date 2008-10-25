/* sgrid_Spectral.c */
/* Wolfgang Tichy 2003 */

#include "sgrid.h"
#include "Spectral.h"


int sgrid_Spectral(void) 
{
  int b;
  printf("Adding Spectral\n");

  AddPar("Spectral_second_deriv_order", "123",
         "order in which mixed second derivs are taken "
         "[123,132,213,231,312,321]");

  for(b=0; b<Geti("nboxes"); b++)
  {
    char str[1000];

    /* whether we use FFTs and when we use them */
    snprintf(str, 999, "box%d_TransformType1", b);
    AddPar(str, "512 NUMREC_FFT", "FFT-type in X-direction "
                "[(min # of points to try to use FFT) (FFT-type)]");
    snprintf(str, 999, "box%d_TransformType2", b);
    AddPar(str, "512 NUMREC_FFT", "FFT-type in Y-direction " 
                "[(min # of points to try to use FFT) (FFT-type)]");
    snprintf(str, 999, "box%d_TransformType3", b);
    AddPar(str, "512 NUMREC_FFT", "FFT-type in Z-direction "
                "[(min # of points to try to use FFT) (FFT-type)]");
  }

  return 0;
}
