/* sgrid_Spectral.c */
/* Wolfgang Tichy 2003 */

#include "sgrid.h"
#include "Spectral.h"


int sgrid_Spectral(void) 
{
  int b;
  printf("Adding Spectral\n");

  /* functions */
  AddFun(PRE_GRID, init_FFTW3_plans, 
         "initialize global FFTW3 plans in FFTs_for_sgrid.c");
  AddFun(POST_FINALIZE_GRID, free_FFTW3_plans, 
         "free the global FFTW3 plans in FFTs_for_sgrid.c");

  /* parameters */
  AddPar("FFTW3_planner_flags", "FFTW_ESTIMATE FFTW_UNALIGNED",
         "FFTW planner flags we use. Possible values are in the pars FFTW_*"
         "(see function add_special_FFTW3_pars)");
  add_special_FFTW3_pars();
  AddPar("Spectral_second_deriv_order", "123",
         "order in which mixed second derivs are taken "
         "[123,132,213,231,312,321]");

  for(b=0; b<Geti("nboxes"); b++)
  {
    char str[1000];

    /* whether we use FFTs and when we use them */
    /* default: use NUMREC_FFT if number of points N>=64 and N=2^n */
    /* NOTE: if FFTW3 is compiled in we should really use:
        box0_TransformType1 = 25 FFTW3  # if basis1=ChebExtrema
        box0_TransformType2 = 24 FFTW3  # if basis2=Fourier
        box0_TransformType3 = 24 FFTW3  # if basis3=Fourier   */
    snprintf(str, 999, "box%d_TransformType1", b);
    AddPar(str, "64 NUMREC_FFT", "FFT-type in X-direction "
                "[(min # of points to try to use FFT) (FFT-type)]");
    snprintf(str, 999, "box%d_TransformType2", b);
    AddPar(str, "64 NUMREC_FFT", "FFT-type in Y-direction " 
                "[(min # of points to try to use FFT) (FFT-type)]");
    snprintf(str, 999, "box%d_TransformType3", b);
    AddPar(str, "64 NUMREC_FFT", "FFT-type in Z-direction "
                "[(min # of points to try to use FFT) (FFT-type)]");
  }

  return 0;
}
