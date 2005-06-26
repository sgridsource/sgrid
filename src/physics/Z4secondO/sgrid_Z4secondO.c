/* sgrid_Z4secondO.c */
/* Wolfgang Tichy  */

#include "sgrid.h"
#include "Z4secondO.h"


int sgrid_Z4secondO(void) 
{
  if (!Getv("physics", "Z4secondO")) return 0;
  printf("Adding Z4secondO\n");

  /* functions */
  AddFun(POST_INITIALDATA, Z4secondO_startup, 
	 "Initialize Z4secondO system from ADM initial data");

  /* variables */
  AddVar("Z4secondO_Theta", "",     "Z4 Theta");
  AddVar("Z4secondO_Z",     "i",    "Z4 Z_i");

  /* derivatives which need to be precomputed before each evo steo */
  /* Note: ADMvars_ddg, ADMvars_dg, ADMvars_dA are overwritten and used 
           in place of Z4secondO_ddg, Z4secondO_dg, Z4secondO_dK */
  AddVar("Z4secondO_dTheta", "i",    "1st deriv of Theta");
  AddVar("Z4secondO_dZ",     "ij",   "1st deriv of Z_i");
  AddVar("Z4secondO_dalp",   "i",    "1st deriv of Lapse");
  AddVar("Z4secondO_ddalp",  "(ij)", "2nd deriv of Lapse");
  AddVar("Z4secondO_dbeta",  "Ij",   "1st deriv of shift");
  AddVar("Z4secondO_ddbeta", "I(jk)","2nd deriv of shift");
    
  /* parameters */
  AddPar("Z4secondO_useDD", "yes",
         "wether we use the DD ops to compute second derivs [no,yes]");
  
  AddPar("Z4secondO_RtoRminusHfactor", "0",
         "how much of Hamiltonian H we subtract from Ricci scalar R");
  AddPar("Z4secondO_kappa1", "0", "exponential damping factor for Theta");
  AddPar("Z4secondO_kappa2", "0", "exponential damping factor for Z_i");
         
  AddPar("Z4secondO_lapsepsipower",  "0", "power of psi in lapse equation");
  AddPar("Z4secondO_lapseharmonicf", "2", "2 for 1+log, 1 for harmonic");
  AddPar("Z4secondO_lapseharmonicm", "0", "m in -m Theta term of lapse Eqn");
  AddPar("Z4secondO_subtractK0",     "no", 
	 "for harmonic and 1+log, subtract initial traceK or not");

  AddPar("Z4secondO_shiftpsipower",   "0", "power of 1/psi in shift equation");
  AddPar("Z4secondO_shiftalphapower", "0", "power of lapse in shift equation");
  AddPar("Z4secondO_shiftgammacoeff", "0.75", "coefficient in shift equation");
  AddPar("Z4secondO_shiftdriver",     "0", "coefficient of diffusion term");
         
  AddPar("Z4secondO_shift",  "constant",
         "shift equation [constant,gamma0,gamma2]");
  AddPar("Z4secondOy_lapse", "constant", 
	 "lapse equation [constant,1+log,harmonic;withshift]");
  return 0;
}
