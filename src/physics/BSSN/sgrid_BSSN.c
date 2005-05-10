/* sgrid_BSSN.c */
/* Bernd Bruegmann 6/02 */
/* Wolfgang Tichy  4/2004 */

#include "sgrid.h"
#include "BSSN.h"




void sgrid_BSSN(void) 
{
  if (!Getv("physics", "BSSN")) return;
  printf("Adding BSSN\n");

  /* functions */
  AddFun(POST_INITIALDATA, BSSN_startup, 
	 "initialize BSSN system from adm initial data");
  AddFun(POST_EVOLVE, BSSNtoADM,
         "fill in ADM variables from BSSN variables");
  // AddFun(ANALYZE,     BSSNconstraints,   "compute BSSN constraints");

  /* variables */
  AddVar("BSSN_g",     "(ij)", "conformal metric");
  AddVar("BSSN_phi",   "",      "conformal factor");
  AddVar("BSSN_A",     "(ij)", "tracefree extrinsic curvature");
  AddVar("BSSN_K",     "",      "trace of extrinsic curvature");
  AddVar("BSSN_G",     "I",     "contracted Gamma");
  AddVar("BSSN_alphaDensity", "", "densitized Lapse");

  /* derivatives which need to be precomputed before each evo steo */
  /* Note: ADMvars_ddg, ADMvars_dg, ADMvars_dA are overwritten and used 
           in place of BSSN_ddg, BSSN_dg, BSSN_dA */
  AddVar("BSSN_dphi",   "i",    "1st deriv of conformal factor");
  AddVar("BSSN_ddphi",  "(ij)", "2nd deriv of conformal factor");
  AddVar("BSSN_dK",     "i",    "1st deriv of trace of extrinsic curvature");
  AddVar("BSSN_dG",     "Ij",   "1st deriv of contracted Gamma");
  AddVar("BSSN_dalp",   "i",    "1st deriv of Lapse");
  AddVar("BSSN_ddalp",  "(ij)", "2nd deriv of Lapse");
  AddVar("BSSN_dbeta",  "Ij",   "1st deriv of shift");
  AddVar("BSSN_ddbeta", "I(jk)","2nd deriv of shift");
    
  /* parameters */
  AddPar("BSSN_densitizedLapse", "no", 
  "wether we evolve a densitized lapse instead of ADM alpha "
  "[no,yes,1+log_withoutShift]");
  AddPar("BSSN_alphaDensityWeight",  "1.0", "weight of densitized lapse");
  
  AddPar("BSSN_forceKzero",    "no",  "set K identically to zero");
  AddPar("BSSN_subtractA",     "yes", "set trace of A identically zero");
  AddPar("BSSN_normalizedetg", "no",  "normalize determinant of gamma to one");

  AddPar("BSSN_lapsepsipower",  "0", "power of psi in lapse equation");
  AddPar("BSSN_lapseharmonicf", "2", "2 for 1+log, 1 for harmonic");
  AddPar("BSSN_subtractK0",     "yes", 
	 "for harmonic and 1+log, subtract initial traceK or not");

  AddPar("BSSN_shiftpsipower",   "0", "power of 1/psi in shift equation");
  AddPar("BSSN_shiftalphapower", "0", "power of lapse in shift equation");
  AddPar("BSSN_shiftgammacoeff", "0.75", "coefficient in shift equation");
  AddPar("BSSN_shiftdriver",     "0", "coefficient of diffusion term");
         
  AddPar("BSSN_shift", "constant", "shift equation [constant,gamma0,gamma2]");
  AddPar("BSSN_lapse", "constant", 
	 "lapse equation [constant,1+log,1+log2,harmonic;withshift]");

  AddPar("BSSN_shift_stop_time", "-1.0", 
         "time when shift stops evolving (-1 for don't stop)");
}
