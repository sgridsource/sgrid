/* sgrid_BSSN.c */
/* Wolfgang Tichy  4/2004, Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "BSSN.h"




int sgrid_BSSN(void) 
{
  if (!Getv("physics", "BSSN")) return 0;
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

  /* var we add to RHS of lapse */
  AddVar("BSSN_alphaRHSterm", "", "term we can add to RHS of alpha");

  /* parameters */
  AddPar("BSSN_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("BSSN_reset_doubleCoveredPoints", "no",
         "whether we reset double covered points after each evo step [no,yes]");
  AddPar("BSSN_filter_vars", "no",
         "which vars we filter in each evo substep [no,all]");
  AddPar("BSSN_filter_type", "none",
         "how we filter [old_spec_filters,coordinateDependentFilter,"
         "simple,naive_Ylm,Ylm_lmshift,X2/3,XYZ2/3,X2/3-1]");
  AddPar("BSSN_filter_time", "afterBC", 
         "when we filter in evo substep [afterRHS,afterBC] or [POST_EVOLVE]");
  if(Getv("BSSN_filter_time", "POST_EVOLVE"))
  {
    printf("scheduling BSSN_filter in bin POST_EVOLVE.\n");
    AddFun(POST_EVOLVE, BSSN_filter, "filter all BSSN variables");
  }
  AddPar("BSSN_filter_n2frac", "0.66666666667",
         "fraction of the n2 coeffs to keep, when filtering Y-direc.");
  AddPar("BSSN_filter_n3frac", "0.66666666667",
         "fraction of the n3 coeffs to keep, when filtering Z-direc.");
  AddPar("BSSN_filter_shift2", "0",
         "shift index of last coeff to keep, when filtering Y-direc.");
  AddPar("BSSN_filter_shift3", "0",
         "shift index of last coeff to keep, when filtering Z-direc.");
  AddPar("BSSN_filter_YZregion", "square",
         "region outside which we apply filters [ellipse,square]");
  AddPar("BSSN_densitizedLapse", "no", 
         "whether we evolve a densitized lapse instead of ADM alpha "
         "[no,yes,1+log_withoutShift]");
  AddPar("BSSN_alphaDensityWeight",  "1.0", "weight of densitized lapse");
  
  AddPar("BSSN_forceKzero",    "no",  "set K identically to zero");
  AddPar("BSSN_recomputeAzz",  "no",  "recompute Azz from the other As");
  AddPar("BSSN_recomputegzz",  "no",  "recompute gzz from the other gs");
  AddPar("BSSN_subtractA",     "yes", "set trace of A identically zero");
  AddPar("BSSN_normalizedetg", "no",  "normalize determinant of gamma to one");
  AddPar("BSSN_enforce_AlgConstr", "no", "whether we enforce algebraic "
	 "constraints after evolution step [no,yes]");
  if(Getv("BSSN_enforce_AlgConstr", "yes"))
    AddFun(PRE_POST_EVOLVE, BSSN_enforce_AlgConstr, "enforce algebraic constraints");
  AddPar("BSSN_YoTermFactor",  "0",
         "YoTermFactor = (xi + 2/3) in Yo's paper (gr-qc/0209066)");
  AddPar("BSSN_GReplacedBydg", "yes", "replace G by dg in RHS of evo eqns");
  AddPar("BSSN_freezeGamma",   "no",  "keep BSSN_G constant during to evo");
  AddPar("BSSN_RtoRminusHfactor", "1",
         "how much of Hamiltonian H we subtract from Ricci scalar R");
  AddPar("BSSN_GentlePhiRHS",  "no", "use H as evo eqn as in Gentle et al");

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
	 "lapse equation [constant,1+log,1+log2,harmonic;"
         "withshift,addliealpha,addalphaRHSterm]");
  /* addalphaRHSterm: whether we add BSSN_alphaRHSterm to RHS of lapse eqn */
  AddPar("BSSN_set_alphaRHSterm", "harmonic0",
         "set BSSN_alphaRHSterm [harmonic0] "
         "harmonic0 means set BSSN_alphaRHSterm s.t. lapse RHS=0 at t=0");

  AddPar("BSSN_shift_stop_time", "-1.0", 
         "time when shift stops evolving (-1 for don't stop)");
         
  return 0;
}
