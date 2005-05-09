/* sgrid_ADMvars.c */
/* Bernd Bruegmann 6/02 */

#include "sgrid.h"
#include "ADMvars.h"

tVarList *psiandderivs, *K_initial;



int sgrid_ADMvars() 
{
  if (!Getv("physics", "ADMvars")) return;
  printf("Adding ADMvars\n");

  /* functions */
  AddFun(ANALYZE, computeADMconstraints, "compute ADM constraints");

  /* variables */
  AddVar("g",        "(ij)", "metric");
  AddVar("K",        "(ij)", "extrinsic curvature");
  AddVar("alpha",    "",      "lapse");
  AddVar("beta",     "I",     "shift");
  AddVar("betadot",  "I",     "time derivative of shift");

  AddVar("ham",      "",      "Hamiltonian constraint");
  AddVar("mom",      "i",     "momentum constraint");
  AddVar("trK",      "",
	 "trace of extrinsic curvature tensor (for output)");

  AddPar("ADMvars_normalizedConstraints", "no",
         "whether we compute normalized constraints [no,yes]");
  if(Getv("ADMvars_normalizedConstraints", "yes"))
  {
    AddPar("ADMvars_ConstraintNorm", "TermByTerm",
           "how to normalize [TermByTerm, dumb]");
    AddVar("normham",  "",      "normalized Hamiltonian constraint");
    AddVar("normmom",  "i",     "normalized momentum constraint");
  }

  /* time independent conformal factor */
  AddConstantVar("psi",       "",      "time independent conformal factor");
  AddConstantVar("dpsiopsi",  "i",     "(del_i psi) / psi");
  AddConstantVar("ddpsiopsi", "(ij)", "(del_i del_j psi) / psi");
  
  /* this is perhaps not the ideal place for global var list ... */
  psiandderivs = vlalloc(NULL);
  vlpush(psiandderivs, Ind("psi"));
  vlpush(psiandderivs, Ind("dpsiopsix"));
  vlpush(psiandderivs, Ind("ddpsiopsixx"));

  /* 1+log slicing may require the initial value of the trace of K_ij */
  AddVar("K_initial", "",      "initial value of trace of K_ij");
  AddFun(POST_INITIALDATA, set_K_initial, 
	 "compute and save initial value of the trace of K_ij");
  K_initial = vlalloc(NULL);
  vlpush(K_initial, Ind("K_initial"));

  /* check for NANs */
  AddPar("ADMvars_NANcheck", "no", "check if ham contains NAN or INF");
  AddPar("ADMvars_exitifNAN", "gxx", 
	 "exit if NAN is found [varname,no] (not a list)");
  if (!Getv("ADMvars_exitifNAN", "no")) 
    AddFun(POST_EVOLVE, ExitIfNAN, "exit if NAN");

  /* add variables to take derivs */
  AddVar("ADMvars_dg",  "(ij)k", "first partial derivs of a symmetric tensor");
  AddVar("ADMvars_ddg", "(ij)(kl)", "second part. derivs of a symm. tensor");
  AddVar("ADMvars_dK",  "(ij)k", "first partial derivs of a symmetric tensor");

  return 0;
}
