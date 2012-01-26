/* sgrid_BNSdata.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "BNSdata.h"


int sgrid_BNSdata() 
{
  if (!Getv("physics", "BNSdata")) return 0;
  printf("Adding BNSdata\n");

  /* functions */
  AddFun(PRE_GRID, set_boxsizes, "setup initial box sizes");
  AddFun(POST_GRID, set_sigma_pm_vars, "setup the sigma_{+-} vars from AnsorgNS");
  AddFun(PRE_INITIALDATA, BNSdata_startup, "initialize BNSdata");
  AddFun(INITIALDATA, setBNSdata, "set the BNS data");
  AddFun(ANALYZE, BNSdata_analyze, "compute properties of BNS data");
  AddFun(OUTPUT, BNS_Interpolate_ADMvars,
         "interpolate ADM initial data using BNSdata_Interpolate_pointsfile");

  /* variables */
  AddVar("BNSdata_Psi",     "",     "new conf. factor");
  AddVar("BNSdata_Psi",     "i",    "1st deriv of Psi");
  AddVar("BNSdata_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("BNSdata_B",       "I",    "beta - Omega cross r");
  AddVar("BNSdata_B",       "Ij",   "1st deriv of B^i");
  AddVar("BNSdata_B",       "I(jk)","2nd deriv of B^i");
  AddVar("BNSdata_alphaP",  "",     "lapse times Psi");
  AddVar("BNSdata_alphaP",  "i",    "1st deriv of alphaP");
  AddVar("BNSdata_alphaP",  "(ij)", "2nd deriv of alphaP");
  AddVar("BNSdata_Sigma",   "",     "Sigma is potential for irrot. part of h*u^i");
  AddVar("BNSdata_Sigma",   "i",    "1st deriv of Sigma");
  AddVar("BNSdata_Sigma",   "(ij)", "2nd deriv of Sigma");
  AddVar("BNSdata_SigmaX",    "",     "1st A-deriv of Sigma");
  AddVar("BNSdata_SigmaXX",   "",     "2nd A-deriv of Sigma");
  AddVar("BNSdata_SigmaXXX",  "",     "3rd A-deriv of Sigma");
  AddVar("BNSdata_lSigmaX",   "",     "1st A-deriv of Sigma");
  AddVar("BNSdata_lSigmaXX",  "",     "2nd A-deriv of linearized Sigma");
  AddVar("BNSdata_lSigmaXXX", "",     "3rd A-deriv of linearized Sigma");

  AddVar("BNSdata_VR",  "I", "velocity in rotating frame: "
                             "V^i = u^i/u^0 - l^i,  l^i = HKV");
  AddVar("BNSdata_wB",  "I", "conformal rotational part of h*u^i: "
                             "h u^i = w^i + D^i Sigma - h u^0 beta^i, "
                             "w^i = Psi^{-6} wB^i, D^i Sigma = g^{ij} d_i Sigma. "
                             "Note: v^i = u^i/u^0, i.e. it is not u^i/(alpha u^0)");
  AddVar("BNSdata_q",   "",  "q := P/rho0");
  AddVar("BNSdata_wB",  "Ij","1st deriv of wB");
  AddVar("BNSdata_q",   "i", "1st deriv of q");

  /* star surfaces at current and previous iterations */
  AddVar("BNSdata_surface_sigma_pm", "",
         "previous Coordinates_AnsorgNS_sigma_pm. We save only in box0/3 and "
         "keep previous surfaces at indices i=1,2,... current surf is at i=0");

  /* sometimes we save the old vars before the ell. solve */
  AddVar("BNSdata_Psiold",    "",  "old Psi");
  AddVar("BNSdata_Bold",      "I", "old B");
  AddVar("BNSdata_alphaPold", "",  "old alphaP");
  AddVar("BNSdata_Sigmaold",  "",  "old Sigma");
  AddVar("BNSdata_qold",      "",  "old q");
  /* q we get with corotating formula */
  AddVar("BNSdata_qcorot",    "",  "q we get with corotating formula");
  AddVar("BNSdata_qnocent",   "",  "q we get without centering");

  AddVar("BNSdata_A", "", "store value of A in box4/5");
  AddVar("BNSdata_B", "", "store value of B in box4/5");
  AddVar("BNSdata_phi", "", "store value of phi in box4/5");

  AddVar("BNSdata_temp1", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp2", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp3", "", "temporary variable(e.g. to store derivs)");
  AddVar("BNSdata_temp4", "", "temporary variable(e.g. to store derivs)");
     
  /* parameters */
  AddPar("BNSdata_rotationstate1", "corotation",
         "rotation state of NS1 [corotation,irrotation");
  AddPar("BNSdata_rotationstate2", "corotation",
         "rotation state of NS2 [corotation,irrotation");
  AddPar("BNSdata_wB_factor", "1", "wB = BNSdata_omega cross (r-r_c) * factor"
         " [1,Psi6,h,1/alpha]");
  AddPar("BNSdata_omegax1", "0", "x-comp of angular velocity of NS1");
  AddPar("BNSdata_omegay1", "0", "y-comp of angular velocity of NS1");
  AddPar("BNSdata_omegaz1", "0", "z-comp of angular velocity of NS1");
  AddPar("BNSdata_omegax2", "0", "x-comp of angular velocity of NS2");
  AddPar("BNSdata_omegay2", "0", "y-comp of angular velocity of NS2");
  AddPar("BNSdata_omegaz2", "0", "z-comp of angular velocity of NS2");
  AddPar("BNSdata_m01",   "0.141202", "rest mass of NS1");
  AddPar("BNSdata_m02",   "0.141202", "rest mass of NS2");
  AddPar("BNSdata_iterate_m0", "no", "whether we iterate rest masses [no,yes]");
  if(Getv("BNSdata_iterate_m0", "yes"))
  {
    AddPar("BNSdata_init_m01", "-1", "initial rest mass of NS1, if positive "
           "BNSdata_m01 will be changed during the iterations until we reach "
           "BNSdata_desired_m01");
    AddPar("BNSdata_init_m02", "-1", "initial rest mass of NS2, if positive "
           "BNSdata_m02 will be changed during the iterations until we reach "
           "BNSdata_desired_m02");
    AddPar("BNSdata_m0change", "0.1", "amount by which change BNSdata_m01/2 "
           "during iterations");
    AddPar("BNSdata_desired_m01", "-1", "desired rest mass1 if we iterate it");
    AddPar("BNSdata_desired_m02", "-1", "desired rest mass2 if we iterate it");
    AddPar("BNSdata_desired_kappa", "-1", "desired kappa in EOS");
    AddPar("BNSdata_adjustkappa", "no", 
           "how we adjust kappa [kappa^(-n/2)*m0_EQ_desired_kappa^(-n/2)*m0,no]");
  }
  AddPar("BNSdata_qmax1", "0", "max q of NS1 along x-axis");
  AddPar("BNSdata_qmax2", "0", "max q of NS2 along x-axis");
  AddPar("BNSdata_xmax1", "0", "pos. of max q of NS1 along x-axis");
  AddPar("BNSdata_xmax2", "0", "pos. of max q of NS2 along x-axis");
  AddPar("BNSdata_actual_xmax1", "0", "x-pos. of actual global max of q in NS1");
  AddPar("BNSdata_actual_ymax1", "0", "y-pos. of actual global max of q in NS1");
  AddPar("BNSdata_actual_zmax1", "0", "z-pos. of actual global max of q in NS1");
  AddPar("BNSdata_actual_xmax2", "0", "x-pos. of actual global max of q in NS2");
  AddPar("BNSdata_actual_ymax2", "0", "y-pos. of actual global max of q in NS2");
  AddPar("BNSdata_actual_zmax2", "0", "z-pos. of actual global max of q in NS2");
  AddPar("BNSdata_center_new_q", "no",
         "if and how we center new q on (BNSdata_xmax1/2,0,0) "
         "[no,center_yz,center_xyz,adjust_domainshapes]");
  AddPar("BNSdata_center_new_q_timebin", "before_ell_solve", 
         "when we center q "
         "[before_ell_solve,after_adjusting_Omega_xCM]");
  AddPar("BNSdata_center_new_q_first_at", "0", "first iteration when we "
         "center q [#,-1]. -1 means never");
  AddPar("BNSdata_center_new_q_flag", "no",
         "set automatically in main iteration [no,yes]");
  AddPar("BNSdata_center_new_q_fac", "0.1",
         "by how much we try to move the the qmax location");
  AddPar("BNSdata_Omega_init", "BNSdata_Omega", "ini. orbital angular velocity "
         "[#,estimate,estimate_from_desired_m0,BNSdata_Omega]");
  AddPar("BNSdata_x_CM_init", "BNSdata_x_CM",
         "initial center of mass in x-direction [#,estimate,BNSdata_x_CM]");
  AddPar("BNSdata_Omega", "estimate",  "orbital angular velocity "
         "[#,estimate,estimate_from_desired_m0]");
  AddPar("BNSdata_b",     "1",  "separation parameter (distance~2b)");
  AddPar("BNSdata_n",     "1",  "polytropic index n, Gamma = 1 + 1/n");
  AddPar("BNSdata_kappa", "1",  "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("BNSdata_x_CM",  "estimate",  "center of mass in x-direction");
  AddPar("BNSdata_C1",    "-1", "C1 in q = (C1/F-1)/(n+1) "
         "[needs to be adjusted so that m01 stays the constant]");
  AddPar("BNSdata_C2",    "-1", "C2 in q = (C2/F-1)/(n+1) "
         "[needs to be adjusted so that m02 stays the constant]");
  AddPar("BNSdata_guess", "TOV", "init. guess for Psi, alphaP and "
         "q [initialize_from_checkpoint,TOV,TOVaverage,TOVproduct] "
         "for shift [TaniguchiShift]");
  AddPar("BNSdata_initfile", "", "name of first initialization file to read");
  AddPar("BNSdata_init_q_fromfields", "no", "init q from other fields [no,yes]");
  if(Getd("BNSdata_m02")==0.0)
  {
    AddPar("BNSdata_yshift1", "0", "shift NS1 in y-direction for testing");
    AddPar("BNSdata_adjustdomain01", "yes", "if we adjust domainshapes "
           "after shift [yes,no]");
  }
  AddPar("BNSdata_Sigma_surface_BCs", "AddInnerVolIntToBC RegularityOnAxis",
         "BCs for Sigma on star surfaces "
         "[AddInnerVolIntToBC,AddInnerSumToBC,RegularityOnAxis,"
         "ZeroAtPoint,InnerVolIntZero,InnerSumZero,AtA0B0,AtA0B1,"
         "none,ZeroInOuterBoxes]");
  AddPar("BNSdata_KeepInnerSigma", "no", "keep Sigma in inner boxes [no,yes]");
  AddPar("BNSdata_useDD", "no",
         "whether we use the DD ops to compute second derivs [no,yes]");
  AddPar("BNSdata_itmax", "10", "max. number of iterations in BNSdata_solve");
  AddPar("BNSdata_break_if_err_below_tol", "at_iterationend after_ell_solve",
         "which if clauses we use to break out of main iteration loop "
         "[at_iterationend,after_ell_solve]");
  AddPar("BNSdata_tol",   "1e-6", "tolerance for BNSdata_solve");
  AddPar("BNSdata_Newton_tolFac", "0.1", "tol for Newton is "
         "Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac)");
  AddPar("BNSdata_esw",   "1", "ell. solve weight: after ell. solve new "
         "and old values are averaged according to: "
         "new = esw*new + (1-esw)*old");
  AddPar("BNSdata_esw1",  "1", "second weight esw1 is used if better");
  AddPar("BNSdata_allow_esw1_first_at", "-1", "first iteration when esw=esw1 "
         "will be tried and allowed if better [#,-1]. -1 means never");
  AddPar("BNSdata_Sigma_esw", "0.2", "ell. solve weight for Sigma");
  AddPar("BNSdata_Sigma_esw1",  "1", "second Sigma weight is used if better");
  AddPar("BNSdata_allow_Sigma_esw1_first_at", "-1",
         "first iteration when Sigma_esw=Sigma_esw1 "
         "will be tried and allowed if better [#,-1]. -1 means never");
  AddPar("BNSdata_reset_qmax_xmax_pars_at", "-1", "list of its when we "
         "reset BNSdata_qmax1/2, BNSdata_xmax1/2");
  AddPar("BNSdata_analyze_xmax", "set_BNSdata_xmax",
         "what we do with BNSdata_xmax1/2 inside BNSdata_analyze "
         "[set_BNSdata_xmax,print_xmax]");
  AddPar("BNSdata_adjust", "nothing", "what we adjust (apart from C1/2) "
         "after ell. solve. E.g. \"keep_xout keep_one_xout\" adjusts Omega "
         "and x_CM to keep either xout1 or xout2 in place. "
         "[nothing,WT_L2_method,keep_xout,fix_xout [keep_one_xout,always],"
         "keep_xmax [keep_one_xmax,always,reset_xmax_if_problem],"
         "keep_xfm [keep_one_xmax,always,reset_xmax_if_problem],"
         "forcebalance [always]"
         "]");
  AddPar("BNSdata_adjust_first_at", "0", 
         "first iteration when we use BNSdata_adjust. -1 means never");
  AddPar("BNSdata_adjust_mintol", "1e-10", "always use tol>=mintol in adjust, "
         "e.g. when BNSdata_adjust=keep_xout");
  AddPar("BNSdata_dOmega_fac", "0.1", "dOmega = Omega*dOmega_Fac");
  AddPar("BNSdata_dx_CM_fac",  "0.1", "dx_CM = b*dx_CM_fac");
  AddPar("BNSdata_domainshape_diff_tol", "3e30", 
         "adjust domain if difference in shape is less than this number."
         " >=1e30 means adjust always.");
  AddPar("BNSdata_domainshape_weight", "1",
         "new and old domainshapes are averaged according to: "
         "new = weight*new + (1-weight)*old");
  AddPar("BNSdata_domainshape_filter", "no", 
         "how/if we filter domainshape sigma_pm [no,Bphi2/3] "
         "Bphi2/3 means use 2/3 rule in B- and phi-dir");
  AddPar("BNSdata_EllSolver_method", "BNS_Eqn_Iterator",
         "how we solve for Psi,B^i,alphaP,Sigma "
         "[allatonce, BNS_Eqn_Iterator, BNS_ordered_Eqn_Iterator,"
         " sequence1, sequence2, sequence3]");
  AddPar("BNSdata_CTS_Eqs_Iteration_order", 
         "BNSdata_Psi BNSdata_Bx BNSdata_By BNSdata_Bz BNSdata_alphaP",
         "Order we use in function BNS_ordered_Eqn_Iterator. "
         "This par has to contain all we solve for!");
  AddPar("BNSdata_linSolver", "UMFPACK", 
         "linear solver used [LAPACK,templates_GMRES,bicgstab,UMFPACK]");
  AddPar("BNSdata_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("BNSdata_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("BNSdata_guess", "test", "initial guess [test]");

  /* pick root finders and their options */
  AddPar("BNSdata_sigp_Bphi_FINDER_reset_Coordinates_AnsorgNS_sigma_pm",
         "newton_linesrch_itsP", "[newton_linesrch_itsP,zbrent_itsP]");

  /* The next 2 pars affect how the overall iteration proceeds. E.g. they
     affect how close x_CM is to zero for equal masses. Set:
     BNSdata_Interpolation_smooths = 0 
     BNSdata_adjust_domainshapes_Grid1_To_Grid2_Interpolator = Interpolate_Var_From_Grid1_To_Grid2_wrapper
     to get the old behavior for corotation */
  AddPar("BNSdata_Interpolation_smooths", "0",
         "number of smooths by Interpolate_Var_From_Grid1_To_Grid2");
  AddPar("BNSdata_adjust_domainshapes_Grid1_To_Grid2_Interpolator",
         "Interp_Var_From_Grid1_To_Grid2_pm",
         "how we interpolate from Grid1 to Grid2 in "
         "compute_new_q_and_adjust_domainshapes "
         "[Interp_Var_From_Grid1_To_Grid2_pm,Interpolate_Var_From_Grid1_To_Grid2_wrapper]");
         
  AddPar("BNSdata_grid", "4ABphi_2xyz",
         "what grid we use [SphericalDF, AnsorgNS, 4ABphi_2xyz]");
  AddPar("BNSdata_box0_Amax", "0.85", "max A we use in box0 [0...1]");
  AddPar("BNSdata_box3_Amax", "0.85", "max A we use in box3 [0...1]");
  AddPar("BNSdata_2xyz_n", "6", "if 4ABphi_2xyz: use n*n*n points in box4/5");
  AddPar("BNSdata_regularization", "none",
         "options for 4ABphi_2xyz and AnsorgNS " 
         "[regularity_on_axis,regularity_on_axis_at_interfaces,"
         "regularity_on_axis_at_center]");
  AddPar("BNSdata_uniqueness_on_axis", "no",
         "set vars equal to val at phi=0 for all phi>0 [no,yes]");

  return 0;
}
