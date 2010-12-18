/* BNSdata.h */
/* Wolfgang Tichy, 2007 */


/* main functions */
int BNSdata_setup_boxsizes(tGrid *grid);
int BNSdata_startup(tGrid *grid);
void   BNS_compute_new_centered_q(tGrid *grid);
double BNS_compute_new_centered_q_atXYZ(tGrid *grid, int bi,
                                        double X, double Y, double Z);
int BNSdata_verify_solution(tGrid *grid);
int BNSdata_analyze(tGrid *grid);
int BNSdata_solve(tGrid *grid);
int setBNSdata(tGrid *grid);
void setADMvars(tGrid *grid);
void set_BNSdata_ABphi(tGrid *grid);


/* funtions from mathematica */
void BNS_CTS(tVarList *vlFu, tVarList *vlu, tVarList *vlJdu, 
             tVarList *vldu, tVarList *vlduDerivs, int nonlin);
void BNS_compute_new_q(tGrid *grid, int iq);
double BNS_compute_new_q_atXYZ(tGrid *grid, int bi,
                               double X, double Y, double Z);
void BNS_set_restmassintegrand(tGrid *grid, int iInteg);
void BNS_set_J_ADM_VolInt_integrand(tGrid *grid, int iInteg);
void BNS_set_M_ADM_VolInt_integrand(tGrid *grid, int iInteg);
void set_BNSdata_Sigma_BC(tVarList *vlFu, tVarList *vlu,
                          tVarList *vlJdu, tVarList *vldu,
                          tVarList *vlduDerivs, int nonlin);

/* for solving all ell. eqns together */
void F_BNSdata(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_BNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);

/* for solving the ell. eqns sequentially */
void F_oneComp(tVarList *vlFw, tVarList *vlw,
               tVarList *vlwDerivs, tVarList *vlc2);
void J_oneComp(tVarList *vlJdw, tVarList *vldw,
               tVarList *vldwDerivs, tVarList *vlw);

/* funcs from TOV */
int TOV_init(double Pc, double kappa, double Gam, int pr, double *rf_surf,
             double *m, double *Phi_c, double *Psi_c, double *m0);
int TOV_m_P_Phi_Psi_m0_OF_rf(double rf, double rf_surf,
                          double kappa, double Gam,
                          double Pc, double Phic, double Psic,
                          double *m, double *P, double *Phi, double *Psi,
                          double *m0);

/* from BNSgrid.c */
int set_boxsizes(tGrid *grid);
int set_sigma_pm_vars(tGrid *grid);
void reset_Coordinates_AnsorgNS_sigma_pm(tGrid *grid, tGrid *gridnew,
                                         int innerdom,  int outerdom);
double ADMmass_fromPsi_inbox1_at_A1B0(tGrid *grid, int iADMmass);
double InnerVolumeIntegral(tGrid *grid, int b, int vind);
double VolumeIntegral_inBNSgridBox(tGrid *grid, int b, int vind);
tGrid *make_grid_with_sigma_pm(tGrid *grid, int nAB, int nphi, int nxyz);
void adjust_box4_5_pars(tGrid *grid);
int BNSgrid_Get_BoxAndCoords_of_xyz(tGrid *grid1,
                                    double *X1, double *Y1, double *Z1,
                                    int b, double x, double y, double z);
void Interpolate_Var_From_Grid1_To_Grid2(tGrid *grid1, tGrid *grid2, int vind);
void Interp_Var_From_Grid1_To_Grid2_pm(tGrid *grid1, tGrid *grid2, int vind,
                                       int innerdom);
void Interpolate_Var_From_Grid1_To_Grid2_wrapper(tGrid *grid1, tGrid *grid2,
                                                 int vind, int dummy);
void copy_Var_at_i0_from_Box1_Box2(tGrid *grid, int vind, int b1, int b2);
void BNS_enforce_uniqueness_on_axis(tVarList *vlu);
double BNS_update_q_atXYZ(tGrid *grid2, 
                          int b2, double X2, double Y2, double Z2,
                          double w, tGrid *grid1);
void BNS_update_q(tGrid *grid2, double w, tGrid *grid1);
void BNSgrid_init_Coords(tGrid *grid);
void BNSgrid_copy_DomainShape(tGrid *grid, int ibd);
void BNSgrid_set_Var_equalmasses_sym(tGrid *grid, int ibd, int iv, int sym);
void BNSgrid_set_allVars_onLeft_equalmasses(tGrid *grid);

/* from BNS_Interpolate_ADMvars.c */
int BNS_Interpolate_ADMvars(tGrid *grid);

/* from BNS_BCs.c */
void set_BNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);
void BNSdata_RegularityConditions_for_Var_at_rho_eq_0(tBox *box, double *FPsi,
                        double *Psi, double *Psix, double *Psiy, double *Psiz);
