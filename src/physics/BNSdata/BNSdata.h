/* BNSdata.h */
/* Wolfgang Tichy, 2007 */


/* main functions */
int BNSdata_setup_boxsizes(tGrid *grid);
int BNSdata_startup(tGrid *grid);
int BNSdata_analyze(tGrid *grid);
int BNSdata_solve(tGrid *grid);
int setBNSdata(tGrid *grid);
void setADMvars(tGrid *grid);


/* funtions from mathematica */
void BNS_CTS(tVarList *vlFu, tVarList *vlu, tVarList *vlJdu, 
             tVarList *vldu, tVarList *vlduDerivs, int nonlin);
void BNS_compute_new_q(tGrid *grid);
double BNS_compute_new_q_atXYZ(tGrid *grid, int bi,
                               double X, double Y, double Z);

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
double InnerVolumeIntergral(tGrid *grid, int b, int vind);
void adjust_box4_5_pars(tGrid *grid);
void Interpolate_Var_From_Grid1_To_Grid2(tGrid *grid1, tGrid *grid2, int vind);
double BNS_update_q_atXYZ(tGrid *grid2, 
                          int b2, double X2, double Y2, double Z2,
                          double w, tGrid *grid1);
void BNS_update_q(tGrid *grid2, double w, tGrid *grid1);
