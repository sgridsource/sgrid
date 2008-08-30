/* ScalarOnKerr.h */
/* Wolfgang Tichy  8/2005 */


int ScalarOnKerr_startup(tGrid *grid);
int KerrChecker(tGrid *grid);
double ScalarOnKerr_Source(double t, double x,double y,double z);
void ScalarOnKerr_evolve(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur);
void set_psi_Pi_boundary(tVarList *unew, tVarList *upre, double dt, 
                         tVarList *ucur);
void set_psi_Pi_phi_boundary(tVarList *unew, tVarList *upre, double dt,
                             tVarList *ucur);
void filter_unew(tVarList *unew, tVarList *upre);
void naive_Ylm_filter_unew(tVarList *unew, tVarList *upre);
void filter_unew_radially(tVarList *unew, tVarList *upre);
int ScalarOnKerr_analyze(tGrid *grid);
void Kerr(tGrid *grid, int i_x, int i_g, int i_gup, int i_Gam, int i_G);
void Kerr3d(tGrid *grid, int i_x, int i_alpha, int i_beta, int i_g, 
            int i_K, int i_TrK, int i_gup,
            int i_Gam, int i_G, int i_dalpha, int i_dbeta);
void set_boundary_ofPi(tVarList *unew, tVarList *upre);
void set_Up_Um_onBoundary(tVarList *unew, tVarList *upre,
                          double dt, tVarList *ucur);
void compute_unew_from_Up_Um_onBoundary(tVarList *unew, tVarList *upre,
                                        double dt, tVarList *ucur);
void interpolate_between_boxes(tVarList *unew, tVarList *upre);
int ScalarOnKerr_setup_boxes(tGrid *g);
void ScalarOnKerr_evolve_1stO(tVarList *unew, tVarList *upre, double dt,
                              tVarList *ucur);
void set_Up_Um_U0_onBoundary(tVarList *unew, tVarList *upre, double dt, 
                             tVarList *ucur);
void compute_unew_from_Up_Um_U0_onBoundary(tVarList *unew, tVarList *upre,
                                           double dt, tVarList *ucur);
void ChooseAndApplyFilter(tVarList *unew);
int ScalarOnKerr_filter(tGrid *grid);
