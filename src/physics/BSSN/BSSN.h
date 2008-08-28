/* BSSN.h */
/* Bernd Bruegmann, 6/02 */
/* Wolfgang Tichy  4/2004 */



int BSSN_startup(tGrid *grid);

void BSSN_init(tGrid* grid, int igb, int iK, int ipsi, int igt, int iAt,
	       int iG, int itrK, int iphi,
	       int i_alpha, int i_alphaDensity);

int BSSNtoADM(tGrid *grid);

int filter_VarList(tVarList *vl);
int BSSN_filter(tGrid *grid);
void BSSN_filter_unew(tVarList *unew, tVarList *upre);
void BSSN_naive_Ylm_filter(tVarList *unew, tVarList *upre);
void filter_with2o3rule_inX(tVarList *unew, tVarList *upre);
void BSSN_ChooseAndApplyFilters(tVarList *vl);
