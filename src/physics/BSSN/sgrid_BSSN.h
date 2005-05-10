/* sgrid_BSSN.h */
/* Bernd Bruegmann, 6/02 */
/* Wolfgang Tichy  4/2004 */


void BSSN_rhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur);

void BSSN_init(tL* level, int igb, int iK, int ipsi, int igt, int iAt,
	       int iG, int itrK, int iphi,
	       int i_alpha, int i_alphaDensity);

