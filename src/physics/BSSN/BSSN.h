/* BSSN.h */
/* Bernd Bruegmann, 6/02 */
/* Wolfgang Tichy  4/2004 */



int BSSN_startup(tL *level);

void BSSN_init(tL* level, int igb, int iK, int ipsi, int igt, int iAt,
	       int iG, int itrK, int iphi,
	       int i_alpha, int i_alphaDensity);

int BSSNtoADM(tL *level);

