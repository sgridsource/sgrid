/* BNSdata.h */
/* Wolfgang Tichy, 2007 */


int BNSdata_startup(tGrid *grid);
int BNSdata_analyze(tGrid *grid);
int BNSdata_solve(tGrid *grid);
int setBNSdata(tGrid *grid);
void setADMvars(tGrid *grid);

void F_BNSdata(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_BNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);
void Precon_E(tVarList *vlJdu, tVarList *vldu,
              tVarList *vlduDerivs, tVarList *vlu);
void BNS_CTS(tVarList *vlFu, tVarList *vlu, tVarList *vlJdu, 
             tVarList *vldu, tVarList *vlduDerivs, int nonlin);
