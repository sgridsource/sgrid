/* TestID.h */
/* Wolfgang Tichy 09/2009 */


int TestID(tGrid *grid);
void TestID_1dGaugeWave(tGrid *grid, int i_x, int i_gb, int i_K, int i_psi, 
                       int i_dpsiopsi, int i_ddpsiopsi, 
                       int i_alpha, int i_beta);
void TestID_1dLinearWave(tGrid *grid, int i_x, int i_gb, int i_K, int i_psi, 
                       int i_dpsiopsi, int i_ddpsiopsi, 
                       int i_alpha, int i_beta);
void TestID_RandomNoise(tGrid *grid, int i_x, int i_gb, int i_K, int i_psi, 
                        int i_dpsiopsi, int i_ddpsiopsi, 
                        int i_alpha, int i_beta);
