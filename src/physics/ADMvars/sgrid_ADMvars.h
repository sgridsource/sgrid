/* sgrid_ADMvars.h */
/* Bernd Bruegmann, 6/02 */

extern tVarList *psiandderivs, *K_initial;

/* Helper functions to set partial derivatives of a symmetric tensor */
void FirstDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc);
void FirstAndSecondDerivsOf_Sab(tBox *box, int i_Sab,
                                int i_dSabc, int i_ddSabcd);
void allDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc, int i_ddSabcd);

/* Helper functions to set partial derivatives of a vector */
void FirstDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab);
void FirstAndSecondDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc);
void allDerivsOf_Sa(tBox *box, int i_Sa, int i_dSab, int i_ddSabc);

/* Helper functions to set partial derivatives of a scalar */
void FirstDerivsOf_S(tBox *box, int i_S, int i_dSa);
void FirstAndSecondDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab);
void allDerivsOf_S(tBox *box, int i_S, int i_dSa, int i_ddSab);

int set_K_initial(tGrid *grid);
int CheckIfFinite(tGrid* grid, char *varname);
