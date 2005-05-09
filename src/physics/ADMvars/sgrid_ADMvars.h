/* sgrid_ADMvars.h */
/* Bernd Bruegmann, 6/02 */

extern tVarList *psiandderivs, *K_initial;

/* Helper functions to set partial derivatives of a symmetric tensor */
void FirstDerivsOf_Sab(tBox *box, int i_Sab, int i_dSabc);
void FirstAndSecondDerivsOf_Sab(tBox *box, int i_Sab,
                                int i_dSabc, int i_ddSabcd);
