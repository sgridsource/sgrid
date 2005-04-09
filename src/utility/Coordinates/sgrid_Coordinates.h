/* sgrid_Coordinates.h */
/* Wolfgang Tichy, 4/2005 */

#include "Coordinates.h"

/* Coordinates.c */
int init_CoordTransform_And_Derivs(tGrid *grid);

/* cartesianDerivs.c */
void cart_partials(tBox *box, double *u, double *u1, double *u2, double *u3);
