/* sgrid_Coordinates.h */
/* Wolfgang Tichy, 4/2005 */

#include "Coordinates.h"

/* Coordinates.c */
int init_CoordTransform_And_Derivs(tGrid *grid);

/* cartesianDerivs.c */
void cart_partials(tBox *box, double *u, double *u1, double *u2, double *u3);
void cart_partial_all(tBox *box, double *u, double *u1, double *u2, double *u3,
                      double *u11, double *u12, double *u13,
                      double *u22, double *u23, double *u33 );

/* doubleCovering.c */
void reset_doubleCoveredPoints(tVarList *unew);

/* coordFilters.c */
void coordinateDependentFilter(tVarList *unew);

/* from ComplexFunctions.c */
double Arg(double x, double y);
double Retanh(double x, double y);
double Imtanh(double x, double y);
double Argtanh(double x, double y);
double Abstanh(double x, double y);
double Argsech(double x, double y);
double Abssech(double x, double y);
