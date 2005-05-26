/* sgrid_Spectral.h */
/* Wolfgang Tichy */

#include "Spectral.h"

/* Functions from derivs.c*/
void cheb_d1(double a, double b, double *u, double *d1u, int n1);
void cheb_d2(double a, double b, double *u, double *d1u, double *d2u, int n1);
void get_memline(double *u, double *line, int direc, int i1, int i2, int n1, int n2, int n3);
void put_memline(double *u, double *line, int direc, int i1, int i2, int n1, int n2, int n3);
void cheb_Deriv1(tBox *box, int direc, double *u, double *du);
void cheb_allDerivs(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 );

void spec_Deriv1(tBox *box, int direc, double *u, double *du);
void spec_allDerivs(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 );


/* Functions from diffmatrices.c */
void initdiffmatrix(double a, double b, double *D, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_deriv)(double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );
void initdiffmatrix2(double a, double b, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_2ndderiv)(double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );
void diffmat_deriv(double *D, double *u, double *du, int n);

/* Functions from filters.c */
void initfiltermatrix(double *F, int k, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*filter_coeffs)(double *, int, int),
                    void (*eval_onPoints)(double *,double *, int) );
void filtermatrix(double *F, double *u, double *uf, int n);
void spec_filter1(tBox *box, int direc, double *u);
