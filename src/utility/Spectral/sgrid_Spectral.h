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
void spec_Deriv2(tBox *box, int direc, double *u, double *du);


/* Functions from diffmatrices.c */
void initdiffmatrix(double a, double b, double *D, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_deriv)(double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );
void initdiffmatrix2(double a, double b, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_2ndderiv)(double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );

/* Functions from matrices.c */
void matrix_times_vector(double *M, double *u, double *Mu, int n);
void matrix_times_matrix(double *M, double *D, double *MD, int n);
void vector_times_matrix(double *u, double *M, double *uM, int n);

/* Functions from filters.c */
void initfiltermatrix(double *F, int k, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*filter_coeffs)(double *, int, int),
                    void (*eval_onPoints)(double *,double *, int) );
void spec_filter1(tBox *box, int direc, double *u);

/* Functions from spec_coeffs.c */
void initMatrix_ForCoeffs(double *M, int n1,
                          void (*get_coeffs)(double *,double *, int));
void initMatrix_ToEvaluate(double *M, int n1,
                           void (*eval_onPoints)(double *,double *, int));
void spec_analysis1(tBox *box, int direc, double *M, double *u, double *c);
void spec_synthesis1(tBox *box, int direc, double *M, double *u, double *c);
void get_spec_functionpointers(tBox *box, int direc,
     void (**get_coeffs)(double *,double *, int),
     void (**coeffs_of_deriv)(double, double, double *,double *, int),
     void (**coeffs_of_2ndderiv)(double, double, double *,double *, int),
     void (**eval_onPoints)(double *,double *, int),
     void (**filter_coeffs)(double *, int, int),
     double (**basisfunc)(void *aux, double a, double b, int k, int n1, double X) );
void get_spec_functionpointerTO_get_coeffs(tBox *box, int direc,
                               void (**get_coeffs)(double *,double *, int));
void spec_Basis_times_CoeffMatrix(double a, double b, int n,
                                  double *BM, double X,
                    void   (*get_coeffs)(double *,double *, int),
                    double (*basisfunc)(void *aux, double a, double b, int k, int n1, double X));

/* Functions from integrals.c */
void spec_Integral1(tBox *box, int direc, double *u, double *U);
void spec_2dIntegral(tBox *box, int norm, double *u, double *U);
double spec_3dIntegral(tBox *box, double *u, double *U);
void spec_sphericalDF2dIntegral(tBox *box, double *u, double *U);
double spec_sphericalDF3dIntegral(tBox *box, double *u, double *U);

/* Functions from interpolate.c */
void spec_Coeffs(tBox *box, double *u, double *c);
void spec_Eval(tBox *box, double *u, double *c);
double spec_interpolate(tBox *box, double *c, double X, double Y, double Z);
