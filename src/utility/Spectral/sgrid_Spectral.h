/* sgrid_Spectral.h */
/* Wolfgang Tichy */

#include "Spectral.h"
#include "OpenMPsave_loops.h"

/* the kind of FFTs we have */
enum
{
  MATRIX_MULTIPLICATION,
  NUMREC_FFT,
  FFTW3_FFT,
  NTRANSFORMTYPES
};


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
void spec_Int1(tBox *box, int direc, double *u, double *U);


/* Functions from diffmatrices.c */
void initdiffmatrix(tBox *box, int dir, double *D, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_deriv)(void *, double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );
void initdiffmatrix2(tBox *box, int dir, double *DD, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );
void init_fdcentered_diffmatrix(double *x, double *D, int n1,
              void (*fd_deriv)(double *, double *,double *, int) );
void convert_grid_to_fd(tGrid *grid);
void initIntegrationMatrix(double a, double b, double *Int, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*coeffs_of_int)(double, double, double *,double *, int),
                    void (*eval_onPoints)(double *,double *, int) );

/* Functions from matrices.c */
void matrix_times_vector(double *M, double *u, double *Mu, int n);
void matrix_times_matrix(double *M, double *D, double *MD, int n);
void vector_times_matrix(double *u, double *M, double *uM, int n);
double scalarproduct_vectors(double *v, double *w, int n);

/* Functions from filters.c */
void initfiltermatrix(double *F, int k, int n1,
                    void (*get_coeffs)(double *,double *, int),
                    void (*filter_coeffs)(double *, int, int),
                    void (*eval_onPoints)(double *,double *, int) );
void spec_filter1(tBox *box, int direc, double *u);
void spec_filter3d_inbox(tBox *box, int vind, int cind,
                         int nf1, int nf2, int nf3);

/* Functions from spec_coeffs.c */
void initMatrix_ForCoeffs(double *M, int n1,
                          void (*get_coeffs)(double *,double *, int));
void initMatrix_ToEvaluate(double *M, int n1,
                           void (*eval_onPoints)(double *,double *, int));
void spec_analysis1_diffmatrix(tBox *box, int direc, double *M, double *u, double *c);
void spec_synthesis1_diffmatrix(tBox *box, int direc, double *M, double *u, double *c);
void spec_analysis1(tBox *box, int direc, double *u, double *c);
void spec_synthesis1(tBox *box, int direc, double *u, double *c);
void spec_analysis1_inplaneN(tBox *box, int direc, int N, int p,
                             double *u, double *c);
void spec_synthesis1_inplaneN(tBox *box, int direc, int N, int p,
                              double *u, double *c);
void set_TransformType_flags_inbox(tBox *box);
void init_spec_functionpointers(tBox *box);
void get_spec_functionpointers(tBox *box, int direc,
     void (**get_coeffs)(double *,double *, int),
     void (**coeffs_of_deriv)(void *, double, double, double *,double *, int),
     void (**coeffs_of_2ndderiv)(void *, double, double, double *,double *, int),
     void (**coeffs_of_int)(double, double, double *,double *, int),
     void (**eval_onPoints)(double *,double *, int),
     void (**filter_coeffs)(double *, int, int),
     double (**basisfunc)(void *aux, double a, double b, int k, int N, double X) );
void get_spec_functionpointerTO_get_coeffs(tBox *box, int direc,
                               void (**get_coeffs)(double *,double *, int));
void spec_Basis_times_CoeffMatrix(double a, double b, int n,
                                  double *BM, double X,
                    void   (*get_coeffs)(double *,double *, int),
                    double (*basisfunc)(void *aux, double a, double b, int k, int n1, double X));
void spec_Basis_times_CoeffMatrix_direc(tBox *box, int dir, 
                                        double *BM, double X);

/* Functions from integrals.c */
void spec_Integral1(tBox *box, int direc, double *u, double *U);
void spec_2dIntegral(tBox *box, int norm, double *u, double *U);
void spec_SurfaceIntegral(tBox *box, int ig, int norm, double *u, double *U);
double spec_3dIntegral(tBox *box, double *u, double *U);
double BoxVolumeIntegral(tBox *box, int vind);
double GridVolumeIntegral(tGrid *grid, int vind);
void spec_sphericalDF2dIntegral(tBox *box, double *u, double *U);
double spec_sphericalDF3dIntegral(tBox *box, double *u, double *U);
void spec_sphericalDF2dIntegral_at_radial_index_i(tBox *box, double *u, double *U, int i);

/* Functions from interpolate.c */
void spec_Coeffs(tBox *box, double *u, double *c);
void spec_Eval(tBox *box, double *u, double *c);
double spec_interpolate(tBox *box, double *c, double X, double Y, double Z);
void spec_Coeffs_inplaneN(tBox *box, int N, int p, double *u, double *c);
void spec_Eval_inplaneN(tBox *box, int N, int p, double *u, double *c);
double spec_interpolate_inplaneN(tBox *box, int N, int p, double *c,
                                 double X1, double X2);
double spec_interpolate_in_dir_at_i1_i2(tBox *box, int dir, int i1, int i2,
                                        double *c, double X3);
void spec_Coeffs_varlist(tBox *box, tVarList *vlu, tVarList *vlc);
void spec_Eval_varlist(tBox *box, tVarList *vlu, tVarList *vlc);
void spec_interpolate_Var_from_grid2_to_grid1(tGrid *grid1, tGrid *grid2,
                                              int vind, int tempind);

/* from SphericalHarmonics.c */
double *alloc_Plm_Tab(int lmax);
void set_YlmTabs(int lmax, double th, double ph, double *ReYtab,double *ImYtab);
void Ylm_from_Tabs(int lmax, double *ReYtab, double *ImYtab, int l, int m,
                  double *ReYlm, double *ImYlm);
void SphHarm_dphi_forRealFunc(double *c, double *cdphi, int lmax);
void SphHarm_sin_theta_dtheta_forRealFunc(double *c, double *csdth, int lmax);
