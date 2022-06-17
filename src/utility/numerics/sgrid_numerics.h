/* sgrid_numerics.h */
/* Wolfgang Tichy, June 2022 */


/* newton1d_brak.c */
int newton1d_brak(double *x0,
                  void (*fdf)(double x, void *par, double *f, double *df),
                  double x1, double x2, void *par, int maxits, double xacc,
                  int pr);
int rtbisect(double *x0, double (*func)(double,void *par),
             double x1, double x2, void *par, int maxits, double xacc, int pr);

/* rt_brak.c */
int rt_brak(double (*func)(double x, void *par), void *par,
            double *x1, double *x2, int maxits);

/* rtbrent_brak.c */
int rtbrent_brak(double *x0, double (*func)(double,void *par),
                 double x1, double x2, void *par, int maxits, double xacc,
                 int pr);

/* newton1d_fd.c */
int newton1d_fd_region(double *x0, double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);
int find_2roots_region(double x0[2],
                       double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);

/* lu.c */
int lu_decomp(int n, double A[][n], int *p_idx, int *parity);
void lu_solve(int n, double A[][n], int *p_idx, double b[]);
double lu_det(int n, double A[][n], int parity);

/* newton_linesearch.c */
int newton_linesearch(int n, double x[], int *check,
                      void (*vecfuncP)(int, double [], double [], void *par),
                      void *par, int vecfuncP_ilow, int itmax, double tolf);


///* matrix_inv.c */
//int gaussjordan_inv(int n, double a[]);
//int M_to_Minv_gaussjordan(int n, const double M[], double Minv[]);

/* minbrent_brak.c */
int minbrent_brak(double (*f)(double x, void *p), void *par,
                  double ax, double bx, double cx,
                  int maxits, double tol,
                  double *xmin, double *fmin, int pr);
int min_brak(double (*func)(double x, void *p), void *par,
             double *ax, double *bx, double *cx,
             double *fa, double *fb, double *fc, int maxits);

/* WT_Random.c, BUT try to avoid RND */
double RND(void);

/* attenuation.c has Attenuation functions */
double Attenuation01(double x, double s, double p);
double Att_0to1_x_x0_w_q_s(double x, double x0, double w, double q, double s);
double dAtt_0to1_x_x0_w_q_s(double x, double x0, double w, double q, double s);
double Att_1to0_x_x0_w_q_s(double x, double x0, double w, double q, double s);
double dAtt_1to0_x_x0_w_q_s(double x, double x0, double w, double q, double s);


/* GSL_odeint.c */
double odeintegrate(double y[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])),
	int kmax, int *kcount, double *xp, double **yp, double dxsav,
	int *status);

/* NumericUtils_shims.c */
int newton_linesrch_itsP(double x[], int n, int *check,
                     void (*vecfuncP)(int, double [], double [], void *par),
                     void *par, int MAXITS, double TOLF);
int newton_linesrch_its(double x[], int n, int *check,
			void (*vecfunc)(int, double [], double []),
			int MAXITS, double TOLF);
int WT_newton(double *x, int n, int *check,
        void (*F_x)(int, double *x, double *Fx, void *par),
        void (*J_x_dx)(int, double *dx, double *Jdx, double *x, void *par),
        void *par, int MAXITS, double TOLF,
        int (*linSol)(int n, double *b, double *dx,
            void (*Jdx)(int, double *, double *, double *, void *),
            int (*precon)(int n, double *b, double *dx, double *x, void *par),
            double *x, void *par, int itmax, double tol),
        int (*precon)(int n, double *b, double *dx, double *x, void *par),
        int linitmax, double lintolfac);
int zbrac_P(double (*func)(double,void *par), double *x1, double *x2,
            void *par);
int zbrent_itsP(double *x0, double (*func)(double,void *par),
                double x1, double x2, void *par, int ITMAX, double tol);
double brent_with_pointer_to_pars(double ax, double bx, double cx,
             double (*f)(double, void *ppointer), double tol,
	     double *xmin, void *parpointer);
void mnbrak_with_pointer_to_pars(double *ax, double *bx, double *cx,
	    double *fa, double *fb, double *fc,
	    double (*func)(double, void *ppointer), void *parpointer);
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double [], double []));
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
