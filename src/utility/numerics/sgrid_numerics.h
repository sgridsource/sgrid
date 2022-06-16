/* nmesh_numerics.h */
/* Wolfgang Tichy, June 2019 */


/* newton1d_brak.c */
int newton1d_brak(double *x0,
                  void (*fdf)(double x, void *par, double *f, double *df),
                  double x1, double x2, void *par, int maxits, double xacc,
                  int pr);
int rtbisect(double *x0, double (*func)(double,void *par),
             double x1, double x2, void *par, int maxits, double xacc, int pr);

/* rt_brak.c */
int rt_brak(double (*func)(double x, void *par),
            double *x1, double *x2, void *par);

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

/* NumericUtils_shims.c */
int newton_linesrch_itsP(double x[], int n, int *check,
                     void (*vecfuncP)(int, double [], double [], void *par),
                     void *par, int MAXITS, double TOLF);
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
#define zbrac_P rt_brak
int zbrent_itsP(double *x0, double (*func)(double,void *par),
                double x1, double x2, void *par, int ITMAX, double tol);

/* WT_Random.c, BUT try to avoid RND */
double RND(void);
