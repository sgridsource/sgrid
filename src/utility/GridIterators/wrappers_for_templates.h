/* wrappers_for_templates.h */
/* Wolfgang Tichy 8/2007 */

/* headers for functions called in wrappers_for_templates.c */
      
/* header for GMRES in templates */
int gmres_(int *N, double *B, double *X, int *RESTRT, 
            double *WORK, int *LDW, double *H, int *LDH,
            int *ITER, double *RESID,
            int (*matvec)(double *alpha, double *x, double *beta, double *y),
            int (*psolve)(double *x, double *b),
            int *INFO);
int bicgstab_(int *N, double *B, double *X, double *WORK, int *LDW, 
              int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*psolve)(double *x, double *b),
              int *INFO);
int cgs_(int *N, double *B, double *X, double *WORK, int *LDW, 
              int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*psolve)(double *x, double *b),
              int *INFO);
