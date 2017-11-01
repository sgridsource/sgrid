/* wrappers_for_templates.h */
/* Wolfgang Tichy 8/2007 */

/* headers for functions called in wrappers_for_templates.c */
      
/* header for GMRES in templates */
int gmres_(long int *N, double *B, double *X, long int *RESTRT, 
            double *WORK, long int *LDW, double *H, long int *LDH,
            long int *ITER, double *RESID,
            int (*matvec)(double *alpha, double *x, double *beta, double *y),
            int (*psolve)(double *x, double *b),
            long int *INFO);
int bicgstab_(long int *N, double *B, double *X, double *WORK, long int *LDW, 
              long int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*psolve)(double *x, double *b),
              long int *INFO);
int cgs_(long int *N, double *B, double *X, double *WORK, long int *LDW, 
              long int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*psolve)(double *x, double *b),
              long int *INFO);
int qmr_(long int *N, double *B, double *X, double *WORK, long int *LDW, 
              long int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*matvectrans)(double *alpha, double *x, double *beta, double *y),
              int (*psolveQ)(double *x, double *b, char *s, short *slen),
              int (*psolveQtrans)(double *x, double *b, char *s, short *slen),
              long int *INFO);
int bicg_(long int *N, double *B, double *X, double *WORK, long int *LDW, 
              long int *ITER, double *RESID,
              int (*matvec)(double *alpha, double *x, double *beta, double *y),
              int (*matvectrans)(double *alpha, double *x, double *beta, double *y),
              int (*psolve)(double *x, double *b),
              int (*psolvetrans)(double *x, double *b),
              long int *INFO);
