/* derivs.c */
/* Wolfgang Tichy 3/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* 1d Cheb deriv d1u[j] of u[j] defined on [a,b] and j=0, ... , n1-1 */
void cheb_d1(double a, double b, double *u, double *d1u, int n1)
{
  static int linelen=0;
  static double *c=NULL;
  static double *cp=NULL;
    
  /* get memory for Cheb coeffs c[j] and Cheb. coeffs cp[j] of deriv */
  /* static memory for lines */
  if(n1>linelen)
  {
    linelen = n1;
    c = (double*) realloc(c, linelen * sizeof(double));
    cp = (double*) realloc(cp, linelen * sizeof(double));
  }
  if(c==NULL) errorexit("cheb_d1: out of memory for c");
  if(cp==NULL) errorexit("cheb_d1: out of memory for cp");

  cheb_coeffs_fromExtrema(c, u, n1-1);
  cheb_deriv(a, b,  c , cp, n1-1);

  //cheb_filter(cp, 2*(n1-1)/3, n1-1);
  cheb_eval_onExtrema(cp, d1u, n1-1);

  /*
    linelen = 0;
    free(c);
    free(cp);
    c = cp = NULL;
  */
}

/* 1d Cheb derivs d1u[j] and  d2u[j] of u[j] defined on [a,b] 
   and j=0, ... , n1-1 */
void cheb_d2(double a, double b, double *u, double *d1u, double *d2u, int n1)
{
  static int linelen=0;
  static double *c=NULL;
  static double *cp=NULL;
    
  /* get memory for Cheb coeffs */
  /* static memory for lines */
  if(n1>linelen)
  {
    linelen = n1;
    c = (double*) realloc(c, linelen * sizeof(double));
    cp = (double*) realloc(cp, linelen * sizeof(double));
  }
  if(c==NULL) errorexit("cheb_d1: out of memory for c");
  if(cp==NULL) errorexit("cheb_d1: out of memory for cp");

//printf("2. c=%p cp=%p ",c,cp); yo();
  
  cheb_coeffs_fromExtrema(c, u, n1-1);

  cheb_deriv(a, b,  c , cp, n1-1);
  //cheb_filter(cp, 2*(n1-1)/3, n1-1);
  cheb_eval_onExtrema(cp, d1u, n1-1);

  cheb_deriv(a, b,  cp, c , n1-1);
  //cheb_filter(c, 2*(n1-1)/3, n1-1);
  cheb_eval_onExtrema(c, d2u, n1-1);

  /*
    linelen = 0;
    free(c);
    free(cp);
    c = cp = NULL;
  */
}

/* read a line of length nl 
   out of 3d var u in direction direc at pos i1,i2 */
void get_memline(double *u, double *line, int direc, int i1, int i2, 
                 int n1, int n2, int n3)
{
  int i;
  if(direc==1)
  {
    for(i=0; i<n1; i++)
      line[i] = u[Index(i,i1,i2)];
  }
  else if(direc==2)
  {
    for(i=0; i<n2; i++)
      line[i] = u[Index(i1,i,i2)];
  }
  else if(direc==3)
  {
    for(i=0; i<n3; i++)
      line[i] = u[Index(i1,i2,i)];
  }
  else
   errorexit("get_memline: possible values for direction direc are 1,2,3.");
}

/* put a line of length nl
   into 3d var u in direction direc at pos i1,i2 */
void put_memline(double *u, double *line, int direc, int i1, int i2,
                 int n1, int n2, int n3)
{
  int i;
  if(direc==1)
  {
    for(i=0; i<n1; i++)
      u[Index(i,i1,i2)] = line[i];
  }
  else if(direc==2)
  {
    for(i=0; i<n2; i++)
      u[Index(i1,i,i2)] = line[i];
  }
  else if(direc==3)
  {
    for(i=0; i<n3; i++)
      u[Index(i1,i2,i)] = line[i];
  }
  else
   errorexit("put_memline: possible values for direction direc are 1,2,3.");
}

/* compute first deriv of 3d var u in dirction direc on a box*/
void cheb_Deriv1(tBox *box, int direc, double *u, double *du)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *duline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    duline = (double*) realloc(duline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[0], box->bbox[1], uline, duline, box->n1);
        put_memline(du, duline, 1, j, k, box->n1, box->n2, box->n3);        
        */
        int n1=box->n1;
        int n2=box->n2;
        cheb_d1(box->bbox[0], box->bbox[1], 
                u+Index(0,j,k), du+Index(0,j,k), box->n1);
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[2], box->bbox[3], uline, duline, box->n2);
        put_memline(du, duline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[4], box->bbox[5], uline, duline, box->n3);
        put_memline(du, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("cheb_Deriv1: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    uline = duline = NULL;
  */
}

/* compute first deriv of 3d var u in dirction direc on a box*/
void cheb_allDerivs(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 )
{
  static int linelen=0;
  static double *uline=NULL;
  static double *duline=NULL;
  static double *dduline=NULL;
  int i,j,k, m3;

  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    duline = (double*) realloc(duline, linelen * sizeof(double));
    dduline = (double*) realloc(dduline, linelen * sizeof(double));
  }

  /* x-direction */
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        cheb_d2(box->bbox[0], box->bbox[1], uline, duline, dduline, box->n1);
        put_memline(u1, duline, 1, j, k, box->n1, box->n2, box->n3);
        put_memline(u11, dduline, 1, j, k, box->n1, box->n2, box->n3); 
        */
        int n1=box->n1;
        int n2=box->n2;
        cheb_d2(box->bbox[0], box->bbox[1], 
                u+Index(0,j,k), u1+Index(0,j,k), u11+Index(0,j,k), box->n1);
      }
  }
  
  /* y-direction */
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        cheb_d2(box->bbox[2], box->bbox[3], uline, duline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u1, uline, 2, i, k, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[2], box->bbox[3], uline, duline, box->n2);
        put_memline(u12, duline, 2, i, k, box->n1, box->n2, box->n3);
      }
  }
  
  /* z-direction */
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        cheb_d2(box->bbox[4], box->bbox[5], uline, duline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u1, uline, 3, i, j, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[4], box->bbox[5], uline, duline, box->n3);
        put_memline(u13, duline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u2, uline, 3, i, j, box->n1, box->n2, box->n3);
        cheb_d1(box->bbox[4], box->bbox[5], uline, duline, box->n3);
        put_memline(u23, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    free(dduline);
    uline = duline = dduline = NULL;
  */
}



/* Derivs with diff. matrix */

/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_Deriv1(tBox *box, int direc, double *u, double *du)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *duline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    duline = (double*) realloc(duline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D1, uline, duline, box->n1);
        put_memline(du, duline, 1, j, k, box->n1, box->n2, box->n3);        
        */
        int n1=box->n1;
        int n2=box->n2;
        diffmat_deriv(box->D1, u+Index(0,j,k), du+Index(0,j,k), n1);
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D2, uline, duline, box->n2);
        put_memline(du, duline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D3 , uline, duline, box->n3);
        put_memline(du, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("cheb_Deriv1: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    uline = duline = NULL;
  */
}


/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_allDerivs(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 )
{
  static int linelen=0;
  static double *uline=NULL;
  static double *duline=NULL;
  static double *dduline=NULL;
  int i,j,k, m3;

  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    duline = (double*) realloc(duline, linelen * sizeof(double));
    dduline = (double*) realloc(dduline, linelen * sizeof(double));
  }

  /* x-direction */
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D1, uline, duline, box->n1);
        diffmat_deriv(box->DD1, uline, dduline, box->n1);
        put_memline(u1, duline, 1, j, k, box->n1, box->n2, box->n3);
        put_memline(u11, dduline, 1, j, k, box->n1, box->n2, box->n3); 
        */
        int n1=box->n1;
        int n2=box->n2;
        diffmat_deriv(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        diffmat_deriv(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);
      }
  }
  
  /* y-direction */
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D2, uline, duline, box->n2);
        diffmat_deriv(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u1, uline, 2, i, k, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D2, uline, duline, box->n2);
        put_memline(u12, duline, 2, i, k, box->n1, box->n2, box->n3);
      }
  }
  
  /* z-direction */
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D3, uline, duline, box->n3);
        diffmat_deriv(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u1, uline, 3, i, j, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D3, uline, duline, box->n3);
        put_memline(u13, duline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u2, uline, 3, i, j, box->n1, box->n2, box->n3);
        diffmat_deriv(box->D3, uline, duline, box->n3);
        put_memline(u23, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    free(dduline);
    uline = duline = dduline = NULL;
  */
}
