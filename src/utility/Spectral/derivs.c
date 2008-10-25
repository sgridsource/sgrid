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


/*****************************************************************/
/* Derivs with diff. matrix                                      */
/*****************************************************************/

/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_Deriv1_diffmatrix(tBox *box, int direc, double *u, double *du)
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
        matrix_times_vector(box->D1, uline, duline, box->n1);
        put_memline(du, duline, 1, j, k, box->n1, box->n2, box->n3);        
        */
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), du+Index(0,j,k), n1);
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(du, duline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3 , uline, duline, box->n3);
        put_memline(du, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_Deriv1: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    uline = duline = NULL;
  */
}


/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_allDerivs_diffmatrix(tBox *box, double *u,
                    double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 )
{
  static int linelen=0;
  static double *uline=NULL;
  static double *duline=NULL;
  static double *dduline=NULL;
  int i,j,k, m3;
  static int second_deriv_order;

  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    duline = (double*) realloc(duline, linelen * sizeof(double));
    dduline = (double*) realloc(dduline, linelen * sizeof(double));
    
    /* read par which determines order of second derivs once */
    second_deriv_order = Geti("Spectral_second_deriv_order");
  }

  /* if X frist, then Y and then Z derivs are taken */
  if(second_deriv_order == 123)
  {
    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D1, uline, duline, box->n1);
        matrix_times_vector(box->DD1, uline, dduline, box->n1);
        put_memline(u1, duline, 1, j, k, box->n1, box->n2, box->n3);
        put_memline(u11, dduline, 1, j, k, box->n1, box->n2, box->n3); 
        */
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);
      }

    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u1, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u12, duline, 2, i, k, box->n1, box->n2, box->n3);
      }

    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u1, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u13, duline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u2, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u23, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else if(second_deriv_order == 132)
  {
    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D1, uline, duline, box->n1);
        matrix_times_vector(box->DD1, uline, dduline, box->n1);
        put_memline(u1, duline, 1, j, k, box->n1, box->n2, box->n3);
        put_memline(u11, dduline, 1, j, k, box->n1, box->n2, box->n3); 
        */
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);
      }

    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u1, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u13, duline, 3, i, j, box->n1, box->n2, box->n3);
      }

    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u1, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u12, duline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u3, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u23, duline, 2, i, k, box->n1, box->n2, box->n3);
      }

  }
  else if(second_deriv_order == 213)
  {
    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);
      }

    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
      }

    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u1, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u13, duline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u2, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u23, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else if(second_deriv_order == 231)
  {
    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);
      }

    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
        
        get_memline(u2, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        put_memline(u23, duline, 3, i, j, box->n1, box->n2, box->n3);
      }

    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }


  }
  else if(second_deriv_order == 312)
  {
    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
      }

    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }

    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u1, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u12, duline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u3, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u23, duline, 2, i, k, box->n1, box->n2, box->n3);
      }
  }
  else if(second_deriv_order == 321)
  {
    /* z-direction */
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D3, uline, duline, box->n3);
        matrix_times_vector(box->DD3, uline, dduline, box->n3);
        put_memline(u3, duline, 3, i, j, box->n1, box->n2, box->n3);
        put_memline(u33, dduline, 3, i, j, box->n1, box->n2, box->n3);
      }

    /* y-direction */
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        matrix_times_vector(box->DD2, uline, dduline, box->n2);
        put_memline(u2, duline, 2, i, k, box->n1, box->n2, box->n3);
        put_memline(u22, dduline, 2, i, k, box->n1, box->n2, box->n3);

        get_memline(u3, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->D2, uline, duline, box->n2);
        put_memline(u23, duline, 2, i, k, box->n1, box->n2, box->n3);
      }

    /* x-direction */
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }
  }
  else
    errorexit("spec_allDerivs: Spectral_second_deriv_order must be "
              "123, 132, 213, 231, 312 or 321");
  
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    free(dduline);
    uline = duline = dduline = NULL;
  */
}


/* compute 2nd deriv of 3d var u in dirction direc on a box */
void spec_Deriv2_diffmatrix(tBox *box, int direc, double *u, double *du)
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
        matrix_times_vector(box->DD1, uline, duline, box->n1);
        put_memline(du, duline, 1, j, k, box->n1, box->n2, box->n3);        
        */
        int n1=box->n1;
        int n2=box->n2;
        matrix_times_vector(box->DD1, u+Index(0,j,k), du+Index(0,j,k), n1);
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(box->DD2, uline, duline, box->n2);
        put_memline(du, duline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(box->DD3 , uline, duline, box->n3);
        put_memline(du, duline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_Deriv2: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  /*
    linelen=0;
    free(uline);
    free(duline);
    uline = duline = NULL;
  */
}


/*****************************************************************/
/* Derivs with FFTs                                              */
/*****************************************************************/

/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_Deriv1_FFT(tBox *box, int direc, double *u, double *du)
{
  double *uline;
  double *duline;
  int i,j,k, m3;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  /* memory for lines */
  if(direc==1)      m3=n1;
  else if(direc==2) m3=n2;
  else              m3=n3;
  uline = (double*)  calloc(m3, sizeof(double));
  duline = (double*) calloc(m3, sizeof(double));

  if(direc==1)
  {
    double a = box->bbox[0];
    double b = box->bbox[1];
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        box->get_coeffs1(duline, uline, n1-1);
        box->coeffs_of_deriv1(a,b, duline, uline, n1-1);
        box->eval_onPoints1(uline, duline, n1-1);
        put_memline(du, duline, 1, j, k, n1, n2, n3);        
      }
  }
  else if(direc==2)
  {
    double a = box->bbox[2];
    double b = box->bbox[3];
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        box->get_coeffs2(duline, uline, n2-1);
        box->coeffs_of_deriv2(a,b, duline, uline, n2-1);
        box->eval_onPoints2(uline, duline, n2-1);
        put_memline(du, duline, 2, i, k, n1, n2, n3);        
      }
  }
  else if(direc==3)
  {
    double a = box->bbox[4];
    double b = box->bbox[5];
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        box->get_coeffs3(duline, uline, n3-1);
        box->coeffs_of_deriv3(a,b, duline, uline, n3-1);
        box->eval_onPoints3(uline, duline, n3-1);
        put_memline(du, duline, 3, i, j, n1, n2, n3);
      }
  }
  else
    errorexit("spec_Deriv1_FFT: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  free(uline);
  free(duline);
}


/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_allDerivs_FFT(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 )
{
  int second_deriv_order = Geti("Spectral_second_deriv_order");
  double *uline;
  double *duline;
  double *dduline;
  int i,j,k, m3;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double a,b;

  /* memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  uline = (double*)  calloc(m3, sizeof(double));
  duline = (double*) calloc(m3, sizeof(double));
  dduline = (double*) calloc(m3, sizeof(double));
  
  /* if X frist, then Y and then Z derivs are taken */
  if(second_deriv_order == 123)
  {
    /* x-direction */
    a = box->bbox[0];
    b = box->bbox[1];
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        box->get_coeffs1(duline, uline, n1-1);
        box->coeffs_of_deriv1(a,b, duline, uline, n1-1); /* coeffs of du are now in uline */
        box->eval_onPoints1(uline, duline, n1-1);
        put_memline(u1, duline, 1, j, k, n1, n2, n3);

        box->coeffs_of_deriv1(a,b, uline, duline, n1-1); /* coeffs of ddu are now in duline */
        box->eval_onPoints1(duline, dduline, n1-1);
        put_memline(u11, dduline, 1, j, k, n1, n2, n3); 
      }

    /* y-direction */
    a = box->bbox[2];
    b = box->bbox[3];
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);

        get_memline(u1, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u12, duline, 2, i, k, n1, n2, n3);
      }

    /* z-direction */
    a = box->bbox[4];
    b = box->bbox[5];
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
        
        get_memline(u1, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u13, duline, 3, i, j, n1, n2, n3);
        
        get_memline(u2, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u23, duline, 3, i, j, n1, n2, n3);
      }
  }
  else if(second_deriv_order == 132)
  {
    /* x-direction */
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        matrix_times_vector(box->D1, uline, duline, n1);
        matrix_times_vector(box->DD1, uline, dduline, n1);
        put_memline(u1, duline, 1, j, k, n1, n2, n3);
        put_memline(u11, dduline, 1, j, k, n1, n2, n3); 
        */
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);
      }

    /* z-direction */
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
        
        get_memline(u1, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u13, duline, 3, i, j, n1, n2, n3);
      }

    /* y-direction */
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);

        get_memline(u1, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u12, duline, 2, i, k, n1, n2, n3);

        get_memline(u3, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u23, duline, 2, i, k, n1, n2, n3);
      }

  }
  else if(second_deriv_order == 213)
  {
    /* y-direction */
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);
      }

    /* x-direction */
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
      }

    /* z-direction */
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
        
        get_memline(u1, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u13, duline, 3, i, j, n1, n2, n3);
        
        get_memline(u2, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u23, duline, 3, i, j, n1, n2, n3);
      }
  }
  else if(second_deriv_order == 231)
  {
    /* y-direction */
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);
      }

    /* z-direction */
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
        
        get_memline(u2, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        put_memline(u23, duline, 3, i, j, n1, n2, n3);
      }

    /* x-direction */
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }


  }
  else if(second_deriv_order == 312)
  {
    /* z-direction */
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
      }

    /* x-direction */
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }

    /* y-direction */
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);

        get_memline(u1, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u12, duline, 2, i, k, n1, n2, n3);

        get_memline(u3, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u23, duline, 2, i, k, n1, n2, n3);
      }
  }
  else if(second_deriv_order == 321)
  {
    /* z-direction */
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->D3, uline, duline, n3);
        matrix_times_vector(box->DD3, uline, dduline, n3);
        put_memline(u3, duline, 3, i, j, n1, n2, n3);
        put_memline(u33, dduline, 3, i, j, n1, n2, n3);
      }

    /* y-direction */
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        matrix_times_vector(box->DD2, uline, dduline, n2);
        put_memline(u2, duline, 2, i, k, n1, n2, n3);
        put_memline(u22, dduline, 2, i, k, n1, n2, n3);

        get_memline(u3, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->D2, uline, duline, n2);
        put_memline(u23, duline, 2, i, k, n1, n2, n3);
      }

    /* x-direction */
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        matrix_times_vector(box->D1, u+Index(0,j,k), u1+Index(0,j,k), n1);
        matrix_times_vector(box->DD1, u+Index(0,j,k), u11+Index(0,j,k), n1);

        matrix_times_vector(box->D1, u2+Index(0,j,k), u12+Index(0,j,k), n1);
        matrix_times_vector(box->D1, u3+Index(0,j,k), u13+Index(0,j,k), n1);
      }
  }
  else
    errorexit("spec_allDerivs_FFT: Spectral_second_deriv_order must be "
              "123, 132, 213, 231, 312 or 321");
  
  /* free memory for lines */
  free(uline);
  free(duline);
  free(dduline);
}


/* compute 2nd deriv of 3d var u in dirction direc on a box */
void spec_Deriv2_FFT(tBox *box, int direc, double *u, double *du)
{
  double *uline;
  double *duline;
  int i,j,k, m3;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  /* memory for lines */
  if(direc==1)      m3=n1;
  else if(direc==2) m3=n2;
  else              m3=n3;
  uline = (double*)  calloc(m3, sizeof(double));
  duline = (double*) calloc(m3, sizeof(double));

  if(direc==1)
  {
    double a = box->bbox[0];
    double b = box->bbox[1];
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        box->get_coeffs1(duline, uline, n1-1);
        box->coeffs_of_deriv1(a,b, duline, uline, n1-1); /* coeffs of du are now in uline */
        box->coeffs_of_deriv1(a,b, uline, duline, n1-1); /* coeffs of ddu are now in duline */
        box->eval_onPoints1(duline, uline, n1-1); /* ddu is in now in uline */
        put_memline(du, uline, 1, j, k, n1, n2, n3);
      }
  }
  else if(direc==2)
  {
    double a = box->bbox[2];
    double b = box->bbox[3];
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        box->get_coeffs2(duline, uline, n2-1);
        box->coeffs_of_deriv2(a,b, duline, uline, n2-1); /* coeffs of du are now in uline */
        box->coeffs_of_deriv2(a,b, uline, duline, n2-1); /* coeffs of ddu are now in duline */
        box->eval_onPoints2(duline, uline, n2-1); /* ddu is now in uline */
        put_memline(du, duline, 2, i, k, n1, n2, n3);        
      }
  }
  else if(direc==3)
  {
    double a = box->bbox[4];
    double b = box->bbox[5];
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        box->get_coeffs2(duline, uline, n2-1);
        box->coeffs_of_deriv2(a,b, duline, uline, n2-1); /* coeffs of du are now in uline */
        box->coeffs_of_deriv2(a,b, uline, duline, n2-1); /* coeffs of ddu are now in duline */
        box->eval_onPoints2(duline, uline, n2-1); /* ddu is now in uline */
        put_memline(du, duline, 3, i, j, n1, n2, n3);
      }
  }
  else
    errorexit("spec_Deriv2_FFT: possible values for direction direc are 1,2,3.");
  /* free memory for lines */
  free(uline);
  free(duline);
}


/*****************************************************************/
/* pick if we use Derivs with or without FFT                     */
/*****************************************************************/

/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_Deriv1(tBox *box, int direc, double *u, double *du)
{
  spec_Deriv1_diffmatrix(box, direc, u, du);
}

/* compute first deriv of 3d var u in dirction direc on a box*/
void spec_allDerivs(tBox *box, double *u, double *u1, double *u2, double *u3,
                    double *u11,double *u12,double *u13,
                    double *u22,double *u23,double *u33 )
{
  spec_allDerivs_diffmatrix(box, u, u1,u2,u3, u11,u12,u13,u22,u23,u33);
}

/* compute 2nd deriv of 3d var u in dirction direc on a box */
void spec_Deriv2(tBox *box, int direc, double *u, double *du)
{
  spec_Deriv2_diffmatrix(box, direc, u, du);
}
