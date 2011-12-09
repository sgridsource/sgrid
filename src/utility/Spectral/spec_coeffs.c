/* spec_coeffs.c */
/* Wolfgang Tichy 7/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* init a n1*n1 matrix M used to compute coeffs */
void initMatrix_ForCoeffs(double *M, int n1,
                          void (*get_coeffs)(double *,double *, int))
{
  int i,j;
  double *u;
  double *c;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));

  if( !(u && c) ) errorexit("initMatrix_ForCoeffs: out of memory for u, c");


  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    u[j]=1.0;

    get_coeffs(c, u, n1-1);
    
    /* set M */
    for(i=0; i<n1; i++) M[n1*i + j] = c[i];
    
    u[j]=0.0;
  }
  free(u);
  free(c);
}

/* init a n1*n1 matrix M used to compute function u from coeffs */
void initMatrix_ToEvaluate(double *M, int n1,
                           void (*eval_onPoints)(double *,double *, int))
{
  int i,j;
  double *u;
  double *c;

  u = (double *) calloc(n1, sizeof(double));
  c = (double *) calloc(n1, sizeof(double));

  if( !(u && c) ) errorexit("initMatrix_ToEvaluate: out of memory for u, c");


  /* read matrix from functions */
  for(j=0; j<n1; j++)
  {
    c[j]=1.0;

    eval_onPoints(c, u, n1-1);
    
    /* set M */
    for(i=0; i<n1; i++) M[n1*i + j] = u[i];
    
    c[j]=0.0;
  }
  free(u);
  free(c);
}

/**************************************************/
/* analysis and synthesis using diffmatrices only */
/**************************************************/
/* get the coeffs c of 3d var u in direction direc on the whole box */
void spec_analysis1_diffmatrix(tBox *box, int direc, double *M, double *u, double *c)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *cline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    cline = (double*) realloc(cline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        get_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n1);
        put_memline(c, cline, 1, j, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n2);
        put_memline(c, cline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(M, uline, cline, box->n3);
        put_memline(c, cline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_analysis1: possible values for direction direc are 1,2,3.");
}

/* compute the 3d var u from the coeffs c in direction direc  
   on the whole box                                         */
void spec_synthesis1_diffmatrix(tBox *box, int direc, double *M, double *u, double *c)
{
  static int linelen=0;
  static double *uline=NULL;
  static double *cline=NULL;
  int i,j,k, m3;
    
  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    cline = (double*) realloc(cline, linelen * sizeof(double));
  }

  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        get_memline(c, cline, 1, j, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n1);
        put_memline(u, uline, 1, j, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(c, cline, 2, i, k, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n2);
        put_memline(u, uline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(c, cline, 3, i, j, box->n1, box->n2, box->n3);
        matrix_times_vector(M, cline, uline, box->n3);
        put_memline(u, uline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_synthesis1: possible values for direction direc are 1,2,3.");
}

/********************************************************************/
/* analysis and synthesis which may use FFTs                        */
/********************************************************************/
/* get the coeffs c of 3d var u in direction direc on the whole box */
void spec_analysis1(tBox *box, int direc, double *u, double *c)
{
  double *M;
  void (*get_coeffs)(double *,double *, int);
  int j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  /* M and get_coeffs */
  if(direc==1)      { M=box->Mcoeffs1; get_coeffs=box->get_coeffs1; }
  else if(direc==2) { M=box->Mcoeffs2; get_coeffs=box->get_coeffs2; }
  else              { M=box->Mcoeffs3; get_coeffs=box->get_coeffs3; }

  if(direc==1)
  {
    if(box->TransformType1)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(j,k, n2,n3, uline,cline,n1)
      {
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        get_coeffs(cline, uline, n1-1);
        put_memline(c, cline, 1, j, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(j,k, n2,n3, uline,cline,n1)
      {
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        matrix_times_vector(M, uline, cline, n1);
        put_memline(c, cline, 1, j, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else if(direc==2)
  {
    if(box->TransformType2)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,k, n1,n3, uline,cline,n2)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        get_coeffs(cline, uline, n2-1);
        put_memline(c, cline, 2, i, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,k, n1,n3, uline,cline,n2)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(M, uline, cline, n2);
        put_memline(c, cline, 2, i, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else if(direc==3)
  {
    if(box->TransformType3)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,j, n1,n2, uline,cline,n3)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        get_coeffs(cline, uline, n3-1);
        put_memline(c, cline, 3, i, j, n1, n2, n3);
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,j, n1,n2, uline,cline,n3)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(M, uline, cline, n3);
        put_memline(c, cline, 3, i, j, n1, n2, n3);
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else
    errorexit("spec_analysis1: possible values for direction direc are 1,2,3.");
}

/* compute the 3d var u from the coeffs c in direction direc  
   on the whole box                                         */
void spec_synthesis1(tBox *box, int direc, double *u, double *c)
{
  double *M;
  void (*eval_onPoints)(double *,double *, int);
  int j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  /* M and eval_onPoints */
  if(direc==1)      { M=box->Meval1; eval_onPoints=box->eval_onPoints1; }
  else if(direc==2) { M=box->Meval2; eval_onPoints=box->eval_onPoints2; }
  else              { M=box->Meval3; eval_onPoints=box->eval_onPoints3; }

  if(direc==1)
  {
    if(box->TransformType1)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(j,k, n2,n3, uline,cline,n1)
      {
        get_memline(c, cline, 1, j, k, n1, n2, n3);
        eval_onPoints(cline, uline, n1-1);
        put_memline(u, uline, 1, j, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(j,k, n2,n3, uline,cline,n1)
      {
        get_memline(c, cline, 1, j, k, n1, n2, n3);
        matrix_times_vector(M, cline, uline, n1);
        put_memline(u, uline, 1, j, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else if(direc==2)
  {
    if(box->TransformType2)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,k, n1,n3, uline,cline,n2)
      {
        get_memline(c, cline, 2, i, k, n1, n2, n3);
        eval_onPoints(cline, uline, n2-1);
        put_memline(u, uline, 2, i, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,k, n1,n3, uline,cline,n2)
      {
        get_memline(c, cline, 2, i, k, n1, n2, n3);
        matrix_times_vector(M, cline, uline, n2);
        put_memline(u, uline, 2, i, k, n1, n2, n3);        
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else if(direc==3)
  {
    if(box->TransformType3)
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,j, n1,n2, uline,cline,n3)
      {
        get_memline(c, cline, 3, i, j, n1, n2, n3);
        eval_onPoints(cline, uline, n3-1);
        put_memline(u, uline, 3, i, j, n1, n2, n3);
      } end_forLines_free2Lines(uline,cline)
    }
    else
    {
      SGRID_LEVEL3_Pragma(omp parallel for)
      forLines_alloc2Lines(i,j, n1,n2, uline,cline,n3)
      {
        get_memline(c, cline, 3, i, j, n1, n2, n3);
        matrix_times_vector(M, cline, uline, n3);
        put_memline(u, uline, 3, i, j, n1, n2, n3);
      } end_forLines_free2Lines(uline,cline)
    }
  }
  else
    errorexit("spec_synthesis1: possible values for direction direc are 1,2,3.");
}



/* check if N is a power of 2 */
int is_power_of_two(int N)
{
  unsigned int nn, ng;
  
  ng=0;
  nn=N;
  while(nn >>= 1) ng++;
  if(N == (1L << ng)) return 1; /* N is power of 2 */
  return 0;
}

/* set flags which determine if we use FFTs */
void set_TransformType_flags_inbox(tBox *box)
{
  char str[1000];
  char bas[1000];
  int dir;
  
  for(dir=1; dir<=3 ; dir++)
  {
    int nb;
    int *Ttype;

    nb=box->n1;  Ttype=&(box->TransformType1); /* <--for dir=1 */
    if(dir==2) { nb=box->n2; Ttype=&(box->TransformType2); }
    else if(dir==3) { nb=box->n3; Ttype=&(box->TransformType3); }

    snprintf(str, 999, "box%d_TransformType%d", box->b, dir);
    snprintf(bas, 999, "box%d_basis%d", box->b, dir);
    /* printf("str=%s=%s\nbas=%s=%s\n", str, Gets(str), bas, Gets(bas)); */
    if( Getv(str, "NUMREC_FFT") && nb>=Geti(str) && Geti(str)>=0 &&
        ( (Getv(bas, "ChebExtrema") && is_power_of_two(nb-1)) ||
          (Getv(bas, "ChebZeros") && is_power_of_two(nb)) ||
          (Getv(bas, "Fourier") && is_power_of_two(nb))          )  )
      *Ttype=NUMREC_FFT;
    else if( Getv(str, "FFTW3") && nb>=Geti(str) && Geti(str)>=0 &&
             ( Getv(bas, "ChebExtrema") ||
               Getv(bas, "ChebZeros") ||
               Getv(bas, "Fourier")        )  )
      *Ttype=FFTW3_FFT;
    else
      *Ttype=MATRIX_MULTIPLICATION;
  }
}

/* get all relevant function points in one box in direction direc */
/* this is the old void get_spec_functionpointers */
void get_spec_functionpointers_from_pars(tBox *box, int direc,
     void (**get_coeffs)(double *,double *, int),
     void (**coeffs_of_deriv)(double, double, double *,double *, int),
     void (**coeffs_of_2ndderiv)(double, double, double *,double *, int),
     void (**coeffs_of_int)(double, double, double *,double *, int),
     void (**eval_onPoints)(double *,double *, int),
     void (**filter_coeffs)(double *, int, int),
     double (**basisfunc)(void *aux, double a, double b, int k, int N, double X) )
{       
  char str[1000];
  int *Ttype;

  *get_coeffs=NULL;
  *coeffs_of_deriv=NULL;
  *coeffs_of_2ndderiv=NULL;
  *coeffs_of_int=NULL;
  *eval_onPoints=NULL;
  *filter_coeffs=NULL;
  *basisfunc=NULL;

  snprintf(str, 999, "box%d_basis%d", box->b, direc);
  //printf("%s=%s\n", str, Gets(str));
  Ttype=&(box->TransformType1); /* <--for direc=1 */
  if(direc==2)      { Ttype=&(box->TransformType2); }
  else if(direc==3) { Ttype=&(box->TransformType3); }

  if( Getv(str, "ChebExtrema") )
  {
    *get_coeffs = cheb_coeffs_fromExtrema;
    *coeffs_of_deriv = cheb_deriv;
    *coeffs_of_int = cheb_int;
    *eval_onPoints = cheb_eval_onExtrema;
    *filter_coeffs = cheb_filter;
    *basisfunc = cheb_basisfunc;
    if(*Ttype==NUMREC_FFT)
    {
      *get_coeffs = cheb_coeffs_fromExtrema_numrecFFT;
      *eval_onPoints = cheb_eval_onExtrema_numrecFFT;
    }
    else if(*Ttype==FFTW3_FFT)
    {
      *get_coeffs = cheb_coeffs_fromExtrema_FFTW3;
      *eval_onPoints = cheb_eval_onExtrema_FFTW3;
    }
  }
  else if( Getv(str, "Fourier") )
  {
    *get_coeffs = four_coeffs;
    *coeffs_of_deriv = four_deriv;
    *coeffs_of_int = four_int;
    *eval_onPoints = four_eval;
    *filter_coeffs = four_filter;
    *basisfunc = four_basisfunc;
    if(*Ttype==NUMREC_FFT)
    {
      *get_coeffs = four_coeffs_numrecFFT;
      *coeffs_of_int = four_int_numrecFFT;
      *eval_onPoints = four_eval_numrecFFT;
    }
    else if(*Ttype==FFTW3_FFT)
    {
      *get_coeffs = four_coeffs_FFTW3;
      *eval_onPoints = four_eval_FFTW3;
    }
  }
  else if( Getv(str, "fd2_onesided") )
  {
    *get_coeffs = fd2_coeffs;
    *coeffs_of_deriv = fd2_deriv_onesidedBC;
    *eval_onPoints = fd2_eval;
    *filter_coeffs = fd2_filter;
    *coeffs_of_2ndderiv = fd2_2ndderiv_onesidedBC;
  }
  else if( Getv(str, "fd2_periodic") )
  {
    *get_coeffs = fd2_coeffs;
    *coeffs_of_deriv = fd2_deriv_periodic;
    *eval_onPoints = fd2_eval;
    *filter_coeffs = fd2_filter;
    *coeffs_of_2ndderiv = fd2_2ndderiv_periodic;
  }
  else if( Getv(str, "ChebZeros") )
  {
    *get_coeffs = cheb_coeffs_fromZeros;
    *coeffs_of_deriv = cheb_deriv;
    *coeffs_of_int = cheb_int;
    *eval_onPoints = cheb_eval_onZeros;
    *filter_coeffs = cheb_filter;
    *basisfunc = cheb_basisfunc;
    if(*Ttype==NUMREC_FFT)
    {
      *get_coeffs = cheb_coeffs_fromZeros_numrecFFT;
      *eval_onPoints = cheb_eval_onZeros_numrecFFT;
    }
    else if(*Ttype==FFTW3_FFT)
    {
      *get_coeffs = cheb_coeffs_fromZeros_FFTW3;
      *eval_onPoints = cheb_eval_onZeros_FFTW3;
    }
  }
  else
  {
    printf("get_spec_functionpointers_from_pars: %s=%s is unknown!\n", str, Gets(str));
    errorexits("get_spec_functionpointers_from_pars: don't know what to do "
               "with %s" , Gets(str));
  }
}

/* get func pointer to get_coeffs */
void get_spec_functionpointerTO_get_coeffs(tBox *box, int direc,
                               void (**get_coeffs)(double *,double *, int))
{
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int n1, double X)=NULL;
           
  get_spec_functionpointers(box, direc, get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, 
                            &eval_onPoints, &filter_coeffs, &basisfunc);
}

/* get all relevant function pointers and write them into box struct */
void init_spec_functionpointers(tBox *box)
{       
  char str[1000];
  int dir;
  void (*get_coeffs)(double *,double *, int);
  void (*coeffs_of_deriv)(double, double, double *,double *, int);
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int);
  void (*coeffs_of_int)(double, double, double *,double *, int);
  void (*eval_onPoints)(double *,double *, int);
  void (*filter_coeffs)(double *, int, int);
  double (*basisfunc)(void *aux, double a, double b, int k, int N, double X);
                                  
  /* loop over directions */
  for(dir=1; dir<=3; dir++)
  {
    get_coeffs=NULL;
    coeffs_of_deriv=NULL;
    coeffs_of_2ndderiv=NULL;
    coeffs_of_int=NULL;
    eval_onPoints=NULL;
    filter_coeffs=NULL;
    basisfunc=NULL;

    /* get the func pointers in direction dir */
    get_spec_functionpointers_from_pars(box, dir, 
                              &get_coeffs, &coeffs_of_deriv,
                              &coeffs_of_2ndderiv, &coeffs_of_int,
                              &eval_onPoints, &filter_coeffs, &basisfunc);

    /* write into box struct */
    if(dir==1)
    {
      box->get_coeffs1      = get_coeffs;
      box->coeffs_of_deriv1 = coeffs_of_deriv;
      box->coeffs_of_int1   = coeffs_of_int;
      box->eval_onPoints1   = eval_onPoints;
      box->filter_coeffs1   = filter_coeffs;
      box->basis1           = basisfunc;
    }
    else if(dir==2)
    {
      box->get_coeffs2      = get_coeffs;
      box->coeffs_of_deriv2 = coeffs_of_deriv;
      box->coeffs_of_int2   = coeffs_of_int;
      box->eval_onPoints2   = eval_onPoints;
      box->filter_coeffs2   = filter_coeffs;
      box->basis2           = basisfunc;
    }
    else if(dir==3)
    {
      box->get_coeffs3      = get_coeffs;
      box->coeffs_of_deriv3 = coeffs_of_deriv;
      box->coeffs_of_int3   = coeffs_of_int;
      box->eval_onPoints3   = eval_onPoints;
      box->filter_coeffs3   = filter_coeffs;
      box->basis3           = basisfunc;
    }
    else
      errorexit("init_spec_functionpointers: dir must be 1,2 or 3");
  }
}

/* read all relevant function points in one box in direction dir from box struct */
void get_spec_functionpointers(tBox *box, int dir,
     void (**get_coeffs)(double *,double *, int),
     void (**coeffs_of_deriv)(double, double, double *,double *, int),
     void (**coeffs_of_2ndderiv)(double, double, double *,double *, int),
     void (**coeffs_of_int)(double, double, double *,double *, int),
     void (**eval_onPoints)(double *,double *, int),
     void (**filter_coeffs)(double *, int, int),
     double (**basisfunc)(void *aux, double a, double b, int k, int N, double X) )
{       
  /* read from  box struct */
  if(dir==1)
  {
    *get_coeffs      = box->get_coeffs1;
    *coeffs_of_deriv = box->coeffs_of_deriv1;
    *coeffs_of_int   = box->coeffs_of_int1;
    *eval_onPoints   = box->eval_onPoints1;
    *filter_coeffs   = box->filter_coeffs1;
    *basisfunc       = box->basis1;
  }
  else if(dir==2)
  {
    *get_coeffs      = box->get_coeffs2;
    *coeffs_of_deriv = box->coeffs_of_deriv2;
    *coeffs_of_int   = box->coeffs_of_int2;
    *eval_onPoints   = box->eval_onPoints2;
    *filter_coeffs   = box->filter_coeffs2;
    *basisfunc       = box->basis2;
  }
  else if(dir==3)
  {
    *get_coeffs      = box->get_coeffs3;
    *coeffs_of_deriv = box->coeffs_of_deriv3;
    *coeffs_of_int   = box->coeffs_of_int3;
    *eval_onPoints   = box->eval_onPoints3;
    *filter_coeffs   = box->filter_coeffs3;
    *basisfunc       = box->basis3;
  }
  else
    errorexit("get_spec_functionpointers: dir must be 1,2 or 3");
}


/* THIS IS OLD USE: spec_Basis_times_CoeffMatrix_direc INSTEAD!!! */
/* compute B_k(X) M_ki     <-- B_k is basis function k
   Cu_k = M_ki u_i         <-- M_ki is coeff matrix
   u(X) = Cu_k B_k(X) = B_k(X) M_ki u_i = BM_i u_i      */
/* 
USE something like this before calling spec_Basis_times_CoeffMatrix:
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, int n1, double X)=NULL;

  get_spec_functionpointers(box, direc, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  if(basisfunc==NULL)
  {
    printf("spec_Basis_times_CoeffMatrix: box%d, direc=%d\n", box->b, direc);
    errorexiti("spec_Basis_times_CoeffMatrix: "
               "basisfunc=NULL (direc=%d)", direc);
  }

  if(direc==1)      { n = box->n1; a=box->bbox[0]; b=box->bbox[1]; }
  else if(direc==2) { n = box->n2; a=box->bbox[2]; b=box->bbox[3]; }
  else if(direc==3) { n = box->n3; a=box->bbox[4]; b=box->bbox[5]; }
  else
    errorexit("spec_Basis_times_CoeffMatrix: possible values for direction direc are 1,2,3.");
*/
void spec_Basis_times_CoeffMatrix(double a, double b, int n,
                                  double *BM, double X,
                    void   (*get_coeffs)(double *,double *, int),
                    double (*basisfunc)(void *aux, double a, double b, int k, int n1, double X))
{			  /* basisfunc is something like cheb_basisfunc */
  double *M;
  double *B;
  int i;

  M = (double *) calloc(n*n, sizeof(double));
  B = (double *) calloc(n, sizeof(double));

  if( !(M && B) ) errorexit("spec_Basis_times_CoeffMatrix: out of memory for M, B");
    
  /* initialize the matrix M used to compute coeffs */
  initMatrix_ForCoeffs(M, n, get_coeffs);

  /* initialize basis functions at point X */
  for(i = 0; i < n; i++)  B[i] = basisfunc(NULL, a,b, i,n, X); // B[i]=B_i(X)
                                    /* NOTE: ^here it would be good to pass 
                                       in (void *) box instead of NULL, 
                                       but we don't have box in this func */

  /* Cu_k = M_ki u_i   <-- M is coeff matrix
     u(X) = Cu_k B_k(X) = B_k(X) M_ki u_i = BM_i u_i */

  /* get BM_i = B_k(X) M_ki */
  vector_times_matrix(B, M, BM, n);

  /* free memory for matrix M and basis funcs B */
  free(M);
  free(B);
}

/* like spec_Basis_times_CoeffMatrix, but we read all info from box */
void spec_Basis_times_CoeffMatrix_direc(tBox *box, int dir, 
                                        double *BM, double X)
{
  double *B;
  int i;
  int n=max3(box->n1, box->n2, box->n3);

  B = (double *) calloc(n, sizeof(double));
  if(!B) errorexit("spec_Basis_times_CoeffMatrix_inbox: out of memory for B");
    
  /* initialize basis functions at point X */
  /* Cu_k = M_ki u_i   <-- M is coeff matrix
     u(X) = Cu_k B_k(X) = B_k(X) M_ki u_i = BM_i u_i */
  if(dir==1)
  {
    n=box->n1;
    for(i=0; i<n; i++) // B[i]=B_i(X)
      B[i]=box->basis1((void *) box, box->bbox[0],box->bbox[1], i,n, X);
    vector_times_matrix(B, box->Mcoeffs1, BM, n); /* get BM_i = B_k(X) M_ki */
  }
  else if(dir==2)
  {
    n=box->n2;
    for(i=0; i<n; i++) // B[i]=B_i(X)
      B[i]=box->basis2((void *) box, box->bbox[2],box->bbox[3], i,n, X);
    vector_times_matrix(B, box->Mcoeffs2, BM, n); /* get BM_i = B_k(X) M_ki */
  }
  else if(dir==3)
  {
    n=box->n3;
    for(i=0; i<n; i++) // B[i]=B_i(X)
      B[i]=box->basis3((void *) box, box->bbox[4],box->bbox[5], i,n, X);
    vector_times_matrix(B, box->Mcoeffs3, BM, n); /* get BM_i = B_k(X) M_ki */
  }
  else
    errorexit("spec_Basis_times_CoeffMatrix_inbox: dir must be 1,2 or 3");

  /* free memory for basis funcs B */
  free(B);
}
