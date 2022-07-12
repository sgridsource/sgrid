/* explicit_Four_trafos.c */
/* does explicit slow Fourier trafos: */
/* Wolfgang Tichy 4/2004 */

#include "sgrid.h"
#include "Spectral.h"

#ifdef FFTW3
#include "fftw3.h"
#endif

/* define PI */
#define twoPI  6.28318530717958647692528676655901
#define PI     3.14159265358979323846264338327950
#define PIh    1.57079632679489661923132169163975
#define PIq    0.785398163397448309615660845819876

/* numrec funcs */
void numrec_four1(double data[], unsigned long nn, int isign);
void numrec_realft(double data[], unsigned long n, int isign);
void numrec_cosft1(double y[], int n);
void numrec_cosft2(double y[], int n, int isign);

/* some global FFTW3 plan arrays */
#ifdef FFTW3
int FFTW3_Nmax_of_plan_array;
fftw_plan *FFTW3_dft_r2c_1d_plan; /* for FT of real data */
fftw_plan *FFTW3_dft_c2r_1d_plan; /* for inverse FT of complex coeff to real */
fftw_plan *FFTW_REDFT00_1d_plan;  /* like cosft1 in numrec */
fftw_plan *FFTW_REDFT10_1d_plan;  /* like cosft2 in numrec */
fftw_plan *FFTW_REDFT01_1d_plan;  /* like inverse cosft2 in numrec */
#endif



/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 
NOTE: four_coeffs returns c[] that are N times of those of four_coeffs_alt */
void four_coeffs_numrecFFT(double c[], double u[], int n)
{
  int j;
  double temp;
  int N=n+1;

  for(j=0; j<N; j++)  c[j] = u[j];
  numrec_realft(c-1, N, 1);
  temp=c[1];
  for(j=2; j<N; j++)  c[j-1] = c[j];
  c[N-1]=temp;
}


/* find function u from Four coeffs c[0...n], computed with four_coeffs */
void four_eval_numrecFFT(double c[], double u[], int n)
{
  int j;
  int N=n+1;
  double toN = 2.0/N;

  u[0] = toN*c[0];
  u[1] = toN*c[N-1];
  for(j=2; j<N; j++)  u[j] = toN*c[j-1];
  numrec_realft(u-1, N, -1);
}


/* compute Four coeffs of integral cint[0...n] from Four coeffs c[0...n] */
void four_int_numrecFFT(double a, double b, double c[], double cint[], int n)
{
  int j;
  double PI2_con, L;
  int N=n+1;
  double *u = (double*) calloc(N+1, sizeof(double));

  L = b-a;
  PI2_con = 2.0*PI/L;

  /* get terms coming from integrating c[0] */
/*
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] = -0.5*L*c[0]/((double) N);
    cint[2*j]   = -0.5*L*c[0]/((double) N); // WRONG!!!!!
    // integrate the func 1 and find its coeffs instead!!!
    // multiply this by c[0]
    if(N!=4) errorexit("four_int is wrong");
  }
  if( N%2 == 0 ) cint[N-1] = -0.5*L*c[0]/((double) N);
  cint[0] = 0.5*n*L*c[0]/((double) N);  
*/
  for(j=0; j<N; j++) u[j]=j*L/N;
  four_coeffs_numrecFFT(cint, u, n); /* get coeffs of the integral of 1 */
  for(j=0; j<N; j++) cint[j] *= c[0]/N; 

  /* add terms coming from integrating everything but the c[0] term */
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] += -c[2*j] / (PI2_con*j);
    cint[2*j]   +=  c[2*j-1] / (PI2_con*j);
  }
  if( N%2 == 0 ) cint[N-1] += 0.;

  /* free temp mem u */  
  free(u);
}


/* compute Cheb coeffs c[0...n] from function u at the zeros of T_N(x).
   Note N=n+1                                                           */
void cheb_coeffs_fromZeros_numrecFFT(double c[], double u[], int n)
{
  int j;
  int N=n+1;
  double toN=2.0/N;

  for(j=0; j<N; j++)  c[j] = toN*u[j];
  numrec_cosft2(c-1, N, 1);
}

/* compute Cheb coeffs c[0...N] from function u at the extrema of T_N(x). */
void cheb_coeffs_fromExtrema_numrecFFT(double c[], double u[], int N)
{
  int j;
  double toN=2.0/N;

  for(j=0; j<=N; j++)  c[j] = toN*u[j];
  numrec_cosft1(c-1, N);
  c[N] *= 0.5;
}

/* find function u on the zeros of T_N(x).   Note N=n+1 */
void cheb_eval_onZeros_numrecFFT(double c[], double u[], int n)
{
  int j;
  int N=n+1;

  for(j=0; j<N; j++)  u[j] = c[j];
  numrec_cosft2(u-1, N, -1);
}

/* find function u on the extrema of T_N(X) */
void cheb_eval_onExtrema_numrecFFT(double c[], double u[], int N)
{
  int j;

  for(j=0; j<N; j++)  u[j] = c[j];
  u[N] = c[N]*2.0;
  numrec_cosft1(u-1, N);
}



/*********************************************************************/
/* FFTs using same conventions as in numrec                          */
/*********************************************************************/
#define SWAP(a,b) tmpRe=(a);(a)=(b);(b)=tmpRe
void numrec_four1(double data[], unsigned long nn, int isign)
{
  unsigned long n, mmax, m,i,j, istep;
  double wRe,wIm, wpRe,wpIm, theta;
  double wtmp, tmpRe, tmpIm;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2)
  {
    if (j > i)
    {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while(m >= 2 && j > m)
    {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  /*  Danielson-Lanczos algorithm: */
  mmax=2;
  while(n > mmax)
  {
    istep=mmax << 1;
    theta=isign*(twoPI/mmax);
    wtmp=sin(0.5*theta);
    wpRe = -2.*wtmp*wtmp;
    wpIm=sin(theta);
    wRe=1.;
    wIm=0.;
    for (m=1;m<mmax;m+=2)
    {
      for (i=m;i<=n;i+=istep)
      {
        j=i+mmax;
        tmpRe=wRe*data[j]-wIm*data[j+1];
        tmpIm=wRe*data[j+1]+wIm*data[j];
        data[j]=data[i]-tmpRe;
        data[j+1]=data[i+1]-tmpIm;
        data[i] += tmpRe;
        data[i+1] += tmpIm;
      }
      wtmp=wRe;
      wRe=wRe*wpRe-wIm*wpIm+wRe;
      wIm=wIm*wpRe+wtmp*wpIm+wIm;
    }
    mmax=istep;
  }
}
#undef SWAP

/* FFT of real valued func */
void numrec_realft(double data[], unsigned long n, int isign)
{
  unsigned long i, i1,i2,i3,i4, np3;
  double c1=0.5,c2;
  double h1r,h1i, h2r,h2i;
  double wRe,wIm, wpRe,wpIm;
  double wtmp, theta;

  theta=PI/(double) (n>>1);
  if(isign == 1)
  {
    c2 = -0.5;
    numrec_four1(data,n>>1,1);
  }
  else
  {
    c2=0.5;
    theta = -theta;
  }
  wtmp=sin(0.5*theta);
  wpRe = -2.*wtmp*wtmp;
  wpIm=sin(theta);
  wRe=1.+wpRe;
  wIm=wpIm;
  np3=n+3;
  for(i=2;i<=(n>>2);i++)
  {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wRe*h2r-wIm*h2i;
    data[i2]=h1i+wRe*h2i+wIm*h2r;
    data[i3]=h1r-wRe*h2r+wIm*h2i;
    data[i4] = -h1i+wRe*h2i+wIm*h2r;
    wtmp=wRe;
    wRe=wRe*wpRe-wIm*wpIm+wRe;
    wIm=wIm*wpRe+wtmp*wpIm+wIm;
  }
  if(isign == 1)
  {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  }
  else
  {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    numrec_four1(data,n>>1,-1);
  }
}

/* cos FT1 of numrec */
void numrec_cosft1(double y[], int n)
{
  int j, n2;
  double sum, y1,y2;
  double wIm=0.,wpIm, wpRe,wRe=1.;
  double wtmp, theta;

  theta=PI/n;
  wtmp=sin(0.5*theta);
  wpRe = -2.*wtmp*wtmp;
  wpIm=sin(theta);
  sum=0.5*(y[1]-y[n+1]);
  y[1]=0.5*(y[1]+y[n+1]);
  n2=n+2;
  for(j=2;j<=(n>>1);j++)
  {
    wtmp=wRe;
    wRe=wRe*wpRe-wIm*wpIm+wRe;
    wIm=wIm*wpRe+wtmp*wpIm+wIm;
    y1=0.5*(y[j]+y[n2-j]);
    y2=(y[j]-y[n2-j]);
    y[j]=y1-wIm*y2;
    y[n2-j]=y1+wIm*y2;
    sum += wRe*y2;
  }
  numrec_realft(y,n,1);
  y[n+1]=y[2];
  y[2]=sum;
  for(j=4;j<=n;j+=2)
  {
    sum += y[j];
    y[j]=sum;
  }
}

/* cos FT2 of numrec */
void numrec_cosft2(double y[], int n, int isign)
{
  int i;
  double sum,sum1,y1,y2,ytmp;
  double theta,wIm=0.,wIm1,wpIm,wpRe,wRe=1.,wRe1,wtmp;

  theta=0.5*PI/n;
  wRe1=cos(theta);
  wIm1=sin(theta);
  wpRe = -2.*wIm1*wIm1;
  wpIm=sin(2.*theta);
  if(isign == 1)
  {
    for(i=1;i<=n/2;i++)
    {
      y1=0.5*(y[i]+y[n-i+1]);
      y2=wIm1*(y[i]-y[n-i+1]);
      y[i]=y1+y2;
      y[n-i+1]=y1-y2;
      wtmp=wRe1;
      wRe1=wRe1*wpRe-wIm1*wpIm+wRe1;
      wIm1=wIm1*wpRe+wtmp*wpIm+wIm1;
    }
    numrec_realft(y,n,1);
    for(i=3;i<=n;i+=2)
    {
      wtmp=wRe;
      wRe=wRe*wpRe-wIm*wpIm+wRe;
      wIm=wIm*wpRe+wtmp*wpIm+wIm;
      y1=y[i]*wRe-y[i+1]*wIm;
      y2=y[i+1]*wRe+y[i]*wIm;
      y[i]=y1;
      y[i+1]=y2;
    }
    sum=0.5*y[2];
    for(i=n;i>=2;i-=2)
    {
      sum1=sum;
      sum += y[i];
      y[i]=sum1;
    }
  }
  else if(isign == -1)
  {
    ytmp=y[n];
    for(i=n;i>=4;i-=2) y[i]=y[i-2]-y[i];
    y[2]=2.*ytmp;
    for(i=3;i<=n;i+=2)
    {
      wtmp=wRe;
      wRe=wRe*wpRe-wIm*wpIm+wRe;
      wIm=wIm*wpRe+wtmp*wpIm+wIm;
      y1=y[i]*wRe+y[i+1]*wIm;
      y2=y[i+1]*wRe-y[i]*wIm;
      y[i]=y1;
      y[i+1]=y2;
    }
    numrec_realft(y,n,-1);
    for(i=1;i<=n/2;i++)
    {
      y1=y[i]+y[n-i+1];
      y2=(0.5/wIm1)*(y[i]-y[n-i+1]);
      y[i]=0.5*(y1+y2);
      y[n-i+1]=0.5*(y1-y2);
      wtmp=wRe1;
      wRe1=wRe1*wpRe-wIm1*wpIm+wRe1;
      wIm1=wIm1*wpRe+wtmp*wpIm+wIm1;
    }
  }
}


/*********************************************************************/
/* FFTW3 transforms                                                  */
/*********************************************************************/
#ifdef FFTW3
/* function to add special FFTW3 pars */
void add_special_FFTW3_pars(void)
{
  char str[1000];
  unsigned int flag;

  /* how we make FFTW3 plans */
  /* pars that contain values of FFTW3 planner flags */
  flag  = FFTW_MEASURE;         snprintf(str, 999, "%d", flag);
  AddPar("FFTW_MEASURE",         str, "value of FFTW_MEASURE");
  flag  = FFTW_DESTROY_INPUT;   snprintf(str, 999, "%d", flag);
  AddPar("FFTW_DESTROY_INPUT",   str, "value of FFTW_DESTROY_INPUT");
  flag  = FFTW_UNALIGNED;       snprintf(str, 999, "%d", flag);
  AddPar("FFTW_UNALIGNED",       str, "value of FFTW_UNALIGNED");
  flag  = FFTW_CONSERVE_MEMORY; snprintf(str, 999, "%d", flag);
  AddPar("FFTW_CONSERVE_MEMORY", str, "value of FFTW_CONSERVE_MEMORY");
  flag  = FFTW_EXHAUSTIVE;      snprintf(str, 999, "%d", flag);
  AddPar("FFTW_EXHAUSTIVE",      str, "value of FFTW_EXHAUSTIVE");
  flag  = FFTW_PRESERVE_INPUT;  snprintf(str, 999, "%d", flag);
  AddPar("FFTW_PRESERVE_INPUT",  str, "value of FFTW_PRESERVE_INPUT");
  flag  = FFTW_PATIENT;         snprintf(str, 999, "%d", flag);
  AddPar("FFTW_PATIENT",         str, "value of FFTW_PATIENT");
  flag  = FFTW_ESTIMATE;        snprintf(str, 999, "%d", flag);
  AddPar("FFTW_ESTIMATE",        str, "value of FFTW_ESTIMATE");
  //flag  = FFTW_WISDOM_ONLY;     snprintf(str, 999, "%d", flag);
  //AddPar("FFTW_WISDOM_ONLY",     str, "value of FFTW_WISDOM_ONLY");
}

/* initialize global plan arrays for the first time */
int init_FFTW3_plans(tGrid* grid)
{
  FFTW3_Nmax_of_plan_array = 0;
  FFTW3_dft_r2c_1d_plan = NULL;
  FFTW3_dft_c2r_1d_plan = NULL;
  FFTW_REDFT00_1d_plan = NULL;
  FFTW_REDFT10_1d_plan = NULL;
  FFTW_REDFT01_1d_plan = NULL;
  return reinit_FFTW3_plans(grid);
}
/* re-initialize global plan arrays */
int reinit_FFTW3_plans(tGrid* grid)
{
  int b;
  int N, Nmax=1;
  unsigned flags;
  char *flagstr;

  /* find max N on grid */
  for(b=0; b<Geti("nboxes"); b++)
  {
    int n1,n2,n3, m3;
    char str[1000];
    snprintf(str, 999, "box%d_n1", b);  n1=Geti(str);
    snprintf(str, 999, "box%d_n2", b);  n2=Geti(str);
    snprintf(str, 999, "box%d_n3", b);  n3=Geti(str);    
    m3=max3(n1,n2,n3)*4 + 1;
    if(m3>Nmax) Nmax=m3;
    // printf("m3=%d Nmax=%d", m3, Nmax);
  }

  /* allocate arrays of plans for all cases */
  if(FFTW3_Nmax_of_plan_array==0)
  {
    FFTW3_dft_r2c_1d_plan = (fftw_plan *) malloc(sizeof(fftw_plan) * (Nmax+1));
    FFTW3_dft_c2r_1d_plan = (fftw_plan *) malloc(sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT00_1d_plan = (fftw_plan *) malloc(sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT10_1d_plan = (fftw_plan *) malloc(sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT01_1d_plan = (fftw_plan *) malloc(sizeof(fftw_plan) * (Nmax+1));
  }
  else
  {
    FFTW3_dft_r2c_1d_plan = (fftw_plan *) realloc(FFTW3_dft_r2c_1d_plan, sizeof(fftw_plan) * (Nmax+1));
    FFTW3_dft_c2r_1d_plan = (fftw_plan *) realloc(FFTW3_dft_c2r_1d_plan, sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT00_1d_plan = (fftw_plan *) realloc(FFTW_REDFT00_1d_plan, sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT10_1d_plan = (fftw_plan *) realloc(FFTW_REDFT10_1d_plan, sizeof(fftw_plan) * (Nmax+1));
    FFTW_REDFT01_1d_plan = (fftw_plan *) realloc(FFTW_REDFT01_1d_plan, sizeof(fftw_plan) * (Nmax+1));
  }

  /* figure out FFTW planner flags we will use */
  printf("reinit_FFTW3_plans: planning for N=1,...,%d with flags:\n", Nmax);
  flags = 0;
  NextEntry(0);
  while( ( flagstr=NextEntry(Gets("FFTW3_planner_flags")) ) != NULL)
  {
    printf("%s ", flagstr);
    flags = flags | Geti(flagstr); 
  }
  printf("=> flags=%d\n", flags);
  fflush(stdout);

  /* make plans */
  for(N=1; N<=Nmax; N++)
  {
    int Nc = ( N / 2 ) + 1;
    fftw_complex *cco = fftw_malloc(sizeof(fftw_complex) * Nc);
    double *u = malloc(sizeof(double) * N);
    double *c = malloc(sizeof(double) * N);

    FFTW3_dft_r2c_1d_plan[N] = 
      fftw_plan_dft_r2c_1d(N, u, cco, flags);
    FFTW3_dft_c2r_1d_plan[N] =
      fftw_plan_dft_c2r_1d(N, cco, u, flags);
    FFTW_REDFT00_1d_plan[N] =
      fftw_plan_r2r_1d(N, u, c, FFTW_REDFT00, flags);
    FFTW_REDFT10_1d_plan[N] =
      fftw_plan_r2r_1d(N, u, c, FFTW_REDFT10, flags);
    FFTW_REDFT01_1d_plan[N] =
      fftw_plan_r2r_1d(N, u, c, FFTW_REDFT01, flags);

    free(c);
    free(u);
    fftw_free(cco);
  }
  FFTW3_Nmax_of_plan_array = Nmax;
  return 0;
} 
/* free the plans */
int free_FFTW3_plans(tGrid* grid)
{
  int N;
  /* deallocate the plans and all their associated data */
  for(N=1; N<=FFTW3_Nmax_of_plan_array; N++)
  {
    fftw_destroy_plan(FFTW3_dft_r2c_1d_plan[N]);
    fftw_destroy_plan(FFTW3_dft_c2r_1d_plan[N]);
    fftw_destroy_plan(FFTW_REDFT00_1d_plan[N]);
    fftw_destroy_plan(FFTW_REDFT10_1d_plan[N]);
    fftw_destroy_plan(FFTW_REDFT01_1d_plan[N]);
  }

  /* free arrays of plans */
  printf("free_FFTW3_plans: freeing global plan arrays.\n");
  free(FFTW3_dft_r2c_1d_plan);
  free(FFTW3_dft_c2r_1d_plan);
  free(FFTW_REDFT00_1d_plan);
  free(FFTW_REDFT10_1d_plan);
  free(FFTW_REDFT01_1d_plan);
  FFTW3_Nmax_of_plan_array = 0;
  return 0;
} 

/* compute Four coeffs c[0...n] from function u 
   at x_k = k/N, k=0,...,N-1 , N=n+1 
NOTE: four_coeffs returns c[] that are N times of those of four_coeffs_alt */
void four_coeffs_FFTW3(double *c, double *u, int n)
{
  int j, s;
  int N=n+1;
  int Nc = ( N / 2 ) + 1;
  fftw_complex *cco = fftw_malloc(sizeof(fftw_complex) * Nc);
  double *co = (double *) cco;

  /* execute right plan */
  fftw_execute_dft_r2c(FFTW3_dft_r2c_1d_plan[N], u, cco);

  /* convert */
  for(s=1, j=2; j<=N; j++, s=-s)  c[j-1] = s*co[j];
  c[0] = co[0];

  fftw_free(cco);
}

/* find function u from Four coeffs c[0...n], computed with four_coeffs */
void four_eval_FFTW3(double *c, double *u, int n)
{
  int j, s;
  int N=n+1;
  double ooN=1./N;
  int Nc = ( N / 2 ) + 1;
  fftw_complex *cco = fftw_malloc(sizeof(fftw_complex) * Nc);
  double *co = (double *) cco;

  /* convert */
  co[Nc*2-1] = 0.;
  for(s=1, j=2; j<=N; j++, s=-s)  co[j] = s*c[j-1]*ooN;
  co[1] = 0.;
  co[0] = c[0]*ooN;

  /* execute right plan */
  fftw_execute_dft_c2r(FFTW3_dft_c2r_1d_plan[N], cco, u);

  fftw_free(cco);
}

/* compute Cheb coeffs c[0...n] from function u at the zeros of T_N(x).
   Note N=n+1                                                           */
void cheb_coeffs_fromZeros_FFTW3(double *c, double *u, int n)
{
  int j;
  int N=n+1;
  double ooN=1./N;

  /* execute right plan */
  fftw_execute_r2r(FFTW_REDFT10_1d_plan[N], u, c);

  /* convert */
  for(j=0; j<N; j++)  c[j] *= ooN;
}

/* compute Cheb coeffs c[0...N] from function u at the extrema of T_N(x). */
void cheb_coeffs_fromExtrema_FFTW3(double *c, double *u, int N)
{
  int j;
  double ooN=1./N;

  /* execute right plan */
  fftw_execute_r2r(FFTW_REDFT00_1d_plan[N+1], u, c);

  /* convert */
  c[N] *= 0.5;
  for(j=0; j<=N; j++)  c[j] *= ooN;
}

/* find function u on the zeros of T_N(x).   Note N=n+1 */
void cheb_eval_onZeros_FFTW3(double *c, double *u, int n)
{
  int j;
  int N=n+1;

  /* execute right plan */
  fftw_execute_r2r(FFTW_REDFT01_1d_plan[N], c, u);

  /* convert */
  for(j=0; j<N; j++)  u[j] *= 0.5;
}

/* find function u on the extrema of T_N(X) */
void cheb_eval_onExtrema_FFTW3(double *c, double *u, int N)
{
  int j;
  double c_N = c[N];

  /* convert */
  c[N] *= 2.;

  /* execute right plan */
  fftw_execute_r2r(FFTW_REDFT00_1d_plan[N+1], c, u);

  /* convert */
  c[N] = c_N;
  for(j=0; j<=N; j++)  u[j] *= 0.5;
}

#else
void four_FFTW3_error(void)
{
  errorexit("four_coeffs_FFTW3/four_eval_FFTW3: in order to compile with fftw3\n"
            "use MyConfig with something like:\n"
            "DFLAGS += -DFFTW3\n"
            "FFTW3DIR = /opt/fftw-3.1.3\n"
            "SPECIALINCS += -I$(FFTW3DIR)/include\n"
            "SPECIALLIBS += -L$(FFTW3DIR)/lib -lfftw3\n"
            "SPECIALLIBS += -lfftw3");
}
void add_special_FFTW3_pars(void) {return;}
int init_FFTW3_plans(tGrid* grid) {return 0;}
int reinit_FFTW3_plans(tGrid* grid) {return 0;}
int free_FFTW3_plans(tGrid* grid) {return 0;}
void four_coeffs_FFTW3(double *c, double *u, int n) {four_FFTW3_error();}
void four_eval_FFTW3(double *c, double *u, int n) {four_FFTW3_error();}
void cheb_coeffs_fromZeros_FFTW3(double *c, double *u, int n) {four_FFTW3_error();}
void cheb_coeffs_fromExtrema_FFTW3(double *c, double *u, int N) {four_FFTW3_error();}
void cheb_eval_onZeros_FFTW3(double *c, double *u, int n) {four_FFTW3_error();}
void cheb_eval_onExtrema_FFTW3(double *c, double *u, int N) {four_FFTW3_error();}
#endif
