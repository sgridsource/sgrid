/* integrals.c */
/* Wolfgang Tichy 8/2005 */

#include "sgrid.h"
#include "Spectral.h"



/* compute U = integral of 3d var u over full interval 
   in direction direc on a box                         */
void spec_Integral1(tBox *box, int direc, double *u, double *U)
{
  int i,j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int N, double X)=NULL;
  int imethod=-42, chebmeth=1, fourmeth=2;

  get_spec_functionpointers(box, direc, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);

  if( get_coeffs==cheb_coeffs_fromExtrema ||
      get_coeffs==cheb_coeffs_fromZeros ||
      get_coeffs==cheb_coeffs_fromExtrema_numrecFFT ||
      get_coeffs==cheb_coeffs_fromZeros_numrecFFT ||
      get_coeffs==cheb_coeffs_fromExtrema_FFTW3 ||
      get_coeffs==cheb_coeffs_fromZeros_FFTW3          ) imethod = chebmeth;
  if( get_coeffs==four_coeffs || 
      get_coeffs==four_coeffs_numrecFFT ||
      get_coeffs==four_coeffs_FFTW3        ) imethod = fourmeth;

  if(direc==1)
  {
    /* write spec coeffs into U */
    spec_analysis1(box, direc, u, U);

    /* write cheb-integral from a to b into U */
    if(imethod==chebmeth)
      for(k = 0; k < n3; k++)
        for(j = 0; j < n2; j++)
        {
          double sum;
          double L = box->bbox[1] - box->bbox[0];

          sum = 0.5*U[Index(0,j,k)]*L;
          for(i = 2; i < n1; i+=2)
            sum += U[Index(i,j,k)]* L/(1.0-i*i);

          for(i = 0; i < n1; i++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(imethod==fourmeth)
      for(k = 0; k < n3; k++)
        for(j = 0; j < n2; j++)
        {
          double U0LoN = U[Index(0,j,k)]*(box->bbox[1]-box->bbox[0])/n1;

          for(i = 0; i < n1; i++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else if(direc==2)
  {
    /* write spec coeffs into U */
    spec_analysis1(box, direc, u, U);

    /* write cheb-integral from a to b into U */
    if(imethod==chebmeth)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[3] - box->bbox[2];

          sum = 0.5*U[Index(i,0,k)]*L;
          for(j = 2; j < n2; j+=2)
            sum += U[Index(i,j,k)]* L/(1.0-j*j);

          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(imethod==fourmeth)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double U0LoN = U[Index(i,0,k)]*(box->bbox[3]-box->bbox[2])/n2;

          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else if(direc==3)
  {
    /* write spec coeffs into U */
    spec_analysis1(box, direc, u, U);

    /* write cheb-integral from a to b into U */
    if(imethod==chebmeth)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[5] - box->bbox[4];

          sum = 0.5*U[Index(i,j,0)]*L;
          for(k = 2; k < n3; k+=2)
            sum += U[Index(i,j,k)]* L/(1.0-k*k);

          for(k = 0; k < n3; k++)
            U[Index(i,j,k)] = sum;
        }
    /* write four-integral from a to b into U */
    else if(imethod==fourmeth)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++)
        {
          double U0LoN = U[Index(i,j,0)]*(box->bbox[5]-box->bbox[4])/n3;

          for(k = 0; k < n3; k++)
            U[Index(i,j,k)] = U0LoN;
        }
    else errorexit("spec_Integral1: don't know how to integrate");
  }
  else
    errorexit("spec_Integral1: possible values for direction direc are 1,2,3.");
}


/* compute surface integral over surfaces normal to norm */
void spec_2dIntegral(tBox *box, int norm, double *u, double *U)
{
  if(norm==1)
  {
    spec_Integral1(box, 2, u, U);
    spec_Integral1(box, 3, U, U);
  }
  else if(norm==2)
  {
    spec_Integral1(box, 1, u, U);
    spec_Integral1(box, 3, U, U);
  }
  else if(norm==3)
  {
    spec_Integral1(box, 1, u, U);
    spec_Integral1(box, 2, U, U);
  }
  else
    errorexit("spec_2dIntegral: possible values for normal norm are 1,2,3.");
}


/* compute volume integral over 3d var u */
double spec_3dIntegral(tBox *box, double *u, double *U)
{
  spec_Integral1(box, 1, u, U);
  spec_Integral1(box, 2, U, U);
  spec_Integral1(box, 3, U, U);
  
  return U[0];
}


/* compute volume integral \int dx dy dz v(x,y,z)  of var v with
   index vind in a box. Here volume Jacobian is included. */
double VolumeIntegral_inBox(tBox *box, int vind)
{
  int idXdx = Ind("dXdx");
  int idYdx = Ind("dYdx");
  int idZdx = Ind("dZdx");
  double *dXdx  = box->v[idXdx];
  double *dXdy  = box->v[idXdx+1];
  double *dXdz  = box->v[idXdx+2];
  double *dYdx  = box->v[idYdx];
  double *dYdy  = box->v[idYdx+1];
  double *dYdz  = box->v[idYdx+2];
  double *dZdx  = box->v[idZdx];
  double *dZdy  = box->v[idZdx+1];
  double *dZdz  = box->v[idZdx+2];
  double *var   = box->v[vind];
  double *Integ = dmalloc(box->nnodes);
  double VolInt;
  int i;

  if( box->x_of_X[1] != NULL ) /* not Cartesian coords */
  {
    forallpoints(box, i)
    {
      double jac, det;

      if(dXdx!=NULL)
        det = dXdx[i]*dYdy[i]*dZdz[i] + dXdy[i]*dYdz[i]*dZdx[i] +
              dXdz[i]*dYdx[i]*dZdy[i] - dXdz[i]*dYdy[i]*dZdx[i] -
              dXdy[i]*dYdx[i]*dZdz[i] - dXdx[i]*dYdz[i]*dZdy[i];
      else
        errorexit("VolumeIntegral_inBox: implement dXdx==NULL case");

      if(det!=0.0) jac = 1.0/fabs(det);
      /* if det=0 jac should really be infinite, but we hope that the integrand
         goes to zero quickly enough that jac=0 makes no difference! */
      else jac = 0.0;

      Integ[i] = var[i] * jac;
    }
    /* integrate with jac */
    VolInt = spec_3dIntegral(box, Integ, Integ);
  }
  else /* integrate without jac */
    VolInt = spec_3dIntegral(box, var, Integ);
//printf("box%d %s VolInt=%g\n", box->b, VarName(vind), VolInt);

  free(Integ);
  return VolInt;
}


/* compute U = 2d integral of var u over theta and phi */
/* Note: U(r) = \int_0^{pi) dtheta \int_0^{2pi) dphi  
                u(r,theta,phi) |sin(theta)| r^2        */
void spec_sphericalDF2dIntegral(tBox *box, double *u, double *U)
{
  int i,j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int N, double X)=NULL;
  double *pX = box->v[Ind("X")];
  int imethod=-42, chebmeth=1, fourmeth=2;

  /* do phi integral */
  spec_Integral1(box, 3, u, U);

  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  if( get_coeffs==cheb_coeffs_fromExtrema ||
      get_coeffs==cheb_coeffs_fromZeros ||
      get_coeffs==cheb_coeffs_fromExtrema_numrecFFT ||
      get_coeffs==cheb_coeffs_fromZeros_numrecFFT   ||
      get_coeffs==cheb_coeffs_fromExtrema_FFTW3 ||
      get_coeffs==cheb_coeffs_fromZeros_FFTW3          ) imethod = chebmeth;
  if( get_coeffs==four_coeffs ||
      get_coeffs==four_coeffs_numrecFFT ||
      get_coeffs==four_coeffs_FFTW3        ) imethod = fourmeth;

  {
    /* write spec coeffs into U */
    spec_analysis1(box, 2, U, U);

    /* write four-integral from a to b into U */
    if(imethod==fourmeth)
      for (k = 0; k < n3; k++)
        for (i = 0; i < n1; i++)
        {
          double sum;
          double L = box->bbox[3] - box->bbox[2];
          int n;
          int N = box->n2;
          /* double theta = thm + PI/((1+N%2)*N); */
          double d = 1.0/(2.0*(1+N%2)*N);
          double PI2 = 2.0*PI;
          double Re_c_n, Im_c_n;

          /* sum over all theta-integrated terms */
          sum = (1.0/PI) * U[Index(i,0,k)];
          sum += 0.5*( sin(PI2*d) * U[Index(i,1,k)]
                      +cos(PI2*d) * U[Index(i,2,k)] );
          for(n=2;n<N/2;n+=2)
          {
            Re_c_n = U[Index(i, 2*n-1, k)]; /* c[2*n-1]; */
            Im_c_n = U[Index(i, 2*n, k)];   /* c[2*n];   */
            sum += 2.0*( cos(PI2*n*d)/((1-n*n)*PI) * Re_c_n 
                        +sin(PI2*n*d)/((n*n-1)*PI) * Im_c_n );
          }
          if( N%4 == 0 )
            sum += cos(PI2*(N/2)*d)/((1-N*N/4)*PI) * U[Index(i,N-1,k)];

          /* adjust sum for L and N to obtain integral over theta */
          sum *= L/N;

          /* write integral into U along direction 2, and multiply by r^2  */
          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = sum*pX[Index(i,j,k)]*pX[Index(i,j,k)];
        }
    else errorexit("spec_sphericalDF2dIntegral: you need Fourier in direction 2");
  }
}


/* compute volume integral over var u if we use sphericalDF */
double spec_sphericalDF3dIntegral(tBox *box, double *u, double *U)
{
  spec_sphericalDF2dIntegral(box, u, U);
  spec_Integral1(box, 1, U, U);
  
  return U[0];
}


/* same as spec_sphericalDF2dIntegral, but do it only at radial index i */
/* Note: U(i) = \int_0^{pi) dtheta \int_0^{2pi) dphi  
                u(r,theta,phi) |sin(theta)| r(i)^2        */
void spec_sphericalDF2dIntegral_at_radial_index_i(tBox *box, 
                                                  double *u, double *U, int i)
{
  int j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(void *, double, double, double *,double *, int)=NULL;
  void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(void *aux, double a, double b, int k, int N, double X)=NULL;
  double *pX = box->v[Ind("X")];
  int imethod=-42, chebmeth=1, fourmeth=2;

  /* do phi integral */
  spec_Integral1(box, 3, u, U);

  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &coeffs_of_int, &eval_onPoints,
                            &filter_coeffs, &basisfunc);
  if( get_coeffs==cheb_coeffs_fromExtrema ||
      get_coeffs==cheb_coeffs_fromZeros ||
      get_coeffs==cheb_coeffs_fromExtrema_numrecFFT ||
      get_coeffs==cheb_coeffs_fromZeros_numrecFFT   ||
      get_coeffs==cheb_coeffs_fromExtrema_FFTW3 ||
      get_coeffs==cheb_coeffs_fromZeros_FFTW3          ) imethod = chebmeth;
  if( get_coeffs==four_coeffs ||
      get_coeffs==four_coeffs_numrecFFT ||
      get_coeffs==four_coeffs_FFTW3        ) imethod = fourmeth;

  {
    /* write spec coeffs into U */
    spec_analysis1(box, 2, U, U);

    /* write four-integral from a to b into U */
    if(imethod==fourmeth)
      for (k = 0; k < n3; k++)
        {
          double sum;
          double L = box->bbox[3] - box->bbox[2];
          int n;
          int N = box->n2;
          /* double theta = thm + PI/((1+N%2)*N); */
          double d = 1.0/(2.0*(1+N%2)*N);
          double PI2 = 2.0*PI;
          double Re_c_n, Im_c_n;

          /* sum over all theta-integrated terms */
          sum = (1.0/PI) * U[Index(i,0,k)];
          sum += 0.5*( sin(PI2*d) * U[Index(i,1,k)]
                      +cos(PI2*d) * U[Index(i,2,k)] );
          for(n=2;n<N/2;n+=2)
          {
            Re_c_n = U[Index(i, 2*n-1, k)]; /* c[2*n-1]; */
            Im_c_n = U[Index(i, 2*n, k)];   /* c[2*n];   */
            sum += 2.0*( cos(PI2*n*d)/((1-n*n)*PI) * Re_c_n 
                        +sin(PI2*n*d)/((n*n-1)*PI) * Im_c_n );
          }
          if( N%4 == 0 )
            sum += cos(PI2*(N/2)*d)/((1-N*N/4)*PI) * U[Index(i,N-1,k)];

          /* adjust sum for L and N to obtain integral over theta */
          sum *= L/N;

          /* write integral into U along direction 2, and multiply by r^2  */
          for(j = 0; j < n2; j++)
            U[Index(i,j,k)] = sum*pX[Index(i,j,k)]*pX[Index(i,j,k)];
        }
    else errorexit("spec_sphericalDF2dIntegral: you need Fourier in direction 2");
  }
}


/* compute the indefinite integral U of 3d var u in direction direc on a box */
void spec_Int1(tBox *box, int direc, double *u, double *U)
{
  double *uline;
  double *Uline;
  int i,j,k, m3;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  /* memory for lines */
  if(direc==1)      m3=n1;
  else if(direc==2) m3=n2;
  else              m3=n3;
  uline = (double*) malloc(m3 * sizeof(double));
  Uline = (double*) malloc(m3 * sizeof(double));

  if(direc==1)
  {
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
      {
        /* 
        get_memline(u, uline, 1, j, k, n1, n2, n3);
        matrix_times_vector(box->Int1, uline, Uline, n1);
        put_memline(U, Uline, 1, j, k, n1, n2, n3);        
        */
        matrix_times_vector(box->Int1, u+Index(0,j,k), U+Index(0,j,k), n1);
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < n3; k++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 2, i, k, n1, n2, n3);
        matrix_times_vector(box->Int2, uline, Uline, n2);
        put_memline(U, Uline, 2, i, k, n1, n2, n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
      {
        get_memline(u, uline, 3, i, j, n1, n2, n3);
        matrix_times_vector(box->Int3 , uline, Uline, n3);
        put_memline(U, Uline, 3, i, j, n1, n2, n3);
      }
  }
  else
    errorexit("spec_Int1: possible values for direction direc are 1,2,3.");

  /* free lines */
  free(Uline);
  free(uline);
}


// ???: NOT TESTED YET !!!
/* compute the definite integral from a to b of the 3d var u 
   in direction direc on a box, save result in U */
void spec_definiteInt1(tBox *box, int direc, double a, double b, 
                       double *u, double *U)
{
  int linelen=0;
  int i,j,k;
  double Ua,Ub;
  double *Uline;
  double *BM[2];

  /* memory for lines and BM */
  linelen=max3(box->n1, box->n2, box->n3);
  Uline = (double *) malloc(linelen * sizeof(double));
  BM[0] = (double *) malloc(linelen * sizeof(double));
  BM[1] = (double *) malloc(linelen * sizeof(double));

  /* obtain BM vectors for interpolation along direc */
  spec_Basis_times_CoeffMatrix_direc(box, direc, BM[0], a);
  spec_Basis_times_CoeffMatrix_direc(box, direc, BM[1], b);

  /* get indefinite integral in U */
  spec_Int1(box, direc, u, U);

  /* now write result for def. int. into the entire Uline */
  if(direc==1)
  {
    for (k = 0; k < box->n3; k++)
      for (j = 0; j < box->n2; j++)
      {
        /* 
        get_memline(U, Uline, 1, j, k, box->n1, box->n2, box->n3);
        Ua = scalarproduct_vectors(BM[0], Uline, box->n1);
        Ub = scalarproduct_vectors(BM[1], Uline, box->n1);
        for (i = 0; i < n1; i++) Uline[i] = Ub-Ua;
        put_memline(U, Uline, 1, j, k, box->n1, box->n2, box->n3);        
        */
        int n1=box->n1;
        int n2=box->n2;
        Ua = scalarproduct_vectors(BM[0], U+Index(0,j,k), n1);
        Ub = scalarproduct_vectors(BM[1], U+Index(0,j,k), n1);
        for (i = 0; i < n1; i++) U[Index(i,j,k)] = Ub-Ua;
      }
  }
  else if(direc==2)
  {
    for (k = 0; k < box->n3; k++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(U, Uline, 2, i, k, box->n1, box->n2, box->n3);
        Ua = scalarproduct_vectors(BM[0], Uline, box->n2);
        Ub = scalarproduct_vectors(BM[1], Uline, box->n2);
        for (j = 0; j < box->n2; j++) Uline[j] = Ub-Ua;
        put_memline(U, Uline, 2, i, k, box->n1, box->n2, box->n3);        
      }
  }
  else if(direc==3)
  {
    for (j = 0; j < box->n2; j++)
      for (i = 0; i < box->n1; i++)
      {
        get_memline(U, Uline, 3, i, j, box->n1, box->n2, box->n3);
        Ua = scalarproduct_vectors(BM[0], Uline, box->n3);
        Ub = scalarproduct_vectors(BM[1], Uline, box->n3);
        for (k = 0; k < box->n3; k++) Uline[k] = Ub-Ua;
        put_memline(U, Uline, 3, i, j, box->n1, box->n2, box->n3);
      }
  }
  else
    errorexit("spec_definiteInt1: possible values for direction direc are 1,2,3.");

  free(Uline);  free(BM[0]);  free(BM[1]);
}
