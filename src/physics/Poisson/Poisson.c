/* Poisson.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "Poisson.h"

#define Power pow



void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);
void ABphi_of_xyz(tBox *box, double *A, double *B, double *phi,
                  double x, double y, double z);


/* initialize Poisson */
int Poisson_startup(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
*/
  printf("Initializing Poisson:\n");

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("Poisson_Psi", 0, 1, 1.0);
  VarNameSetBoundaryInfo("Poisson_Chi", 0, 1, 1.0);

  /* enable all vars */
  enablevar(grid, Ind("Poisson_Psi"));
  enablevar(grid, Ind("Poisson_rh1"));
  enablevar(grid, Ind("Poisson_Psix"));
  enablevar(grid, Ind("Poisson_Psixx"));
  enablevar(grid, Ind("Poisson_Err_Psi"));
  enablevar(grid, Ind("Poisson_Chi"));
  enablevar(grid, Ind("Poisson_rh2"));
  enablevar(grid, Ind("Poisson_Chix"));
  enablevar(grid, Ind("Poisson_Chixx"));
  enablevar(grid, Ind("Poisson_Err_Chi"));
  enablevar(grid, Ind("Poisson_temp1"));
  enablevar(grid, Ind("Poisson_temp2"));
  enablevar(grid, Ind("Poisson_temp3"));
  enablevar(grid, Ind("Poisson_temp4"));
  enablevar(grid, Ind("Poisson_temp5"));
  enablevar(grid, Ind("Poisson_temp6"));
  
  /* set initial values in all in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *Psi = box->v[Ind("Poisson_Psi")];
    double *Chi = box->v[Ind("Poisson_Chi")];
    double *rh1 = box->v[Ind("Poisson_rh1")];
    double *rh2 = box->v[Ind("Poisson_rh2")];

    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      /* set Psi and Chi */
      if(Getv("Poisson_grid", "SphericalDF"))
      {
        Psi[i] = x*y*z; // 1.0/sqrt(x*x + y*y + z*z);
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
      }
      else if(Getv("Poisson_grid", "AnsorgNS") || 
              Getv("Poisson_grid", "4ABphi_2xyz"))
      {
        Psi[i] = 0.0;
        Chi[i] = 0.0;
        rh1[i] = 0.0;
        rh2[i] = 0.0;
        if(b==0 || b==5)  rh1[i] = -4.0*PI;
        if(b==3 || b==4)  rh2[i] = -8.0*PI;
/*
        Psi[i] = b*b;
        Psi[i] = 1*exp(-0.5*(
                 (x-1)*(x-1)/(1*1) + (y-0)*(y-0)/(1.5*1.5)
               + (z-0.5)*(z-0.5)/(2*2)   ) );
        Psi[i] = exp(-0.5*x*x);
*/
        Psi[i] = 1.0/sqrt(x*x + y*y + z*z+1);
        Psi[i] = 1.0/sqrt(x*x + (y-1)*(y-1) + (z-0.2)*(z-0.2)+1);
        //Psi[i] = 0.0001*(b+1)/sqrt(x*x + y*y + z*z+1);
        //if(b==0 || b==3)
        //  Psi[i] = (b-1)*x;
      }
    }
  }

  return 0;
}


/* Solve the Equations */
int Poisson_solve(tGrid *grid)
{
  int    itmax = Geti("Poisson_itmax");
  double tol   = Getd("Poisson_tol");
  int    linSolver_itmax  = Geti("Poisson_linSolver_itmax");
  double linSolver_tolFac = Getd("Poisson_linSolver_tolFac");
  double normresnonlin;

  tVarList *vlu, *vlFu, *vluDerivs;
  tVarList *vldu, *vlr, *vlduDerivs;
  tVarList *vlrhs;

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vlFu = vlalloc(grid);
  vluDerivs= vlalloc(grid);
  vlrhs = vlalloc(grid);

  /* add Poisson_Psi and Poisson_Chi to vlu */
  vlpush(vlu, Ind("Poisson_Psi"));
  vlpush(vlu, Ind("Poisson_Chi"));

  /* add derivs of Poisson_Psi and Poisson_Chi to vluDerivs */
  vlpush(vluDerivs, Ind("Poisson_Psix"));
  vlpush(vluDerivs, Ind("Poisson_Psixx"));
  vlpush(vluDerivs, Ind("Poisson_Chix"));
  vlpush(vluDerivs, Ind("Poisson_Chixx"));

  /* add Poisson_rh1 and Poisson_rh2 to vlrhs */
  vlpush(vlrhs, Ind("Poisson_rh1"));
  vlpush(vlrhs, Ind("Poisson_rh2"));

  /* add Poisson_Err_Psi and Poisson_Err_Chi to vlFu */
  vlpush(vlFu, Ind("Poisson_Err_Psi"));
  vlpush(vlFu, Ind("Poisson_Err_Chi"));
    
  /* enable vlu, vluDerivs, vlrhs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 
  enablevarlist(vlrhs);
  enablevarlist(vlFu);

  /* now duplicate vlu and vluDerivs for linarized Eqs. */
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlr        = AddDuplicateEnable(vlFu, "_l");
  vlduDerivs = vluDerivs; /* maybe: vlduDerivs=AddDuplicateEnable(vluDerivs, "_l"); */

  /* call Newton solver */
/*
{
double Y;
double X=1;
for(Y=0;Y<=1;Y+=0.1)
{
printf("X=%g Y=%g: ",X,Y);
printf("x0=%g y0=%g  ",
grid->box[0]->x_of_X[1]((void *) grid->box[0], 0, X, Y, 0),
grid->box[0]->x_of_X[2]((void *) grid->box[0], 0, X, Y, 0));

printf("x1=%g y1=%g\n",
grid->box[1]->x_of_X[1]((void *) grid->box[1], 0, X, Y, 0),
grid->box[1]->x_of_X[2]((void *) grid->box[1], 0, X, Y, 0));
}
}
*/
/*
{
double X,Y,Z;
double x,y,z;
int i;
tBox *box=grid->box[0];
x=0.6; y=0.0; z=0.0;
printf("nearestXYZ=%f ", nearestXYZ_of_xyz(box, &i, &X,&Y,&Z, x,y,z));
printf("(X,Y,Z)=(%.12f,%.12f,%.12f)\n", X,Y,Z);
Y+=1e-10;
ABphi_of_xyz(box, &X,&Y,&Z, x,y,z);
printf("(x,y,z)=(%f,%f,%f)   (X,Y,Z)=(%.12f,%.12f,%.12f)\n", x,y,z, X,Y,Z);
}
*/
//printf("calling write_grid(grid)\n");
//write_grid(grid);
F_Poisson(vlFu, vlu, vluDerivs, vlrhs);
printf("calling write_grid(grid)\n");
write_grid(grid);
//exit(11);
vlu->n=1; vlFu->n=1; vluDerivs->n=1; vlrhs->n=1;
vldu->n=1; vlr->n=1; vlduDerivs->n=1; vlrhs->n=1;

/*
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, vlrhs,
         itmax, tol, &normresnonlin, 1,
         bicgstab, Precon_I, vldu, vlr, vlduDerivs, vlrhs,
         linSolver_itmax, linSolver_tolFac);
*/

  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, vlrhs,
         itmax, tol, &normresnonlin, 1,
         LinSolve_withLAPACK, Precon_I, vldu, vlr, vlduDerivs, vlrhs,
         linSolver_itmax, linSolver_tolFac);

/*
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, vlrhs,
         1, tol, &normresnonlin, 1,
         bicgstab, Precon_I, vldu, vlr, vlduDerivs, vlrhs,
         2, linSolver_tolFac);
*/

  /* free varlists */     
  VLDisableFree(vldu);
  /* VLDisableFree(vlduDerivs); */
  VLDisableFree(vlr);     
  vlfree(vlu);
  vlfree(vluDerivs);
  vlfree(vlrhs);
  vlfree(vlFu);
        
  return 0;
}


/* compute absolute error in ANALYSIS */
int Poisson_analyze(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
*/
  printf("Poisson: computing absolute error\n");
  
  /* set initial data in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *Psi    = box->v[Ind("Poisson_Psi")];
    double *PsiErr = box->v[Ind("Poisson_Err_Psi")];
    double *Chi    = box->v[Ind("Poisson_Chi")];
    double *ChiErr = box->v[Ind("Poisson_Err_Chi")];


    /* subtract true values */
    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      if(Getv("Poisson_grid", "SphericalDF"))
      {
        //PsiErr[i] = Psi[i]-1.0/sqrt(x*x + y*y + z*z);
        //ChiErr[i] = Chi[i]-2.0/sqrt(x*x + y*y + z*z);
      }
      else if(Getv("Poisson_grid", "AnsorgNS") || 
              Getv("Poisson_grid", "4ABphi_2xyz"))
      {
/*
        if(b==1||b==2||b==3)
          PsiErr[i] = Psi[i]-1.0/sqrt((x-1)*(x-1) + y*y + z*z);
        if(b==0||b==1||b==2)
          ChiErr[i] = Chi[i]-2.0/sqrt((x+1)*(x+1) + y*y + z*z);
*/
      }
    }
  }
  return 0;
}


void filterSing(tBox *box, double *u)
{
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  int filt2=n2-1;
  int filt3=n3-1; // n3-1; // n3/2;
  int i,j,k;
  void (*get_coeffs)(double *,double *, int)=NULL;
  void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
  void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
  void (*eval_onPoints)(double *,double *, int)=NULL;
  void (*filter_coeffs)(double *, int, int)=NULL;
  double (*basisfunc)(double a, double b, int k, int N, double X)=NULL;
  static int linelen=0;
  static double *uline=NULL;
  static double *ufline=NULL;
  int m3;

  get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints, 
                            &filter_coeffs, &basisfunc);
  initfiltermatrix(box->F2, n2+1-filt2, n2, 
                   get_coeffs, filter_coeffs, eval_onPoints);

  get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                            &coeffs_of_2ndderiv, &eval_onPoints, 
                            &filter_coeffs, &basisfunc);
  initfiltermatrix(box->F3, n3+1-filt3, n3, 
                   get_coeffs, filter_coeffs, eval_onPoints);

  /* static memory for lines */
  m3=max3(box->n1, box->n2, box->n3);
  if(m3>linelen)
  {
    linelen = m3;
    uline = (double*) realloc(uline, linelen * sizeof(double));
    ufline = (double*) realloc(ufline, linelen * sizeof(double));
  }

  i=box->n1-1;
//  /* filter along phi-direc if A=1 for B=Bmin and B=Bmax */
//  for(j=0;j<n2;j=j+n2-1)
  /* filter along phi-direc if A=1 for all B */
  for(j=0;j<n2;j++)
  {
    get_memline(u, uline,  3, i, j, n1, n2, n3);
    matrix_times_vector(box->F3 , uline, ufline, box->n3);
    put_memline(u, ufline, 3, i, j, n1, n2, n3);
  }
  /* filter along B-direc if A=1 for all phi */
  for(k=0;k<n3;k++)
  {
    get_memline(u, uline,  2, i, k, n1, n2, n3);
    matrix_times_vector(box->F2, uline, ufline, box->n2);
    put_memline(u, ufline, 2, i, k, n1, n2, n3);
  }
  /* filter along phi-direc if A=1 for B=Bmin and B=Bmax */
  //  for(j=0;j<n2;j=j+n2-1)

}                




/* evaluate Poisson eqn for vlu */
void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlrhs)
{
  tGrid *grid = vlu->grid;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[vlFu->index[0]];
    // double *Psi  = box->v[vlu->index[0]];
    double *rh1  = box->v[vlrhs->index[0]];
    double *Psixx = box->v[vluDerivs->index[3]];
    double *Psiyy = box->v[vluDerivs->index[6]];
    double *Psizz = box->v[vluDerivs->index[8]];
    double *FChi = box->v[vlFu->index[1]];
    // double *Chi  = box->v[vlu->index[1]];
    double *rh2  = box->v[vlrhs->index[1]];
    double *Chixx = box->v[vluDerivs->index[12]];
    double *Chiyy = box->v[vluDerivs->index[15]];
    double *Chizz = box->v[vluDerivs->index[17]];
    int i;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes"))
    {
      allDerivsOf_S(box, vlu->index[0], vluDerivs->index[0],
                    vluDerivs->index[3]);
      allDerivsOf_S(box, vlu->index[1], vluDerivs->index[9],
                    vluDerivs->index[12]);
    }
    else
    {
      FirstAndSecondDerivsOf_S(box, vlu->index[0], vluDerivs->index[0],
                               vluDerivs->index[3]);
      FirstAndSecondDerivsOf_S(box, vlu->index[1], vluDerivs->index[9],
                               vluDerivs->index[12]);
    }
    
    if(b==0 & 0)
    {
      double *Psix = box->v[vluDerivs->index[0]];
      double *Psiy = box->v[vluDerivs->index[1]];
      double *Psiz = box->v[vluDerivs->index[2]];
      double *Psixy = box->v[vluDerivs->index[4]];
      double *Psixz = box->v[vluDerivs->index[5]];
      double *Psiyz = box->v[vluDerivs->index[7]];
      // //FirstDerivsOf_S(box, vlu->index[0], vluDerivs->index[0]);
      filterSing(box, Psix); // ???
      filterSing(box, Psiy);
      filterSing(box, Psiz);
      FirstDerivsOf_S(box, vluDerivs->index[2],vluDerivs->index[6]);
      FirstDerivsOf_S(box, vluDerivs->index[1],vluDerivs->index[5]);
      FirstDerivsOf_S(box, vluDerivs->index[0],vluDerivs->index[3]);
      //filterSing(box, Psixx);
      //filterSing(box, Psiyy);
      //filterSing(box, Psizz);
      ////filterSing(box, Psiyz);
    }

//    forallpoints(box, i)  FPsi[i] = Psi[i]*Psi[i] - 9.0;
    forallpoints(box, i)
    {
      FPsi[i] = Psixx[i] + Psiyy[i] + Psizz[i] - rh1[i];
      FChi[i] = Chixx[i] + Chiyy[i] + Chizz[i] - rh2[i];
    }

/*
double *Psi  = box->v[vlu->index[0]];
double *Psix = box->v[vluDerivs->index[0]];
double *Psiy = box->v[vluDerivs->index[1]];
double *Psiz = box->v[vluDerivs->index[2]];
double *Psixy = box->v[vluDerivs->index[4]];
double *Psixz = box->v[vluDerivs->index[5]];
double *Psiyz = box->v[vluDerivs->index[7]];

double *M = box->v[vluDerivs->index[4]];
double *c = box->v[vluDerivs->index[5]];
double *d = box->v[vluDerivs->index[7]];

spec_Coeffs(box, Psi, c);

double X,Y,Z;
double x,y,z;


if(b<4)
{
x=1.2; y=-0.5; z=0.5;
printf("nearestXYZ=%f ", nearestXYZ_of_xyz(box, &i, &X,&Y,&Z, x,y,z));
printf("(X,Y,Z)=(%.12f,%.12f,%.12f)\n", X,Y,Z);
//X=0.1; Y=0.3; 
Z=PI-PI/4;
XYZ_of_xyz(box, &X,&Y,&Z, x,y,z);
printf("(x,y,z)=(%f,%f,%f)   (X,Y,Z)=(%.12f,%.12f,%.12f)\n", x,y,z, X,Y,Z);
printf("##### Psi=%.12f\n", spec_interpolate(box, c, X,Y,Z));
}
*/

    /* BCs */
    set_BCs(vlFu, vlu, vluDerivs, 1);
  }
}


/* evaluate linearized Laplace */
void J_Poisson(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  tGrid *grid = vldu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *JdPsi = box->v[vlJdu->index[0]];
    // double *dPsi  = box->v[vldu->index[0]];
    double *dPsixx = box->v[vlduDerivs->index[3]];
    double *dPsiyy = box->v[vlduDerivs->index[6]];
    double *dPsizz = box->v[vlduDerivs->index[8]];
    double *JdChi = box->v[vlJdu->index[1]];
    // double *dChi  = box->v[vldu->index[1]];
    double *dChixx = box->v[vlduDerivs->index[12]];
    double *dChiyy = box->v[vlduDerivs->index[15]];
    double *dChizz = box->v[vlduDerivs->index[17]];
    int i;

    /* compute the derivs */
    if(Getv("Poisson_useDD", "yes"))
    {
      allDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                    vlduDerivs->index[3]);
      allDerivsOf_S(box, vldu->index[1], vlduDerivs->index[9],
                    vlduDerivs->index[12]);
    }
    else
    {
      FirstAndSecondDerivsOf_S(box, vldu->index[0], vlduDerivs->index[0],
                               vlduDerivs->index[3]);
      FirstAndSecondDerivsOf_S(box, vldu->index[1], vlduDerivs->index[9],
                               vlduDerivs->index[12]);
    }
//    forallpoints(box, i)  JdPsi[i] = 2.0*Psi[i]*dPsi[i];
    forallpoints(box, i)
    {
      JdPsi[i] = dPsixx[i] + dPsiyy[i] + dPsizz[i];
      JdChi[i] = dChixx[i] + dChiyy[i] + dChizz[i];
    }

    /* BCs */
    set_BCs(vlJdu, vldu, vlduDerivs, 0);
  }
}


/* Einheitsmatrix als Precon */ 
void Precon_I(tVarList *vlJdu, tVarList *vldu,
              tVarList *vlduDerivs, tVarList *vlu)
{
  tGrid *grid = vldu->grid;
  int i,j,b;
  	
  for(j = 0; j < vldu->n; j++)
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *Jdu = box->v[vlJdu->index[j]];
      double *du  = box->v[vldu->index[j]];

      forallpoints(box, i)  Jdu[i] = du[i];
    }
}


void set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  	
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *FPsi = box->v[vlFu->index[0]];
    double *Psi  = box->v[vlu->index[0]];
    double *Psix = box->v[vluDerivs->index[0]];
    double *Psiy = box->v[vluDerivs->index[1]];
    double *Psiz = box->v[vluDerivs->index[2]];
    double *FChi = box->v[vlFu->index[1]];
    double *Chi  = box->v[vlu->index[1]];
    double *Chix = box->v[vluDerivs->index[9]];
    double *Chiy = box->v[vluDerivs->index[10]];
    double *Chiz = box->v[vluDerivs->index[11]];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int i,j,k;

    /* BCs */
    if(Getv("Poisson_grid", "SphericalDF"))
    {
      forplane1(i,j,k, n1,n2,n3, 0)
      {
        FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*nonlin;
        FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 2.0*nonlin;
      }

      forplane1(i,j,k, n1,n2,n3, n1-1)
      {
        FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*nonlin;
        FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - 1.0*nonlin;
      }
    }
    else if (Getv("Poisson_grid", "AnsorgNS"))
    {
      double *P;
      double *dP[4];
      double *C;
      double *dC[4];
      double *BM;

      if(b==0)  /* in box0 */
      {
        BM = (double *) calloc(n1, sizeof(double));

        /* values at A=0 are equal in box0 and box1 */
        P = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C = grid->box[1]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0.0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double *line = (double *) calloc(n1, sizeof(double));
          double U0;

          /* find values U0 in box0 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];

          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
          free(line);
        }
        free(BM);
      }
      else if(b==3)  /* in box3 */
      {
        BM = (double *) calloc(n1, sizeof(double));

        /* values at A=0 are equal in box3 and box2 */
        P = grid->box[2]->v[vlu->index[0]]; /* values in box2 */
        C = grid->box[2]->v[vlu->index[1]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, n1, BM, 0.0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          int l;
          double *line = (double *) calloc(n1, sizeof(double));
          double U0;

          /* find values U0 in box3 at A=0*/
          get_memline(Psi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
//printf("U0=%f P[Index(i,j,k)]=%f\n", U0, P[Index(i,j,k)]);
          get_memline(Chi, line, 1, j,k, n1,n2,n3);
          for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
          FChi[Index(i,j,k)] = U0 - C[Index(i,j,k)];
          free(line);
        }
        free(BM);
      }
      else if(b==1)
      {
        BM = (double *) calloc(max3(grid->box[0]->n1, box->n1, box->n2), 
                               sizeof(double));

        /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
        dP[1] = grid->box[2]->v[vluDerivs->index[0]];
        dC[1] = grid->box[2]->v[vluDerivs->index[9]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          FChi[Index(i,j,k)] = Chix[Index(i,j,k)] - dC[1][Index(i,j,k)];
        }

        /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[0]->v[vluDerivs->index[0]];
        dP[2] = grid->box[0]->v[vluDerivs->index[1]];
        dP[3] = grid->box[0]->v[vluDerivs->index[2]];
        dC[1] = grid->box[0]->v[vluDerivs->index[9]];
        dC[2] = grid->box[0]->v[vluDerivs->index[10]];
        dC[3] = grid->box[0]->v[vluDerivs->index[11]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, grid->box[0]->n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];
          double DP[4], DC[4];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box0 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[0]->n1;
            double *line = (double *) calloc(n1, sizeof(double));

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DC[m] += BM[l]*line[l];
            free(line);
          }
          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - DC[1]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - DC[2]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - DC[3]);
        }
        
        /* Psi=Chi=0 at infinity */
        // fixme: B=0 is not on grid!!!
        for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
        {
          FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
        }
        free(BM);
      }
      else if(b==2)
      {
        BM = (double *) calloc(max3(grid->box[3]->n1, box->n1, box->n2), 
                               sizeof(double));

        /* values at A=1 are equal in box1 and box2 */
        P  = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C  = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[3]->v[vluDerivs->index[0]];
        dP[2] = grid->box[3]->v[vluDerivs->index[1]];
        dP[3] = grid->box[3]->v[vluDerivs->index[2]];
        dC[1] = grid->box[3]->v[vluDerivs->index[9]];
        dC[2] = grid->box[3]->v[vluDerivs->index[10]];
        dC[3] = grid->box[3]->v[vluDerivs->index[11]];
        spec_Basis_times_CoeffMatrix(0.0,1.0, grid->box[3]->n1, BM, 0,
                                     cheb_coeffs_fromZeros, cheb_basisfunc);
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];
          double DP[4], DC[4];
          int m,l;
          /* find derivs of Psi, Chi at A=0 in box3 
             and store them in DP[m], DC[m] */
          for(m=1; m<=3; m++)
          {
            int n1 = grid->box[3]->n1;
            double *line = (double *) calloc(n1, sizeof(double));

            DP[m] = 0.0;
            get_memline(dP[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
            DC[m] = 0.0;
            get_memline(dC[m], line, 1, j,k, n1,n2,n3);
            for(l=0; l<n1; l++)  DC[m] += BM[l]*line[l];
            free(line);
          }
          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - DC[1]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - DC[2]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - DC[3]);
        }

        /* Psi=Chi=0 at infinity */
        // fixme: B=0 is not on grid!!!
        for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
        {
          FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
          FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
        }

        free(BM);
      }
      else errorexiti("b=%d should be impossible!", b);
    } /* end: else if (Getv("Poisson_grid", "AnsorgNS")) */
    else if (Getv("Poisson_grid", "4ABphi_2xyz"))
    {
      double *P;
      double *dP[4];
      double *C;
      double *dC[4];
      double *X, *Y, *Z;
      double *Pcoeffs, *Ccoeffs;
      double Pinterp, Cinterp;
      double x,y,z;

      /* special rho=0 case??? */
      if(b==0 || b==1 || b==2 || b==3)
      {
        int pl;
        char str[1000];
        snprintf(str, 999, "box%d_basis2", b);
        if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* loop over rho=0 boundary */
          for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
            forplane2(i,j,k, n1,n2,n3, pl)
            {
              if(k>0) /* phi>0: impose u_phi_phi=0 */
              {
                FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                FChi[Index(i,j,k)] = Chi_phi_phi[Index(i,j,k)];
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              {
                double Psi_rho = Psiy[Index(i,j,k)];
                double Chi_rho = Chiy[Index(i,j,k)];
                double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                double Chi_rho_phi_phi = Chi_y_phi_phi[Index(i,j,k)];
                FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                FChi[Index(i,j,k)] = Chi_rho + Chi_rho_phi_phi;
              }
            }
        }
        /* same as before, but also interpolate to rho=0 */
        else if(Getv("Poisson_4ABphi_2xyz", "regularity_on_axis"))
        {
          double *Psi_phi_phi = box->v[Ind("Poisson_temp1")];
          double *Chi_phi_phi = box->v[Ind("Poisson_temp2")];
          double *Psi_y_phi_phi = box->v[Ind("Poisson_temp3")];
          double *Chi_y_phi_phi = box->v[Ind("Poisson_temp4")];
          double *temp5 = box->v[Ind("Poisson_temp5")];
          double *temp6 = box->v[Ind("Poisson_temp6")];
          double *line = (double *) calloc(n2, sizeof(double));
          double *BM[2];
          BM[0] = (double *) calloc(n2, sizeof(double));
          BM[1] = (double *) calloc(n2, sizeof(double));

          /* get u_phi_phi */
          spec_Deriv2(box, 3, Psi, Psi_phi_phi);
          spec_Deriv2(box, 3, Chi, Chi_phi_phi);
          
          /* get u_rho_phi_phi at phi=0 */
          /* d/drho = dx^i/drho d/dx^i, 
             dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
             ==> d/drho u = d/dy u  at phi=0 */           
          /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
          cart_partials(box, Psi_phi_phi, temp5, Psi_y_phi_phi, temp6);
          cart_partials(box, Chi_phi_phi, temp5, Chi_y_phi_phi, temp6);

          /* obtain BM vectors for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[0], 0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM[1], 1,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);

          /* loop over rho~0 boundary */
          for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
          {
            int l;
            double U0, V0;

            forplane2(i,j,k, n1,n2,n3, pl)
            {
              if(k>0) /* phi>0: impose u_phi_phi=0 */
              {
                /* find value Psi_phi_phi at B=0 or 1 */
                get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  U0 += BM[j>0][l]*line[l];
                FPsi[Index(i,j,k)] = U0;

                /* find value Chi_phi_phi at B=0 or 1 */
                get_memline(Chi_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  U0 += BM[j>0][l]*line[l];
                FChi[Index(i,j,k)] = U0;
              }
              else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
              { /* Psi_rho = Psiy  
                   Chi_rho = Chiy
                   Psi_rho_phi_phi = Psi_y_phi_phi
                   Chi_rho_phi_phi = Chi_y_phi_phi */
                /* find value Psi_rho at B=0 or 1 */
                get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Psi_rho_phi_phi at B=0 or 1 */
                get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  V0 += BM[j>0][l]*line[l];

                FPsi[Index(i,j,k)] = U0 + V0;

                /* find value Chi_rho at B=0 or 1 */
                get_memline(Chiy, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  U0 += BM[j>0][l]*line[l];

                /* find value Chi_rho_phi_phi at B=0 or 1 */
                get_memline(Chi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n1; l++)  V0 += BM[j>0][l]*line[l];

                FChi[Index(i,j,k)] = U0 + V0;
              }
            }
          }
          free(BM[0]);
          free(BM[1]);
          free(line);
        }
      } /* end: special rho=0 case??? */

      if(b==0)  /* in box0 */
      {
        /* values at A=0 are equal in box0 and box1 */
        P = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* values at A=Amin are interpolated from box5 */
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[5]->v[vlu->index[0]]; /* values in box5 */
        C = grid->box[5]->v[vlu->index[1]];
        Pcoeffs = grid->box[5]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[5]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[5], P, Pcoeffs);
        spec_Coeffs(grid->box[5], C, Ccoeffs);
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
        {
          int ind=Index(i,j,k);
          x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          Pinterp = spec_interpolate(grid->box[5], Pcoeffs, x,y,z);
          Cinterp = spec_interpolate(grid->box[5], Ccoeffs, x,y,z);
          FPsi[ind] = Psi[ind] - Pinterp;
          FChi[ind] = Chi[ind] - Cinterp;
        }
      }
      else if(b==3)  /* in box3 */
      {
        /* values at A=0 are equal in box3 and box2 */
        P = grid->box[2]->v[vlu->index[0]]; /* values in box2 */
        C = grid->box[2]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }
        /* values at A=Amin are interpolated from box4 */
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[4]->v[vlu->index[0]]; /* values in box4 */
        C = grid->box[4]->v[vlu->index[1]];
        Pcoeffs = grid->box[4]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[4]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[4], P, Pcoeffs);
        spec_Coeffs(grid->box[4], C, Ccoeffs);
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
        {
          int ind=Index(i,j,k);
          x = box->x_of_X[1]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          y = box->x_of_X[2]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          z = box->x_of_X[3]((void *) box, ind, X[ind],Y[ind],Z[ind]); 
          Pinterp = spec_interpolate(grid->box[4], Pcoeffs, x,y,z);
          Cinterp = spec_interpolate(grid->box[4], Ccoeffs, x,y,z);
          FPsi[ind] = Psi[ind] - Pinterp;
          FChi[ind] = Chi[ind] - Cinterp;
        }
      }
      else if(b==1)  /* in box1 */
      {
        /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
        dP[1] = grid->box[2]->v[vluDerivs->index[0]];
        dC[1] = grid->box[2]->v[vluDerivs->index[9]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          FChi[Index(i,j,k)] = Chix[Index(i,j,k)] - dC[1][Index(i,j,k)];
        }

        /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[0]->v[vluDerivs->index[0]];
        dP[2] = grid->box[0]->v[vluDerivs->index[1]];
        dP[3] = grid->box[0]->v[vluDerivs->index[2]];
        dC[1] = grid->box[0]->v[vluDerivs->index[9]];
        dC[2] = grid->box[0]->v[vluDerivs->index[10]];
        dC[3] = grid->box[0]->v[vluDerivs->index[11]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];

          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - dC[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - dC[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - dC[3][Index(i,j,k)]);
        }
        
        /* Psi=Chi=0 at infinity */
        if(Getv("box1_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *BM = (double *) calloc(n2, sizeof(double));
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
          }
          free(line);
          free(BM);
        }
      }
      else if(b==2)  /* in box2 */
      {
        /* values at A=1 are equal in box1 and box2 */
        P  = grid->box[1]->v[vlu->index[0]]; /* values in box1 */
        C  = grid->box[1]->v[vlu->index[1]];
        forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
        {
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          FChi[Index(i,j,k)] = Chi[Index(i,j,k)] - C[Index(i,j,k)];
        }

        /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
        /* Below we use the approximate normal vec 
           ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
        dP[1] = grid->box[3]->v[vluDerivs->index[0]];
        dP[2] = grid->box[3]->v[vluDerivs->index[1]];
        dP[3] = grid->box[3]->v[vluDerivs->index[2]];
        dC[1] = grid->box[3]->v[vluDerivs->index[9]];
        dC[2] = grid->box[3]->v[vluDerivs->index[10]];
        dC[3] = grid->box[3]->v[vluDerivs->index[11]];
        forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
        {
          double B   = box->v[Ind("Y")][Index(i,j,k)];
          double phi = box->v[Ind("Z")][Index(i,j,k)];

          FPsi[Index(i,j,k)] = 
           cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          FChi[Index(i,j,k)] = 
           cos(PI*B)         * (Chix[Index(i,j,k)] - dC[1][Index(i,j,k)]) +
           sin(PI*B)*cos(phi)* (Chiy[Index(i,j,k)] - dC[2][Index(i,j,k)]) +
           sin(PI*B)*sin(phi)* (Chiz[Index(i,j,k)] - dC[3][Index(i,j,k)]);
        }

        /* Psi=Chi=0 at infinity */
        if(Getv("box2_basis2", "ChebExtrema"))
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
          {
            FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)];
            FChi[Index(n1-1,0,k)] = Chi[Index(n1-1,0,k)];
          }
        else // fix: B=0 is not on grid for ChebZeros!!!
        {
          int l;
          double U0;
          double *BM = (double *) calloc(n2, sizeof(double));
          double *line = (double *) calloc(n2, sizeof(double));

          /* obtain BM vector for interpolation along B */
          spec_Basis_times_CoeffMatrix(box->bbox[2],box->bbox[3], n2, BM, 0.0,
                                       cheb_coeffs_fromZeros, cheb_basisfunc);
          for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
          {
            /* find value of Psi at A=1, B=0 */
            get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(n1-1,0,k)] = U0;

            /* find value of Chi at A=1, B=0 */
            get_memline(Chi, line, 2, n1-1,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FChi[Index(n1-1,0,k)] = U0;
          }
          free(line);
          free(BM);
        }
      }
      else if(b==5)  /* in box5 */
      {
        /* values at border are interpolated from box0 */
        double A,B,phi;
        int pl, k_phi;
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[0]->v[vlu->index[0]]; /* values in box0 */
        C = grid->box[0]->v[vlu->index[1]];
        Pcoeffs = grid->box[0]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[0]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[0], P, Pcoeffs);
        spec_Coeffs(grid->box[0], C, Ccoeffs);
        for(pl=0; pl<n1; pl=pl+n1-1)
        {
          int i0;
          int ind=Index(pl,0,0);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);

            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
        for(pl=0; pl<n2; pl=pl+n2-1)
        {
          int i0;
          int ind=Index(0,pl,0);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
        for(pl=0; pl<n3; pl=pl+n3-1)
        {
          int i0;
          int ind=Index(0,0,pl);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[0], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
      }
      else if(b==4)  /* in box4 */
      {
        /* values at border are interpolated from box3 */
        double A,B,phi;
        int pl, k_phi;
        X = box->v[Ind("X")];
        Y = box->v[Ind("Y")];
        Z = box->v[Ind("Z")];
        P = grid->box[3]->v[vlu->index[0]]; /* values in box3 */
        C = grid->box[3]->v[vlu->index[1]];
        Pcoeffs = grid->box[3]->v[Ind("Poisson_temp1")];
        Ccoeffs = grid->box[3]->v[Ind("Poisson_temp2")];
        spec_Coeffs(grid->box[3], P, Pcoeffs);
        spec_Coeffs(grid->box[3], C, Ccoeffs);
        for(pl=0; pl<n1; pl=pl+n1-1)
        {
          int i0;
          int ind=Index(pl,0,0);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
        for(pl=0; pl<n2; pl=pl+n2-1)
        {
          int i0;
          int ind=Index(0,pl,0);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
        for(pl=0; pl<n3; pl=pl+n3-1)
        {
          int i0;
          int ind=Index(0,0,pl);
          phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
          k_phi = grid->box[0]->n3 * phi/(2.0*PI);
          nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
                                    X[ind],Y[ind],Z[ind], 3, k_phi);
          forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
          {
            ind=Index(i,j,k);
            ABphi_of_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
                           
            Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
            Cinterp = spec_interpolate(grid->box[3], Ccoeffs, A,B,phi);
            FPsi[ind] = Psi[ind] - Pinterp;
            FChi[ind] = Chi[ind] - Cinterp;
          }
        }
      }
      else errorexiti("b=%d should be impossible!", b);

    } /* end: else if (Getv("Poisson_grid", "4ABphi_2xyz")) */

  }
}

/* compute A,B,phi from x,y,z */
void ABphi_of_xyz(tBox *box, double *A, double *B, double *phi,
                  double x, double y, double z)
{
  if(*A==0.0) *A+=1e-10;
  if(*A==1.0) *A-=1e-10;
  if(*B==0.0) *B+=1e-10;
  if(*B==1.0) *B-=1e-10;
  *phi = Arg(y,z);   if(*phi<0) *phi = 2.0*PI + *phi;
  XYZ_of_xyz(box, A,B,phi, x,y,z);
  if(*A<0.0) *A=0.0;
  if(*A>1.0) *A=1.0;
  if(*B<0.0) *B=0.0;
  if(*B>1.0) *B=1.0;
}
