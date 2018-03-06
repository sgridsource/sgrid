/* TestDerivs.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "TestDerivs.h"

#define Power pow

/* initialize TestDerivs */
int TestDerivs_startup(tGrid *grid)
{
  int b;
  double A         = Getd("TestDerivs_A");
  double sigmax    = Getd("TestDerivs_sigmax");
  double sigmay    = Getd("TestDerivs_sigmay");
  double sigmaz    = Getd("TestDerivs_sigmaz");
  double x0        = Getd("TestDerivs_x0");
  double y0        = Getd("TestDerivs_y0");
  double z0        = Getd("TestDerivs_z0");
  
  printf("Initializing TestDerivs:\n");

//  //populate_bfaces(grid);
//  double q[3]={65.623508,8.236062,2.648371};
//  double X[3]={1.,1.,1.};
//  b_XYZ_of_xyz(grid,&X[0],&X[1],&X[2],q[0],q[1],q[2]);
//  abort();

  
  /* enable all vars */
  enablevar(grid, Ind("TestDerivs_u"));
  enablevar(grid, Ind("TestDerivs_Err_dux"));
  enablevar(grid, Ind("TestDerivs_Err_dduxx"));
  
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
    double *u = box->v[Ind("TestDerivs_u")];

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
      /* Gaussian */
      u[i] = A*exp(-0.5*( 
               (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
             + (z-z0)*(z-z0)/(sigmaz*sigmaz)   ) );
    }
  }

//convert_grid_to_fd_onesidedBC(grid);
//forallboxes(grid,b)
//{
//tBox *box = grid->box[b];
//int n2 = box->n2;
//int n3 = box->n3;
////initdiffmatrix(box->bbox[2], box->bbox[3], box->D2, box->DD2, n2,
////               fd2_coeffs, fd2_deriv_periodic, fd2_eval);
//initdiffmatrix(box->bbox[4], box->bbox[5], box->D3, box->DD3, n3,
//               fd2_coeffs, fd2_deriv_periodic, fd2_eval);
//}

  return 0;
}


/* compute absolute error in ANALYSIS */
int TestDerivs_analyze(tGrid *grid)
{
  int b;
  double A         = Getd("TestDerivs_A");
  double sigmax    = Getd("TestDerivs_sigmax");
  double sigmay    = Getd("TestDerivs_sigmay");
  double sigmaz    = Getd("TestDerivs_sigmaz");
  double x0        = Getd("TestDerivs_x0");
  double y0        = Getd("TestDerivs_y0");
  double z0        = Getd("TestDerivs_z0");
  
  printf("TestDerivs: computing absolute error in derivatives\n");
  
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
    double *u  = box->v[Ind("TestDerivs_u")];
    double *ux = box->v[Ind("TestDerivs_Err_dux")];
    double *uy = box->v[Ind("TestDerivs_Err_duy")];
    double *uz = box->v[Ind("TestDerivs_Err_duz")];
    double *uxx= box->v[Ind("TestDerivs_Err_dduxx")];
    double *uxy= box->v[Ind("TestDerivs_Err_dduxx")+1];
    double *uxz= box->v[Ind("TestDerivs_Err_dduxx")+2];
    double *uyy= box->v[Ind("TestDerivs_Err_dduxx")+3];
    double *uyz= box->v[Ind("TestDerivs_Err_dduxx")+4];
    double *uzz= box->v[Ind("TestDerivs_Err_dduxx")+5];

    /* compute the derivs */
    if(Getv("TestDerivs_useDD", "yes"))
      allDerivsOf_S(box, Ind("TestDerivs_u"), Ind("TestDerivs_Err_dux"), 
                    Ind("TestDerivs_Err_dduxx"));
    else
      FirstAndSecondDerivsOf_S(box, Ind("TestDerivs_u"), Ind("TestDerivs_Err_dux"),
                               Ind("TestDerivs_Err_dduxx"));

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
      if(fabs(x)+fabs(y)+fabs(z)>1e299) continue; /* give up if x,y,z is inf */
      /*
      u[i] = A*exp(-0.5*( 
               (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
             + (z-z0)*(z-z0)/(sigmaz*sigmaz)   ) );

      $Assuumptions = sigmax>0 && sigmay>0 && sigmaz>0
      u=A*Exp[-( (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
                +(z-z0)*(z-z0)/(sigmaz*sigmaz))/2]

      CForm[Simplify[D[u,x]/u]]
       = (-x + x0)/Power(sigmax,2)
      CForm[Simplify[D[u,y]/u]]
       = (-y + y0)/Power(sigmay,2)
      CForm[Simplify[D[u,z]/u]]
       = (-z + z0)/Power(sigmaz,2)
       
      CForm[Simplify[D[u,x,x]/u]]
       = (-Power(sigmax,2) + Power(x - x0,2))/Power(sigmax,4)
      CForm[Simplify[D[u,x,y]/u]]
       = ((x - x0)*(y - y0))/(Power(sigmax,2)*Power(sigmay,2))
      CForm[Simplify[D[u,x,z]/u]]
       = ((x - x0)*(z - z0))/(Power(sigmax,2)*Power(sigmaz,2))
      CForm[Simplify[D[u,y,y]/u]]
       = (-Power(sigmay,2) + Power(y - y0,2))/Power(sigmay,4)
      CForm[Simplify[D[u,y,z]/u]]
       = ((y - y0)*(z - z0))/(Power(sigmay,2)*Power(sigmaz,2))
      CForm[Simplify[D[u,z,z]/u]]
       = (-Power(sigmaz,2) + Power(z - z0,2))/Power(sigmaz,4)
      */
//printf("z0=%f ",z0);
      ux[i] -= (-x + x0)/Power(sigmax,2) * u[i];
      uy[i] -= (-y + y0)/Power(sigmay,2) * u[i];
      uz[i] -= (-z + z0)/Power(sigmaz,2) * u[i];
      uxx[i]-= (-Power(sigmax,2) + Power(x - x0,2))/Power(sigmax,4) * u[i];
      uxy[i]-= ((x - x0)*(y - y0))/(Power(sigmax,2)*Power(sigmay,2)) * u[i];
      uxz[i]-= ((x - x0)*(z - z0))/(Power(sigmax,2)*Power(sigmaz,2)) * u[i];
      uyy[i]-= (-Power(sigmay,2) + Power(y - y0,2))/Power(sigmay,4) * u[i];
      uyz[i]-= ((y - y0)*(z - z0))/(Power(sigmay,2)*Power(sigmaz,2)) * u[i];
      uzz[i]-= (-Power(sigmaz,2) + Power(z - z0,2))/Power(sigmaz,4) * u[i];
    }
/*
spec_Deriv2(box, 1, u, ux);
spec_Deriv2(box, 2, u, uy);
spec_Deriv2(box, 3, u, uz);
*/
  }
  return 0;
}
