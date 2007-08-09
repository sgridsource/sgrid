/* Poisson.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "Poisson.h"

#define Power pow

// remove this
int linSol(tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	     int imax, double tol, double *res,
	     void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));

/* initialize Poisson */
int Poisson_startup(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
  double sigmax    = Getd("Poisson_sigmax");
  double sigmay    = Getd("Poisson_sigmay");
  double sigmaz    = Getd("Poisson_sigmaz");
  double x0        = Getd("Poisson_x0");
  double y0        = Getd("Poisson_y0");
  double z0        = Getd("Poisson_z0");
*/
  printf("Initializing Poisson:\n");

  /* set boundary information for ScalarWave evolution: 
     farlimit, falloff, propagation speed   */
  VarNameSetBoundaryInfo("Poisson_Psi", 0, 1, 1.0);

  /* enable all vars */
  enablevar(grid, Ind("Poisson_Psi"));
  enablevar(grid, Ind("Poisson_Psix"));
  enablevar(grid, Ind("Poisson_Psixx"));
  enablevar(grid, Ind("Poisson_Err_Psi"));
  
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
    double *Psi = box->v[Ind("Poisson_Psi")];

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
      /* set Psi */
      Psi[i] = 2.0;
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

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vlFu = vlalloc(grid);
  vluDerivs= vlalloc(grid);
  
  /* add Poisson_Psi to vlu */
  vlpush(vlu, Ind("Poisson_Psi"));

  /* add derivs of Poisson_Psi to vluDerivs */
  vlpush(vluDerivs, Ind("Poisson_Psix"));
  vlpush(vluDerivs, Ind("Poisson_Psixx"));

  /* enable vlu, vluDerivs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 

  /* now duplicate vlu for result of F(u) */
  vlFu = AddDuplicateEnable(vlu,  "res");

  /* now duplicate vlu and vluDerivs for linarized Eqs. */
  vldu       = AddDuplicateEnable(vlu,  "l");
  vlr        = AddDuplicateEnable(vlFu, "l");
  vlduDerivs = AddDuplicateEnable(vluDerivs, "l");

  /* call Newton solver */
  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, 0,
         itmax, tol, &normresnonlin, 1,
         bicgstab, Precon_I, vldu, vlr, vlduDerivs, vlu,
         linSolver_itmax, linSolver_tolFac);

  Newton(F_Poisson, J_Poisson, vlu, vlFu, vluDerivs, 0,
         itmax, tol, &normresnonlin, 1,
         linSol, Precon_I, vldu, vlr, vlduDerivs, vlu,
         linSolver_itmax, linSolver_tolFac);

  /* free varlists */     
  VLDisableFree(vldu);
  VLDisableFree(vlduDerivs);
  VLDisableFree(vlr);     
  vlfree(vlu);
  vlfree(vluDerivs);
  vlfree(vlFu);
        
  return 0;
}


/* compute absolute error in ANALYSIS */
int Poisson_analyze(tGrid *grid)
{
  int b;
/*
  double A         = Getd("Poisson_A");
  double sigmax    = Getd("Poisson_sigmax");
  double sigmay    = Getd("Poisson_sigmay");
  double sigmaz    = Getd("Poisson_sigmaz");
  double x0        = Getd("Poisson_x0");
  double y0        = Getd("Poisson_y0");
  double z0        = Getd("Poisson_z0");
*/
  printf("Poisson: computing absolute error in derivatives\n");
  
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

      PsiErr[i] = Psi[i]-0;
    }
  }
  return 0;
}


/* evaluate Poisson eqn for vlu */
void F_Poisson(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2)
{
  tGrid *grid = vlu->grid;
  int i,j,b;
  	
  for(j = 0; j < vlu->n; j++)
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *Fu = box->v[vlFu->index[j]];
      double *u  = box->v[vlu->index[j]];

//    /* compute the derivs */
//    if(Getv("Poisson_useDD", "yes"))
//      allDerivsOf_S(box, Ind("Poisson_u"), Ind("Poisson_Err_dux"), 
//                    Ind("Poisson_Err_dduxx"));
//    else
//      FirstAndSecondDerivsOf_S(box, Ind("Poisson_u"), Ind("Poisson_Err_dux"),
//                               Ind("Poisson_Err_dduxx"));

      forallpoints(box, i)  Fu[i] = u[i]*u[i] - 9.0;
    }
}

/* evaluate linearized Laplace */
void J_Poisson(tVarList *vlJdu, tVarList *vldu,
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
      double *u   = box->v[vlu->index[j]];

      forallpoints(box, i)  Jdu[i] = 2.0*u[i]*du[i];
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

int linSol(tVarList *x, tVarList *b, tVarList *r, tVarList *c1,tVarList *c2,
	     int imax, double tol, double *res,
	     void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int i,j,bi;
  	
//  lop(r, , c);

  for(j = 0; j < b->n; j++)
    forallboxes(grid, bi)
    {
      tBox *box = grid->box[bi];
      double *px = box->v[x->index[j]];
      double *pb = box->v[b->index[j]];
      double *pc2= box->v[c2->index[j]];

      forallpoints(box, i)  px[i] = pb[i]/(2*pc2[i]);
    }
  return 666;
}
