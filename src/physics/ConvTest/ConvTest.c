/* ConvTest.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "ConvTest.h"



/* evolution step for testing */
void ConvTest_evo(tVarList *vlunew, 
		  tVarList *vlupre, double dt, tVarList *vlucur)
{
  tGrid *grid = vlucur->grid;
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *ucur = vlldataptr(vlucur, box, 0);
    double *upre = vlldataptr(vlupre, box, 0);
    double *unew = vlldataptr(vlunew, box, 0);
    int i;
    double *u1 = box->v[Ind("ConvTest_u1")];
    double *u2 = box->v[Ind("ConvTest_u2")];
    double *u3 = box->v[Ind("ConvTest_u3")];
    double *x = box->v[Ind("x")];
    double *y = box->v[Ind("y")];
    double *z = box->v[Ind("z")];

    cart_partials(box, ucur, u1, u2, u3);

    forallpoints(box, i)
    {
      double r = sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      double urhs = ( u1[i]*x[i] + u2[i]*y[i] + u3[i]*z[i] )/r;
      
      if(dt!=0.0)
        unew[i] = upre[i] + dt * urhs;
      else
        unew[i] = urhs;
    }
    /* set boundary to zero */
    {
      int j,k;
      int n1=box->n1;
      int n2=box->n2;
      int n3=box->n3;

      for(k=0; k<n3; k++)
        for(j=0; j<n2; j++)
          unew[Index(n1-1,j,k)] = 0;
    }
  }
}


/* initialize test */
int ConvTest_startup(tGrid *grid) 
{
  tVarList *vlu;
  int b;
  
  enablevar(grid, Ind("ConvTest_u"));
  enablevar(grid, Ind("ConvTest_u1"));
  enablevar(grid, Ind("ConvTest_u2"));
  enablevar(grid, Ind("ConvTest_u3"));
  enablevar(grid, Ind("ConvTest_err"));
  
  vlu = vlalloc(grid);
  vlpush(vlu, Ind("ConvTest_u"));

  for (b = 0; b < grid->nboxes; b++)
  {  
    tBox *box = grid->box[b];
    int i,j,k;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *X = box->v[Ind("X")];
    double *Y = box->v[Ind("Y")];
    double *Z = box->v[Ind("Z")];
    double *x = box->v[Ind("x")];
    double *y = box->v[Ind("y")];
    double *z = box->v[Ind("z")];

    double *u = box->v[Ind("ConvTest_u")];
    double *u1 = box->v[Ind("ConvTest_u1")];
    double *u2 = box->v[Ind("ConvTest_u2")];
    double *u3 = box->v[Ind("ConvTest_u3")];

    printf("Testing convergence with advection eqn:\n");

    forallpoints(box, i)
    {
     double xs=x[i];
     double ys=y[i];
     double zs=z[i];
     double rs, f;

     rs = sqrt(xs*xs + ys*ys + zs*zs) + grid->time;
     f  = exp( -10.0*( rs*rs - 2*1.2*rs*xs/X[i] + 1.2*1.2 ) );
     
     u[i] = f;
    }
    cart_partials(box, u, u1, u2, u3);

    for(k=0; k<n3; k++)
      for(j=0; j<n2; j++)
        u[Index(n1-1,j,k)] = 0;
  }
  evolve_vlregister(vlu);
  evolve_rhsregister(ConvTest_evo);
  return 0;
}


/* compute normalized absolute error in ANALYSIS */
int ConvTest_analyze(tGrid *grid)
{
  int b;
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *X = box->v[Ind("X")];
    double *Y = box->v[Ind("Y")];
    double *Z = box->v[Ind("Z")];
    double *x = box->v[Ind("x")];
    double *y = box->v[Ind("y")];
    double *z = box->v[Ind("z")];
    double *u = box->v[Ind("ConvTest_u")];
    double *err = box->v[Ind("ConvTest_err")];
    int i;

    forallpoints(box, i)
    {
      double xs=x[i];
      double ys=y[i];
      double zs=z[i];
      double rs, f;

      rs = sqrt(xs*xs + ys*ys + zs*zs) + grid->time;
      f  = exp( -10.0*( rs*rs - 2*1.2*rs*xs/X[i] + 1.2*1.2 ) );
      //printf("f=%g t=%g| ", f, grid->time);
      /* computer error */
      err[i] = u[i] - f;
      //err[i] = f;
    }
  }
  return 0;
}


int ConvTest_filter(tGrid *grid)
{
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int j,k;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *u = box->v[Ind("ConvTest_u")];

    /* use filters */
    spec_filter1(box, 1, u);
    spec_filter1(box, 2, u);
    spec_filter1(box, 3, u);

    for(k=0; k<n3; k++)
      for(j=0; j<n2; j++)
        u[Index(n1-1,j,k)] = 0;
  }
  return 0;
}
