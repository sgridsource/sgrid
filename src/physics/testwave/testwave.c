/* testwave.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "testwave.h"



/* evolution step for testing */
void testwave_evo(tVarList *vlunew, 
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
    double *kcur = vlldataptr(vlucur, box, 1);
    double *kpre = vlldataptr(vlupre, box, 1);
    double *knew = vlldataptr(vlunew, box, 1);
    int i;
    double *u11 = box->v[Ind("testwave_u11")];
    double *u22 = box->v[Ind("testwave_u22")];
    double *u33 = box->v[Ind("testwave_u33")];
    double *dum = box->v[Ind("testwave_dum")];

    //spec_allDerivs(box, ucur, dum, dum, dum, u11, dum, dum, u22, dum, u33);
    cart_partials(box, ucur, dum, u22, u33);
    cart_partials(box, dum, u11, u22, u33);

    forallpoints(box, i)
    {
      double krhs = u11[i] + 0*u22[i] + 0*u33[i]
                    - 0.0 * ucur[i]*ucur[i]*ucur[i]/(1.0+ucur[i]*ucur[i]);
      double urhs = kcur[i];
      
      if(dt!=0.0)
      {
        knew[i] = kpre[i] + dt * krhs;
        unew[i] = upre[i] + dt * urhs;
        //knew[i] = kpre[i];
        //unew[i] = upre[i] - dt * u11[i];
      }
      else
      {
        knew[i] = krhs;
        unew[i] = urhs;
        //knew[i] = 0;
        //unew[i] = -u11[i];
      }
    }
    /* set boundary to zero */
    {
      int i,j,k;
      int n1=box->n1;
      int n2=box->n2;
      int n3=box->n3;

      for(k=0; k<n3; k++)
        for(j=0; j<n2; j++)
        {
          unew[Index(0,j,k)] = unew[Index(n1-1,j,k)] = 0;
          knew[Index(0,j,k)] = knew[Index(n1-1,j,k)] = 0;
          //unew[Index(0,j,k)] = 0;
        }

      for(k=0; k<n3; k++)
        for(i=0; i<n1; i++)
        {
          unew[Index(i,0,k)] = unew[Index(i,n2-1,k)] = 0;
          knew[Index(i,0,k)] = knew[Index(i,n2-1,k)] = 0;
        }

      for(j=0; j<n2; j++)
        for(i=0; i<n1; i++)
        {
          unew[Index(i,j,0)] = unew[Index(i,j,n3-1)] = 0;
          knew[Index(i,j,0)] = knew[Index(i,j,n3-1)] = 0;
        }

//      cheb_coeffs_fromExtrema(dum, unew+Index(0,1,1), n1-1);
//      cheb_filter(dum, 2*(n1-1)/3, n1-1);
//      cheb_eval_onExtrema(dum, unew+Index(0,1,1), n1-1);

//      cheb_coeffs_fromExtrema(dum, knew+Index(0,1,1), n1-1);
//      cheb_filter(dum, 2*(n1-1)/3, n1-1);
//      cheb_eval_onExtrema(dum, knew+Index(0,1,1), n1-1);

//      cheb_coeffs_fromExtrema(u11+Index(0,1,1), unew+Index(0,1,1), n1-1);
    }
  }
}


/* initialize test */
int testwave_startup(tGrid *grid) 
{
  tVarList *vlu;
  int b;
  
  enablevar(grid, Ind("testwave_u"));
  enablevar(grid, Ind("testwave_k"));
  enablevar(grid, Ind("testwave_u11"));
  enablevar(grid, Ind("testwave_u22"));
  enablevar(grid, Ind("testwave_u33"));
  enablevar(grid, Ind("testwave_dum"));

  vlu = vlalloc(grid);
  vlpush(vlu, Ind("testwave_u"));
  vlpush(vlu, Ind("testwave_k"));

  for (b = 0; b < grid->nboxes; b++)
  {  
    tBox *box = grid->box[b];
    int i,j,k, ijk;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *pu = box->v[Ind("testwave_u")];
    double *pk = box->v[Ind("testwave_k")];

    double *u11 = box->v[Ind("testwave_u11")];
    double *u22 = box->v[Ind("testwave_u22")];
    double *u33 = box->v[Ind("testwave_u33")];
    double *dum = box->v[Ind("testwave_dum")];

    printf("Testing wave eqn:\n");
    
    forallijk(i,j,k)
    {
     double X=-cos(i*PI/(n1-1))-0.2;
     double Y=-cos(j*PI/(n2-1));
     double Z=-cos(k*PI/(n3-1));
     double f=exp( -10*(X*X + Y*Y + Z*Z) );

     ijk=Index(i,j,k);
     pu[ijk] = f;
     pk[ijk] = -(-2*10*X*f); // = - d/dx u
     //pk[ijk] = -2*10*f + (2*10*X)*(2*10*X)*f;
     //pk[ijk] = -2*10*X*f;
     
     //pu[ijk] = 1;
     //pk[ijk] = 0;
    }
    spec_allDerivs(box, pu, dum, dum, dum,
                   u11, dum, dum, u22, dum, u33 );

//    cheb_coeffs_fromExtrema(u11, pu+Index(0,1,1), n1-1);
//    cheb_deriv(-1,1, u11, dum, n1-1);
//    cheb_eval_onExtrema(dum, u11+Index(0,1,1), n1-1);

//    cheb_coeffs_fromExtrema(dum, pu+Index(0,1,1), n1-1);
//    cheb_eval_onExtrema(dum, u11+Index(0,1,1), n1-1);

//    for(i=0; i<n1; i++)
//    {
//      u11[Index(i,1,1)] =
//      cheb_eval(-1,1, dum, n1-1, -cos(i*PI/(n1-1)) );
//    }

      for(k=0; k<n3; k++)
        for(j=0; j<n2; j++)
        {
          pu[Index(0,j,k)] = pu[Index(n1-1,j,k)] = 0;
          pk[Index(0,j,k)] = pk[Index(n1-1,j,k)] = 0;
        }
      //cheb_coeffs_fromExtrema(u11+Index(0,1,1), pu+Index(0,1,1), n1-1);

  }
  evolve_vlregister(vlu);
  evolve_rhsregister(testwave_evo);
  return 0;
}


/* compute normalized absolute error in ANALYSIS */
int testwave_analyze(tGrid *grid)
{
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *u = box->v[Ind("testwave_u")];
    int i;
    double xxx;

    /* computer error */
    forallpoints(box, i) 
      xxx = fabs(u[i]/exp(-grid->time) - 1);
  }
  return 0;
}


int testwave_filter(tGrid *grid)
{
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int i,j,k, ijk;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *pu = box->v[Ind("testwave_u")];
    double *pk = box->v[Ind("testwave_k")];
    double *dum = box->v[Ind("testwave_dum")];

    pu[Index(0,1,1)] = pu[Index(n1-1,1,1)] = 0;
    pk[Index(0,1,1)] = pk[Index(n1-1,1,1)] = 0;    

//    cheb_coeffs_fromExtrema(dum, pu+Index(0,1,1), n1-1);
//    cheb_filter(dum, 2*(n1-1)/3, n1-1);
//    cheb_eval_onExtrema(dum, pu+Index(0,1,1), n1-1);

//    cheb_coeffs_fromExtrema(dum, pk+Index(0,1,1), n1-1);
//    cheb_filter(dum, 2*(n1-1)/3, n1-1);
//    cheb_eval_onExtrema(dum, pk+Index(0,1,1), n1-1);
  }
  return 0;
}
