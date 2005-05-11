/* Dissipation.c */
/* Wolfgang Tichy 11.5.2005 */

#include "sgrid.h"


/* add diffusion terms to all variables: 
   u_n = u_n + dt * [ (DissFac) dt^DissOrder \nabla^2 u_c ] */

/* Note: the result will be written into vlu_n  */
void addDissipation(tVarList *vlu_n, tVarList *vlu_p, double dt, 
                    tVarList *vlu_c)
{
  tGrid *grid = vlu_c->grid;
  double Lap_u;
  double *u_n;
  double *u_c;
  int bi, j, ijk;
  double DissFac     = Getd("evolve_DissipationFactor");
  double DissOrder   = Getd("evolve_Dissipation_dt_Order");
  double dt_to_Order = pow(grid->dt, DissOrder);
  int dt_nonzero;
    
  if(DissFac<=0) return;

  if(dt!=0.0) dt_nonzero=1;
  else        dt_nonzero=0;

  if( (grid->iteration == 0) )
  {
    printf("addDissipation: Numerical Dissipation is active.\n");
    printf("evolve_DissipationFactor = %e\n", DissFac);
    printf("evolve_Dissipation_dt_Order = %e\n", DissOrder);
    printf("dt_nonzero:%d\n", dt_nonzero);
  }
  
 
  forallboxes(grid,bi)  
  {
    tBox *box=grid->box[bi];
    /* abuse ADMvars_dgxxi and ADMvars_ddgxxij as temp storage */
    double *dux = box->v[Ind("ADMvars_dgxxx")];
    double *duy = box->v[Ind("ADMvars_dgxxy")];
    double *duz = box->v[Ind("ADMvars_dgxxz")];
    double *duxx = box->v[Ind("ADMvars_ddgxxxx")];
    double *duxy = box->v[Ind("ADMvars_ddgxxxy")];
    double *duxz = box->v[Ind("ADMvars_ddgxxxz")];
    double *duyy = box->v[Ind("ADMvars_ddgxxyy")];
    double *duyz = box->v[Ind("ADMvars_ddgxxyz")];
    double *duzz = box->v[Ind("ADMvars_ddgxxzz")];
      
    /* loop over all vars */
    for(j = 0; j < vlu_c->n; j++)
    {
       u_c = box->v[vlu_c->index[j]];
       u_n = box->v[vlu_n->index[j]];

       cart_partials(box, u_c, dux, duy, duz);
       cart_partials(box, dux, duxx, duxy, duxz);
       cart_partials(box, duy, duxy, duyy, duyz);
       cart_partials(box, duz, duxz, duyz, duzz);
        
       forallpoints(box,ijk)   
       {
          /* compute \nabla^2 u_c (with metric ( 1   0   0 
                                                 0   1   0 
                                                 0   0   1  )  ) */
          Lap_u = duxx[ijk] + duyy[ijk] + duzz[ijk];

          /* add Dissipation terms to u_n.   Note: Lap_u = \nabla^2 u_c */
          if(dt_nonzero)
            u_n[ijk] += dt * DissFac * dt_to_Order * Lap_u;
          else
            u_n[ijk] += DissFac * dt_to_Order * Lap_u;
       }
    }
  } /* end of box loop */
}
