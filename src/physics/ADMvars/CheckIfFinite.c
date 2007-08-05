/* CheckIfFinite.c */
/* Wolfgang Tichy 3/2003 */
/* check if ham contains NANs or INFs */

#include <math.h>
#include "sgrid.h"
#include "ADMvars.h"




int CheckIfFinite(tGrid* grid, char *varname)
{
  double *var;
  int ivar;
  double num=0.0;
  int ccc, ccc_old;
  double *x; 
  double *y;  
  double *z;
  int b;
  int messageflag=0;

  /* The stuff inside the brackets below is only needed because 
     errorexit(""); does not work properly!                          */
  /* alternative method that does not find out position of NAN
     incidentally, it makes sure every processor knows about it
  */
  {
    double *sum;
    tVarList *vl = VLPtrEnable1(grid, varname);

    // TODO: FIX this!!!!!
    //bampi_allreduce_sum(vl, &sum);  

    if(!finite(*sum)) num=0.1;

    // TODO: FIX this!!!!!
    //free(sum);
num = 0.1;
    
    vlfree(vl);
    
    if(num==0.0) return 0;
    
num = 0.0;
  }

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];

    /* printf("Checking for INF or NAN in %s\n", varname); */
    ivar= Ind(varname);
    var = box->v[ivar]; 
 
    forallpoints(box,ccc)
    {
      if(var==NULL) 
      {
        x=box->v[Ind("x")];
        y=box->v[Ind("y")];
        z=box->v[Ind("z")];
        printf("pointer to %s is NULL at ccc=%d:  x=%f y=%f z=%f\n",
	     VarName(ivar),ccc,x[ccc],y[ccc],z[ccc]); 
        continue;
      }
        
      if( !finite(var[ccc]) ) 
      {
        if(messageflag==0)
        {
	x=box->v[Ind("x")];
	y=box->v[Ind("y")];
	z=box->v[Ind("z")];
	printf("NAN/INF: %s=%g at ccc=%d: box%d x=%g y=%g z=%g\n", 
	       VarName(ivar), var[ccc], ccc, b, x[ccc], y[ccc], z[ccc]);
        }
        ccc_old=ccc;
        messageflag++;
        num++;
      }
      else
      {
        if(messageflag>1)
          printf("NAN/INF: %s=%g til    %d: box%d x=%g y=%g z=%g\n",
                 VarName(ivar), var[ccc_old], ccc_old, b,
                 x[ccc_old], y[ccc_old], z[ccc_old]);
                       
        messageflag=0;
      }
    }
  }
  if(num>0)
  {
    printf("%d NAN/INFs were detected at time %g, iteration %d.\n",
           ((int) num), grid->time, grid->iteration);
    fflush(stdout);

    if(num<1) return -1;
  }

  return ((int) num);
}




/* check for NANs, alternative version, to be registered 
   currently registered in POST_EVOLVE so that NANs are not output
   register in PRE_EVOLVE so that the previous output step has completed

   if useful, extend to list of variables in ADMvars_exitifNAN but this is supposed
   to be cheap enough to be called after every evolution step
*/
int ExitIfNAN(tGrid* grid)
{
  char *varname = Gets("ADMvars_exitifNAN");

  if( CheckIfFinite(grid, varname) !=0 )
  {
    printf("ExitIfNAN: %s is not finite after evolve!\n"  
	   "  iteration %d\n"
           "  crashtime %f, outdir %s\n",
	   varname, grid->iteration+1, 
	   (grid->iteration+1) * grid->dt, Gets("outdir"));
    errorexit("Too bad.");
  }
  
  return 0;
}
