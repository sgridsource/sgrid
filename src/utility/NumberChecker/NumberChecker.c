/* NumberChecker.c */
/* Wolfgang Tichy 9/2008 */
/* check if a var contains NANs or INFs or is too large */

#include <math.h>
#include "sgrid.h"
#include "NumberChecker.h"



int NumberChecker_CheckIfFinite(tGrid* grid, char *varname)
{
  double *var;
  int ivar;
  double num=0.0;
  int ijk, ijk_old, b_old;
  tGrid* grid_old;
  double *x; 
  double *y;  
  double *z;
  double *X;
  double *Y;
  double *Z;
  int b;
  int messageflag=0;
  double maxNumber = Getd("NumberChecker_INF");

  /* The stuff inside the brackets below is only needed because 
     errorexit(""); does not work properly!                          */
  /* alternative method that does not find out position of NAN
     incidentally, it makes sure every processor knows about it
  */
  {
    double *sum;
    //tVarList *vl = VLPtrEnable1(grid, varname);

    // TODO: FIX this!!!!!
    //bampi_allreduce_sum(vl, &sum);  

    //if(!finite(*sum)) num=0.1;

    // TODO: FIX this!!!!!
    //free(sum);
num = 0.1;
    
    //vlfree(vl);
    
    if(num==0.0) return 0;
    
num = 0.0;
  }

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    x=box->v[Ind("x")];
    y=box->v[Ind("y")];
    z=box->v[Ind("z")];
    X=box->v[Ind("X")];
    Y=box->v[Ind("Y")];
    Z=box->v[Ind("Z")];

    /* printf("Checking for INF or NAN in %s\n", varname); */
    ivar= Ind(varname);
    var = box->v[ivar]; 
 
    if(var==NULL) 
    {
      printf("pointer to %s is NULL in box%d grid=%p\n",
             VarName(ivar), b, grid); 
      continue;
    }
    else forallpoints(box,ijk)
    {
      if( !finite(var[ijk]) || fabs(var[ijk])>maxNumber ) 
      {
        if(messageflag==0)
          printf("NAN/INF: %s=%g at ijk=%d: x=%g y=%g z=%g "
                 "box%d grid=%p X=%g Y=%g Z=%g\n", 
	         VarName(ivar), var[ijk], ijk, x[ijk], y[ijk], z[ijk],
	         b, grid, X[ijk], Y[ijk], Z[ijk]);
        ijk_old = ijk;
        b_old   = b;
        grid_old= grid;
        messageflag++;
        num++;
      }
      else
      {
        if(messageflag>1)
          printf("NAN/INF: %s=%g til    %d: x=%g y=%g z=%g "
                 "box%d grid=%p X=%g Y=%g Z=%g\n",
                 VarName(ivar), var[ijk_old], ijk_old,
                 x[ijk_old], y[ijk_old], z[ijk_old],
                 b_old, grid_old, X[ijk_old], Y[ijk_old], Z[ijk_old]);
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

   if useful, extend to list of variables in NumberChecker_exitifNAN but this is supposed
   to be cheap enough to be called after every evolution step
*/
int NumberChecker_ExitIfNAN(tGrid* grid)
{
  char *varname = Gets("NumberChecker_exitifNAN");

  if( NumberChecker_CheckIfFinite(grid, varname) !=0 )
  {
    int iterationmax = Geti("iterations");
    double timemax   = Getd("finaltime");
    if(timemax > 0) iterationmax = timemax/grid->dt + 0.5;

    printf("NumberChecker_ExitIfNAN: %s is not finite after evolve!\n"  
	   "  iteration %d\n"
           "  crashtime %g, outdir %s\n",
	   varname, grid->iteration+1, 
	   (grid->iteration+1) * grid->dt, Gets("outdir"));
    /* errorexit("Too bad."); */
    printf("Too bad.\n");
    printf("NumberChecker_ExitIfNAN: setting grid->iteration = iterationmax\n");
    grid->iteration = iterationmax;
  }
  
  return 0;
}
