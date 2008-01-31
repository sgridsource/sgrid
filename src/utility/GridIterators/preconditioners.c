/* preconditioners.c */
/* Wolfgang Tichy 1/2008 */


#include "sgrid.h"
#include "GridIterators.h"


/* Einheitsmatrix als Preconditioner (good for testing) */ 
void Preconditioner_I(tVarList *vlJdu, tVarList *vldu,
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
