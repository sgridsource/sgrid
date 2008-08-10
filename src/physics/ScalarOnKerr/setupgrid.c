/* grid_setup.c */
/* Wolfgang Tichy 2008 */


#include "sgrid.h"
#include "ScalarOnKerr.h"

#define STRLEN 1000

/* setup initial boxsizes */
int ScalarOnKerr_setup_boxes(tGrid *g)
{
  tGrid *grid;
  int nboxes=Geti("nboxes");
  int b;
  double Xshift;

  if(Getv("ScalarOnKerr_overlap_shells", "no")) return 0;

  printf("ScalarOnKerr_setup_boxes: setting box sizes used ...\n");

  /* make a dummy grid with 3 vars */
  grid=make_empty_grid(globalnvariables, 1);
  set_BoxStructures_fromPars(grid, 1);

  /* shift shell radii */
  for(Xshift=0.0, b=1; b<nboxes; b++)
  {
    tBox *box = grid->box[b];
    tBox *lbox= grid->box[b-1];
    int n1 = box->n1;
    int ln1=lbox->n1;
    double *X  = box->v[Ind("X")];
    double *lX =lbox->v[Ind("X")];
    double dX,dlX;
    char str[STRLEN];

    dX  =  X[1]-X[0];
    dlX = lX[ln1-1]-lX[ln1-2];

    if(dX>dlX) Xshift+=lX[ln1-1]-X[1];
    else       Xshift+=lX[ln1-2]-X[0];

    snprintf(str, STRLEN, "box%d_min1", b);  Setd(str, X[0]+Xshift);
    snprintf(str, STRLEN, "box%d_max1", b);  Setd(str, X[n1-1]+Xshift);
  }

  /* free dummy grid */
  free_grid(grid);

  return 0;
}
