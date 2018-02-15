/* quick_out.c */
/* Wolfgang Tichy, 2/14/2018 */
/* some output functions useful for debugging */

#include "sgrid.h"
#include "output.h"



/* add a var and enable in a box, return its index */
int quick_AddVar(tBox *box, char *name)
{
  int inew;
  /* add var name if we don't have it yet, and enable it */
  inew = IndLax(name);
  if(inew<0) AddVarToGrid(box->grid, name, "", "added for quick_Var_output");
  inew = Ind(name);
  enablevar_inbox(box, inew);
  return inew;
}


/* output a var in a box, with normal write_grid, but some output pars
   and time changed, it can also exit after writing if stop=1 */
void quick_Var_output(tBox *box, char *name,
                      double fake_t, int fake_i, int stop)
{
  tGrid *grid = box->grid;
  double t = grid->time;
  int   it = grid->iteration;
  char *output0 = strdup(Gets("0doutput"));
  char *output1 = strdup(Gets("1doutput"));
  char *output2 = strdup(Gets("2doutput"));
  char *output3 = strdup(Gets("3doutput"));
  char *format3 = strdup(Gets("3dformat"));
  char *outiter0 = strdup(Gets("0doutiter"));
  char *outiter1 = strdup(Gets("1doutiter"));
  char *outiter2 = strdup(Gets("2doutiter"));
  char *outiter3 = strdup(Gets("3doutiter"));

  quick_AddVar(box, name);

  /* output all right now, but set time for output to fake_t */
  Sets("0doutput", name);
  Sets("1doutput", name);
  Sets("2doutput", name);
  Sets("3doutput", name);
  Sets("3dformat", "vtk text double arrange_as_1d");
  Sets("0doutiter","1");
  Sets("1doutiter","1");
  Sets("2doutiter","1");
  Sets("3doutiter","1");
  grid->time = fake_t;
  grid->iteration = fake_i;
  write_grid(grid);
  /* restore time and some pars */
  grid->time = t;
  grid->iteration = it;
  Sets("0doutput", output0);
  Sets("1doutput", output1);
  Sets("2doutput", output2);
  Sets("3doutput", output3);
  Sets("3dformat", format3);
  Sets("0doutiter", outiter0);
  Sets("1doutiter", outiter1);
  Sets("2doutiter", outiter2);
  Sets("3doutiter", outiter3);
  free(output0);
  free(output1);
  free(output2);
  free(output3);
  free(format3);
  free(outiter0);
  free(outiter1);
  free(outiter2);
  free(outiter3);
  if(stop) exit(1);
}

/* same as quick_Var_output, but first copy an array into var in box */
void quick_Array_output(tBox *box, double *Ar, char *name,
                        double fake_t, int fake_i, int stop)
{
  int i;
  int inew = quick_AddVar(box, name);

  /* copy ar into var */
  forallpoints(box, i) box->v[inew][i] = Ar[i];

  /* output the arraay */
  quick_Var_output(box, name, fake_t, fake_i, stop);
}
