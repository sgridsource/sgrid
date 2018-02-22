/* quick_out.c */
/* Wolfgang Tichy, 2/14/2018 */
/* some output functions useful for debugging */

#include "sgrid.h"
#include "output.h"



/* add a var and return its index */
int quick_AddVar(tBox *box, char *name)
{
  int inew;
  /* add var name if we don't have it yet */
  inew = IndLax(name);
  if(inew<0) AddVarToGrid(box->grid, name, "", "added for quick_Var_output");
  inew = Ind(name);
  return inew;
}


/* output a var in a box, with normal write_grid, but some output pars,
   time and iteration are changed */
void quick_Vars_output(tBox *box, char *names, double fake_t, int fake_i)
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

  /* output all right now, but set time for output to fake_t */
  Sets("0doutput", names);
  Sets("1doutput", names);
  Sets("2doutput", names);
  Sets("3doutput", names);
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
}

/* same as quick_Vars_output, but first copy an array into var in box */
void quick_Array_output(tBox *box, double *Ar, char *name,
                        double fake_t, int fake_i)
{
  int i, inew, ison=0;

  inew = quick_AddVar(box, name);
  if(box->v[inew]) ison=1;
  else             enablevar_inbox(box, inew);

  /* copy ar into var */
  forallpoints(box, i) box->v[inew][i] = Ar[i];

  /* output the array, and then disable the var */
  quick_Vars_output(box, name, fake_t, fake_i);
  if(!ison) disablevar_inbox(box, inew);
}

/* same as quick_Vars_output, but for an entire varlist */
void quick_VarList_output(tBox *box, tVarList *vl, double fake_t, int fake_i)
{
  int nl = 3;
  char *names = calloc(nl, sizeof(char));
  int vind;

  for(vind=0; vind<vl->n; vind++)
  {
    char *varname = VarName(vl->index[vind]);
    nl += strlen(varname)+1;
    names = realloc(names, nl);
    if(vind>0) strcat(names, " ");
    strcat(names, varname);
  }
  quick_Vars_output(box, names, fake_t, fake_i);
  free(names);
}
