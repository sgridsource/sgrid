/* output.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 12/99 */

#include "sgrid.h"
#include "output.h"




/* decide whether it is time for output 
   return 0 if not
   return number of output, 1 for t=0, 2 for next, ...
   tricky bit is that we have two criteria, di and dt, and dt could mean 
   also 'output after this time is exceeded' if dt and the timestep don't match
*/
int timeforoutput_di_dt(tGrid *grid, int di, double dt) 
{
  /* time for output based on number of iterations */
  if (di > 0 && grid->iteration % di == 0) 
    return grid->iteration/di + 1;

  /* time for output based on time interval, assumes t >= 0 */
  if (dt > 0) {
    int i = (grid->time + dequaleps)/dt;
    if (dequal(grid->time-i*dt, 0))
      return i + 1;
  }
    
  /* not time for output */
  return 0;
}




/* decide whether it is time for output in any format */
int timeforoutput_any(tGrid *grid)
{
  return 
    timeforoutput_di_dt(grid, Geti("0doutiter"), Getd("0douttime")) ||
    timeforoutput_di_dt(grid, Geti("1doutiter"), Getd("1douttime")) ||
    timeforoutput_di_dt(grid, Geti("2doutiter"), Getd("2douttime")) ||
    timeforoutput_di_dt(grid, Geti("3doutiter"), Getd("3douttime"));
}




/* decide whether it is time for output for a particular variable 
   the variable has to be present in any one of the output parameters, and
   the iteration count has to be right
*/
#define NOUTPUT 14
int timeforoutput_index(tGrid *grid, int index)
{
  static int firstcall = 1;
  static int di[NOUTPUT];
  static double dt[NOUTPUT];
  static char output[NOUTPUT][20];
  char *name, s[20];
  int d, n;

  if (firstcall) {
    char *dstring[4] = {"D", "X", "Y", "Z"};
    char *ddstring[6] = {"XY", "XZ", "YZ", "XD", "YD", "ZD"};
    char *rstring[4] = {"max", "min", "norm","norminf"};

    firstcall = 0;

    /* cache */
    for (d = 0; d <= 3; d++) {
      sprintf(s, "%ddoutiter", d);
      di[d] = Geti(s);
      sprintf(s, "%ddouttime", d);
      dt[d] = Getd(s);
      sprintf(output[d], "%ddoutput", d);
    }
    for (d = 0; d <= 3; d++) {
      sprintf(output[d+4], "0doutput%s", rstring[d]);
      di[d+4] = di[0];
      dt[d+4] = dt[0];
    }
    for (d = 0; d <= 3; d++) {
      sprintf(output[d+4], "1doutput%s", dstring[d]);
      di[d+4] = di[1];
      dt[d+4] = dt[1];
    }
    for (d = 0; d < 6; d++) {
      sprintf(output[d+8], "2doutput%s", ddstring[d]);
      di[d+8] = di[2];
      dt[d+8] = dt[2];
    }
  }

  name = VarName(index);

  for (n = 0; n < NOUTPUT; n++)
    if (Getv(output[n], name)) 
      if (timeforoutput_di_dt(grid, di[n], dt[n])) 
	return 1;

  return 0;
}




/* decide whether it is time for output for any member of a variable list */
int timeforoutput(tGrid *grid, tVarList *vl)
{
  int i;

  for (i = 0; i < vl->n; i++)
    if (timeforoutput_index(grid, vl->index[i])) return 1;
  return 0;
}




/* master function to do all the writing */
int write_grid(tGrid *grid)
{
  static int firstcall = 1;
  static int di[4];
  static double dt[4];
  static char *ou[4], *ou0d[4], *ou1d[4], *ou2d[6];
  static int all[4];
  char s[20];
  int b, d;
  char *dstring[4]  = {"D", "X", "Y", "Z"};
  char *ddstring[6] = {"XY", "XZ", "YZ", "XD", "YD", "ZD"};
  char *rstring[4]  = {"max", "min", "norm","norminf"};

  if (firstcall) {
    firstcall = 0;

    /* cache */
    for (d = 0; d <= 3; d++) {
      sprintf(s, "%ddoutiter", d);
      di[d] = Geti(s);
      sprintf(s, "%ddouttime", d);
      dt[d] = Getd(s);
      sprintf(s, "%ddoutput", d);
      ou[d] = Gets(s);
      sprintf(s, "%ddoutputall", d);
      all[d] = Getv(s, "yes");
      sprintf(s, "1doutput%s", dstring[d]);
      ou1d[d] = Gets(s);
    }
    for (d = 0; d < 4; d++) {
      sprintf(s, "0doutput%s", rstring[d]);
      ou0d[d] = Gets(s);

    }
    for (d = 0; d < 6; d++) {
      sprintf(s, "2doutput%s", ddstring[d]);
      ou2d[d] = Gets(s);
    }
    if (Getv("outputall", "yes"))
      for (d = 0; d <= 3; d++) all[d] = 1;
  }

  /* 0d output */
  d = 0;
  if (timeforoutput_di_dt(grid, di[d], dt[d]))
  {
    printf("0d output not implemented");
  }

  /* 1d output */
  d = 1;
  if (timeforoutput_di_dt(grid, di[d], dt[d]))
  {
    for (b = 0; b < grid->nboxes; b++)
    {
      tBox *box = grid->box[b];
      char str[1000];
      int start=0;

      //printf("1dout ... |%s|\n", ou[1]);
      while(sscanf(ou[1]+start, "%s", str)==1)
      {
        start += strlen(str);
        if(ou[1][start]==' ') start++;
        //printf("1dout |%s|\n",str);
        xgraph_out1_boxvar(box, str);
      }
    }
  }

  /* 2d output */
  d = 2;
  if (timeforoutput_di_dt(grid, di[d], dt[d]))
  {
    for (b = 0; b < grid->nboxes; b++)
    {
      tBox *box = grid->box[b];
      char str[1000];
      int start=0;

      //printf("2dout ... |%s|\n", ou[2]);
      while(sscanf(ou[2]+start, "%s", str)==1)
      {
        start += strlen(str);
        if(ou[2][start]==' ') start++;
        //printf("2dout |%s|\n", str);
        gnuplot_out2d_boxvar(box, str);
      }
    }
  }

  /* 3d output */
  d = 3;
  if(timeforoutput_di_dt(grid, di[d], dt[d]))
  {
    printf("3d output not implemented");
  }
  return 0;
}

