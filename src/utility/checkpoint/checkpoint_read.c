/* checkpoint_read.c */
/* Wolfgang Tichy 11/2008 */

/* read checkpoint files 
*/

#include "sgrid.h"
#include "checkpoint.h"



/* read checkpoint file 
   - every processor reads its own file
*/
void checkpoint_read(tGrid *g)
{
  FILE *fp;
  int i, l;
  char str[10000];

  /* if there is not a complete set of checkpoint files,
     don't do anything: assume that this run starts from the beginnning
  */
  fp = checkpoint_openfiles("_previous");
  if (!fp) return;

  /* for all processes: read the pars and iteration numbers */
  for (i = 0; i < sgrid_mpi_size(); i++)
  {
    /* read if it is my turn */
    if (i == sgrid_mpi_rank())
      checkpoint_read_ParsAndIterations_local(g, fp);
    /* everyone please wait here */
    sgrid_mpi_barrier();
  }

  /* move all boxes according to what the new pars say */
  /*
  printf("moving boxes where they should be according to the parameters:\n");
  for(l = g->lmax; l >= g->lmin; l--)
    move_box(g, l);
  */

  /* for all processes: read the vars */
  for (i = 0; i < sgrid_mpi_size(); i++)
  {
    /* read if it is my turn */
    if (i == sgrid_mpi_rank())
      checkpoint_read_Vars_local(g, fp);

    /* everyone please wait here */
    sgrid_mpi_barrier();
  }

  /* clean up */
  fclose(fp);
}


/* utility function: strip trailing newline */
void striptrailingnewline(char *s)
{
  int n = strlen(s) - 1;
  if (n >= 0 && s[n] == '\n')
    s[n] = 0;
}


/* read pars, iteration numbers and times from checkpoint file on 
   one processor */
void checkpoint_read_ParsAndIterations_local(tGrid *grid, FILE *fp)
{
  char name[10000], value[10000], str[10000], *currentvalue;
  int i,j,ret;
  int bmin,bmax;
  int numberofpars, parsread;

  if (1) printf("reading parameters and iterations from checkpoint file:\n");

  /* read one line text header, check with "head -1 checkpoint.0" */
  if(fgotonext(fp, "$BEGIN_checkpoint_header:") == EOF)
    errorexit("$BEGIN_checkpoint_header: is missing!");
  fgetparameter(fp, "$bmin", value);  bmin = atoi(value);
  fgetparameter(fp, "$bmax", value);  bmax = atoi(value);
  fgetparameter(fp, "$time", value);      grid->time = atof(value);
  fgetparameter(fp, "$iteration", value); grid->iteration = atoi(value);
  printf("found $bmin = %d , $bmax = %d , $time = %g, $iteration = %d\n",
         bmin, bmax, grid->time, grid->iteration);

  /* read parameter database
     some parameters are used like global variables and change during evolution
  */
  if(fgotonext(fp, "$BEGIN_parameter_database:") == EOF)
    errorexit("$BEGIN_parameter_database: is missing!");

  fgetparameter(fp, "$number_of_parameters", value);
  numberofpars = atoi(value);
  printf("  found %d parameters\n", numberofpars);

  parsread = 0;
  while( (ret=fscanline(fp, str)) != EOF )
  {
    while(ret==0)  ret = fscanline(fp, str);
    if( ret==EOF || strcmp(str,"$END_parameter_database")==0 ) break;

    if(0) printf("  str: %s\n", str);
    extrstr_before_after_EQ(str, name, value);
 
    /* remove spaces at the end of name and the beginning of value */
    strcpy(str, name);
    sscanf(str, "%s", name);
    strcpy(str, value);
    j = strspn(str, " ");  /* strcpy(value, str+j); */
    if(j>0) strcpy(value, str+1); /* <--we remove at most one space */
    if(0) printf("  found parameter %s = ", name);
    if(0) printf("%s\n", value);
    currentvalue = GetsLax(name);
    if(currentvalue==0)
    {
      printf("skipping %s = %s\n", name, value);
      continue;
    }
    if(strcmp(value, currentvalue) != 0)
    {
      /* handle list of parameters that are allowed to keep a new value
         (should be made actual list)
      */
      if(strcmp(name, "checkpoint") &&
	 strcmp(name, "finaltime")  &&
	 strcmp(name, "iterations")  )
      {
	if (1) printf("overwriting %s = %s -> %s\n",
		      name, currentvalue, value);
	Sets(name, value);
      }
    }
    parsread++;
  }
  printf("  read %d parameters\n", parsread);
}



/* read vars from checkpoint file on one processor */
void checkpoint_read_Vars_local(tGrid *grid, FILE *fp)
{
  tBox *box;
  int ilax, i, j, l, nnodes_on_l;
  int lold = -11111;
  int n, nall = 0, nv;
  int pr = 0;
  char varname[10000], value[10000];

  if (1) printf("reading variables from checkpoint file:\n");
  if(fgotonext(fp, "$BEGIN_variables:") == EOF)
    errorexit("$BEGIN_variables: is missing!");  
  
  l = lold;
  /* for all saved variables */
  while(fgetparameter(fp, "$variable", varname) != EOF)
  {
    ilax = IndLax(varname);
    if(ilax>=0)
      i = ilax;
    else
    {
      printf("variable %s in checkpoint file does not exist!\n", varname);
      i = 0;
    }

    fgetparameter(fp, "box", value);
    l = atoi(value);
    box = grid->box[l];
    if(l != lold)
    {
      lold = l;
      n = nv = 0;
    }

    fgetparameter(fp, "box->nnodes", value);
    nnodes_on_l  = atoi(value);
    if(nnodes_on_l != box->nnodes)
      errorexit("box->nnodes in checkpoint file is not what it should be!");

    /* read the trailing newline behind the box->nnodes par */
    fscanline(fp, value);

    /* read that variable */
    if(ilax>=0 && box->v[i] != 0)
    {
      /* read */
      fread(box->v[i], sizeof(double), nnodes_on_l, fp);

      nv++;
      n += nnodes_on_l;
      nall += nnodes_on_l;
      if(pr) printf("read box %d, nnodes=%d, variable=%s, index=%d,\n", 
		     l, nnodes_on_l, VarName(i), i);
    }
    else
    {
      //errorexit("we should overread the next nnodes_on_l doubles in fp");
      printf("skipping the %d double numbers of %s (in box %d) ...\n",
             nnodes_on_l, varname, l);
      for(j=1; j<=nnodes_on_l*sizeof(double)/sizeof(char); j++) fgetc(fp);
    }
  }
  /* info */
  if (1) printf("read boxes %d to %d, %d variables, %d doubles, %d MB\n",
		0, grid->nboxes-1, nv, nall, 
		(int)(nall * sizeof(double) / 1000000));
}
