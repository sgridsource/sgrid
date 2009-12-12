/* checkpoint_write.c */
/* Wolfgang Tichy 11/2008 */

/* write checkpoint files 
*/

#include "sgrid.h"
#include "checkpoint.h"




/* write checkpoint file 
   - every processor writes its own file
   - assumes that every processor can write its own file
   - assumes that there is one global file system and that it is better
     to have the processors write one after another rather than hitting
     the disks with simultaneous requests from each processor
     (make parameter to choose this behaviour)
*/
void checkpoint_write(tGrid *g)
{
  char *filename = checkpoint_filename("", "w");
  char *filename_new;
  FILE *fp;
  int i;

  /* name of backup file */
  filename_new = cmalloc(strlen(filename) + 20);
  sprintf(filename_new, "%s_new", filename);

  /* for all processes */
  for(i = 0; i < sgrid_mpi_size(); i++)
  {
    /* if it is my turn */
    if(i == sgrid_mpi_rank())
    {
      /* open file for writing */
      fp = fopen(filename_new, "wb");
      if (!fp) errorexits("failed opening %s", filename_new);

      /* now we are ready to go! */
      checkpoint_write_local(g, fp);

      /* done */
      fclose(fp);
    }
    /* everyone please wait here */
    sgrid_mpi_barrier();
  }

  /* Move filename_new to filename, but only after we are sure
     all has been written.                                       */
  fp = fopen(filename_new, "rb");
  if(fp)
  {
    fclose(fp);
    system3("mv -f", filename_new, filename); 
  }
  
  /* clean up */
  free(filename_new);
  free(filename);
}




/* write checkpoint file on one processor 
   - should add checks for grid size etc. */
void checkpoint_write_local(tGrid *g, FILE *fp)
{
  tBox *box = g->box[0];
  int i, j, b;
  int n, nall = 0, nv;
  int pr = 0;
  int number_of_variables;

  /* write one line text header, check with "head -1 checkpoint.0" */
  fprintf(fp, "$BEGIN_checkpoint_header:  ");
  fprintf(fp, "$bmin = %d , $bmax = %d , ", 0, g->nboxes-1);
  fprintf(fp, "$time = %g , ", g->time);
  fprintf(fp, "$iteration = %d\n", g->iteration);
  fprintf(fp, "$END_checkpoint_header\n\n");

  /* write parameter database
     some parameters are used like global variables and change during evolution
  */
  fprintf(fp, "$BEGIN_parameter_database:  ");
  fprintf(fp, "$number_of_parameters = %d\n", GetnParameters());
  for (i = 0; i < GetnParameters(); i++)
  {
    fprintf(fp, "%s = ", GetnameInd(i));
    fprintf(fp, "%s\n", GetsInd(i));
  }
  fprintf(fp, "$END_parameter_database\n\n");


  fprintf(fp, "$BEGIN_variables:\n");
  /* for all boxes */
  for (b = 0; b <= g->nboxes-1; b++)
  {
    box = g->box[b];
    n = nv = 0;

    /* write all variables */
    if (Getv("checkpoint_variables", "all"))
    {
      /* for all variables with storage */
      for (i = 0; i < g->nvariables; i++)
      {
	if (box->v[i] != 0)
	{
	  /* write */
	  fprintf(fp, "$variable = %s : box = %d , box->nnodes = %d\n",
	          VarName(i), b, box->nnodes);
	  fwrite(box->v[i], sizeof(double), box->nnodes, fp);
          fprintf(fp, "\n");

	  n += box->nnodes;
	  nv++;
	  if (pr) printf("wrote box %d, nnodes=%d, variable=%s\n", 
			 b, box->nnodes, VarName(i));
	}
      }
    }

    /* write only those variables that are needed
       - needs to be tested, made more flexible!
    */
    else
    {
      tVarList *vl = checkpoint_varlist(g);
      
      /* for all variables in list with storage */
      for (j = 0; j < vl->n; j++)
      {
	i = vl->index[j];
	if (box->v[i] != 0)
	{
	  /* write */
	  fprintf(fp, "$variable = %s : box = %d , box->nnodes = %d\n",
	          VarName(i), b, box->nnodes);
	  fwrite(box->v[i], sizeof(double), box->nnodes, fp);
          fprintf(fp, "\n");

	  n += box->nnodes;
	  nv++;
          if (pr) printf("wrote box %d, nnodes=%d, variable=%s\n",
			 b, box->nnodes, VarName(i));
	}
      }
      vlfree(vl);
    }
    /* info */
    if (1) printf("wrote b %d, %d variables, total of %d points\n", b, nv, n);
    nall += n;
  }
  fprintf(fp, "$END_variables\n\n");

  /* info */
  if (1) printf("wrote boxes %d to %d, %d variables, %d doubles, %d MB\n",
		0, g->nboxes-1, nv, nall,
		(int)(nall * sizeof(double) / 1000000));
}
