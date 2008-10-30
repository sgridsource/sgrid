/* main.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "main.h"



/* initialize libraries 
   the automatically generated file calls the initializers for each module
*/
void initialize_libraries(void)
{
  prdivider(0);
  printf("Initializing libraries\n");

#include "sgrid_automatic_initialize.c"
}




/**************************************************************************/
/* main */
int main(int argc, char **argv) 
{
  tGrid *g;

  read_command_line(argc, argv);
  parse_parameter_file(Gets("parameterfile"));
  initialize_libraries();

  while (iterate_parameters()) {
    g = make_grid(1);
    initialize_grid(g);
    evolve_grid(g);
    finalize_grid(g);
  }
  return 0;
}




/* read command line */
int read_command_line(int argc, char **argv)
{
  int i; 

  if (0) 
    for (i = 0; i < argc; i++)
      printf("argv[%d] = %s\n", i, argv[i]);

  if (argc != 2)
  {
    printf("Welcome to sgrid.\n");
    printf("Usage:  sgrid name.par\n");
    exit(0);
  }

  prdivider(0);
  printf("Welcome to sgrid.\n");
  prdivider(0);

  /* got two parameters */
  if (argc == 2) {
    char *parfile = (char *) calloc(sizeof(char), strlen(argv[1])+40);
    char *outdir  = (char *) calloc(sizeof(char), strlen(argv[1])+40);
    char *outdirp = (char *) calloc(sizeof(char), strlen(argv[1])+40);

    /* determine name of parameter file and output directory */
    strcpy(parfile, argv[1]);
    if (!strstr(parfile, ".par")) strcat(parfile, ".par");
    strcpy(outdir, parfile);
    *strstr(outdir, ".par") = '\0';
    
    /* make output directory, save current one */
    strcpy(outdirp, outdir);
    strcat(outdirp, "_previous");
    system2("rm -rf", outdirp);
    system3("mv", outdir, outdirp); 
    system2("mkdir", outdir);
    system3("cp", parfile, outdir);

    /* first parameter initializes parameter data base */
    printf("Adding command line parameters\n");
    AddPar("outdir", outdir, "output directory");
    AddPar("parameterfile", parfile, 
	   "name of parameter file given on command line");
    // AddPar("trace_memory", "no", "enable memory tracing");
  }

  /* more initialization */
  return 0;
}




/* initialize grid */
int initialize_grid(tGrid *g)
{
  if (1) {
    prdivider(0);
    printf("Initializing grid\n");
  }

  /* some things need to be done before initial data */
  RunFun(PRE_INITIALDATA, g); 

  /* compute initial data */
  RunFun(INITIALDATA, g); 

  /* some things need to be done after initial data */
  RunFun(POST_INITIALDATA, g); 

  /* initial data is just another new time slice */
  RunFun(PRE_POST_EVOLVE, g);
  RunFun(POST_EVOLVE, g); 

  /* initial data complete */
  prdivider(0);
  printf("Done with initialization\n");
  printf(" iteration %d, time=%g\n", g->iteration, g->time);

  /* analyze initial data */
  RunFun(ANALYZE, g);

  /* output for permanent variables */
  RunFun(OUTPUT, g);

  /* post output */
  RunFun(POST_OUTPUT, g);

  /* checkpoint, just in case we need it here already */
  checkpoint(g);

  return 0;
}




/* evolve grid */
int evolve_grid(tGrid *grid)
{
  int iterationmax = Geti("iterations");
  double timemax = Getd("finaltime");

  prdivider(0);

  if (timemax > 0)
    iterationmax = timemax/grid->dt + 0.5;

  if (iterationmax > 0) 
    printf("Evolving grid for %d top grid iterations to time %.3f\n", 
	   iterationmax, iterationmax * grid->dt);

  if (iterationmax <= 0) return 0;

  /* outer most evolution loop */
  while (grid->iteration < iterationmax)
  { 
    /* pre evolve */
    RunFun(PRE_EVOLVE, grid); 

    /* evolve */
    RunFun(EVOLVE, grid); 

    /* post evolve */
    RunFun(PRE_POST_EVOLVE, grid);
    RunFun(POST_EVOLVE, grid); 

    /* evolution step complete */
    grid->iteration++;
    grid->time = grid->iteration * grid->dt;

    /* print info */
    printf(" iteration %d, time=%g\n", grid->iteration, grid->time);
    fflush(stdout); 

    /* analyze */
    RunFun(ANALYZE, grid);

    /* output for permanent variables */
    RunFun(OUTPUT, grid);

    /* post output */
    RunFun(POST_OUTPUT, grid);

    /* checkpoint */
    checkpoint(grid);
  
    /* update since this may change during evolution, say when checkpointing */
    timemax= Getd("finaltime");
    iterationmax= (timemax > 0) ? timemax/grid->dt + 0.5 : Geti("iterations");
  }
  return 0;
}




/* finalize grid */
int finalize_grid(tGrid *g)
{
  prdivider(0);
  free_grid(g);
  return 0;
}
