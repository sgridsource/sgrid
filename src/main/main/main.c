/* main.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "main.h"


/* global var that tells us how often we have restarted sgrid */
int sgrid_restarts;


/* initialize libraries 
   the automatically generated file calls the initializers for each module */
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
  int restart, lnum;
  int largc; char **largv;

  /* init MPI if it is there, otherwise do nothing */
  sgrid_MPI_Init(&argc, &argv);
  if(sgrid_MPI_size()>1) printf("MPI_size=%d: process%d started.\n", 
                                sgrid_MPI_size(), sgrid_MPI_rank());

  /* sgrid can restart itself, if we start it like this:
     sgrid --argsfile args.txt
     The file args.txt contains lines with sgrid args (parfile and options).
     sgrid will go through these lines and run with each line. If there
     is more than one line it will restart. If sgrid is compiled with
     DFLAGS += -DUSEMPI
     and we have several MPI processes, each process will do a different
     line. Lines will be equally divided among MPI processes. E.g. proc0
     will do every third line if there are 3 MPI proceesses. */
  sgrid_restarts=0;
  do
  {
    char buffer[5000];
    restart = get_next_argument_line(argc, argv, buffer, &largc, &largv, &lnum);

    /* do nothing at all if rank>0 and there is no line */
    if( sgrid_MPI_rank()>0 && restart==0 ) break;
    /* continue with next line if this line is not for this rank */
    if( (lnum+sgrid_MPI_rank()) % sgrid_MPI_size() != 0 ) continue;

    if(sgrid_restarts>0) /* check if this a restart, i.e. not first run */
    {
      if(restart==0) break; /* stop if no more command line args */
      prdivider(0);         /* otherwise print info */
      printf("*** sgrid restart %d ***\n", sgrid_restarts);
      printf("This is fragile: All static C-vars are unchanged on restart! Vars and Funs\n"
             "are also unchanged! This is wrong if certain Pars change on restart!\n");
    }
    
    /* start the main parts of sgrid: */
    initTimeIn_s();
    read_command_line(largc, largv);
    parse_parameter_file(Gets("parameterfile"));
    parse_command_line_options();
    make_output_directory();
    initialize_libraries();

    iterate_parameters(0); /* start of new iteration */
    while(iterate_parameters(1))
    {
      RunFun(POST_PARAMETERS, 0); /* hook for funs right after iterate_parameters */
      RunFun(PRE_GRID, 0);  /* hook for special grid preparation */
      g = make_grid(1);
      RunFun(POST_GRID, g); /* hook for special treatment after grid creation */
      initialize_grid(g);
      evolve_grid(g);
      finalize_grid(g);
      RunFun(POST_FINALIZE_GRID, g); /* hook after finalize_grid, e.g. for special cleanup */
      makeparameter("outdir_previous_iteration", "", "outdir of previous iteration");
      Sets("outdir_previous_iteration", Gets("outdir"));
    }
    free_global_parameter_database_contents(); /* free strings in parameter database pdb */
    /* NOTE: currently we do not free vdb in variables.c */
    /*       we also do not free fps in skeleton.c */
    sgrid_restarts++; /* sgrid_restarts counts the actual restarts */
  } while(restart); /* end restarts loop */

  sgrid_MPI_Finalize();
  return 0;
}



/* Get multiple command lines from a file and return the args in 
   largc, largv, also return the number of each line in lnum.
   The function's return value is 1 if we got an argsline from file. */
int get_next_argument_line(int argc, char **argv, char *buffer,
                           int *largc, char ***largv, int *lnum)
{
  static FILE *fp=NULL;
  static arglinenum=0; /* contains line number of args we read */
  char *argsfile;
  char line[5000];
  char *str = buffer;
  int ret;

  /* default is to just return argc,argv in largc,largv */
  *largc = argc;
  *largv = argv;
  *lnum = arglinenum;

  /* do nothing if we are not in multiple command line mode */
  if(argc!=3) return 0;
  if(strcmp(argv[1], "--argsfile")!=0) return 0;

  argsfile = argv[2];
  printf("Reading multiple command line arguments from file %s\n", argsfile);

  /* open argsfile  */
  if(fp==NULL) fp = fopen(argsfile, "r");
  if(fp==NULL)
  {
    printf("Could not open file \"%s\"\n", argsfile);
    errorexit("");
  }

  /* get one line from argsfile */
  do
  {
    ret = fscanline(fp, line);
    if(ret==EOF) { fclose(fp); fp=NULL; arglinenum=0; return 0;}
    if(ret>0) if(line[0]=='#') ret=0; /* skip lines starting with # */
  } while(ret==0);

  /* str = argv[0] + line */
  strcpy(str, argv[0]);
  strcat(str, " ");
  strcat(str, line);

  /* set largc,largv pointers from str */
  *largc = construct_argv(str, largv);
  arglinenum++; /* inc arglinenum so that we know how many we got */

  return 1;
}


/* read command line */
int read_command_line(int argc, char **argv)
{
  int i; 

  if (0) 
    for (i = 0; i < argc; i++)
      printf("argv[%d] = %s\n", i, argv[i]);

  if (argc < 2)
  {
    printf("Welcome to sgrid.\n");
    printf("Usage:  sgrid name.par\n");
    printf("or:     sgrid name.par options and extra arguments\n");
    printf("or:     sgrid --argsfile args.txt\n");
    printf("\n");
    printf("options: --keep_previous           do not touch name_previous\n");
    printf("         --modify-par:\"P=v\"        set par P to value v\n");
    printf("         --argsfile                read sgrid args from file\n");
    printf(" all options must start with --\n"); 
    exit(0);
  }

  prdivider(0);
  printf("Welcome to sgrid.\n");
  prdivider(0);

  /* got two parameters */
  if (argc >= 2)
  {
    int nopts, nargs;
    char argi[1000];
    char descr[1000];
    char options[1000];
    char *parfile = (char *) calloc(strlen(argv[1])+40, sizeof(char));
    char *outdir  = (char *) calloc(strlen(argv[1])+40, sizeof(char));
    int parnamelen;

    /* determine name of parameter file and output directory */
    strcpy(parfile, argv[1]);
    parnamelen = strlen(parfile);
    if(!strstr(parfile, ".par") || parnamelen<5) strcat(parfile, ".par");
    parnamelen = strlen(parfile);
    if(!strstr(parfile+parnamelen-4, ".par")) strcat(parfile, ".par");
    strcpy(outdir, parfile);
    outdir[strlen(outdir)-4]=0; /* remove .par */

    /* first parameter initializes parameter data base */
    printf("Adding command line parameters\n");
    AddPar("outdir", outdir, "output directory");
    AddPar("parameterfile", parfile, 
	   "name of parameter file given on command line");
    // AddPar("trace_memory", "no", "enable memory tracing");

    /* add other args */
    nargs=nopts=0;
    options[0]=0;
    for(i=2; i<argc; i++)
    {
      if(argv[i][0]=='-')
      {
        if(strlen(argv[i])==1) errorexit("- is not a valid option");
        if(argv[i][1]!='-') errorexit("all options must start with --");
        /* save all options in options */
        if(nopts>0) strncat(options, " ", 999);
        strncat(options, argv[i], 999);
        nopts++;
      }
      else
      {
        snprintf(argi, 999, "sgrid_arg%d", nargs+2);
        snprintf(descr, 999, "sgrid command line argument%d", nargs+2);
        AddPar(argi, argv[i], descr);
        nargs++;
      }
    }
    /* add sgrid command line options */
    if(nopts>0) AddPar("sgrid_options", options, "sgrid command line options");

    free(parfile);
    free(outdir);
  }

  /* more initialization */

  return 0;
}

int make_output_directory(void)
{
  char *outdir  = Gets("outdir");
  char *outdirp = (char *) calloc(strlen(outdir)+40, sizeof(char));

  /* set outdirp to outdir_previous */
  strcpy(outdirp, outdir);
  strcat(outdirp, "_previous");

  /* check if a shell is available to execute commands later */
  /* NOTE: system2 and system3 are smart enough to do "mkdir", "rm -rf"
           and "mv" even without a shell by using POSIX calls. */
  if(system(NULL)==0)
  {
    printf("WARNING: system(NULL)=0 => cannot execute shell commands!\n");
    printf("         Consider using system_emu.\n");
  }

  /* check if we remove outdir_previous */
  if(!GetvLax("sgrid_options", "--keep_previous"))
  {
    /* remove outdir_previous and move outdir to outdir_previous */
    system2("rm -rf", outdirp);
    system3("mv", outdir, outdirp);
  }

  /* make output directory, save parfile */
  system2("mkdir", outdir);
  /* system3("cp", Gets("parameterfile"), outdir); */
  copy_file_into_dir(Gets("parameterfile"), outdir);

  /* redirect stdout and stderr for MPI jobs 
     all output is collected in outdir/stdout.001 etc  */
  if (sgrid_MPI_rank()>0)
  {
    char s[1000], f[100];
    sprintf(f, "%%s/stdout.%%0%dd", (int) log10(sgrid_MPI_size())+1);
    sprintf(s, f, outdir, sgrid_MPI_rank());  
    freopen(s, "w", stdout);   
    freopen(s, "w", stderr);
  }

  free(outdirp);
  return 0;
}


/* go through options and act accordingly */
int parse_command_line_options(void)
{
  char *optionstr;
  char *str1;
  char *str2;
  char *par;
  char *val;
  
  /* get length of sgrid_options string */
  if(GetsLax("sgrid_options")==0) return 0;
  printf("Parsing command line options\n");

  /* parse for all --modify-par: */
  optionstr = strdup((GetsLax("sgrid_options")));
  str1 = optionstr;
  while( (str1=strstr(str1, "--modify-par:"))!=NULL )
  {
    str1+=13;
    str2=strstr(str1, "=");
    if(str2==NULL) break;
    str2[0]=0; /* replace = with 0 */
    par=str1;
    val=str2+1;
    str2=strstr(val, " --");
    if(str2!=NULL) str2[0]=0; /* now val is 0 terminated for sure */
    if(0) printf("par=%s|val=%s|\n", par,val);
    
    /* set par to new value */
    AddOrModifyPar(par, val, "set with --modify-par option");

    /* move forward in str1 */
    if(str2!=NULL) str2[0]=' '; /* restore space before -- */
    str1=val;
    if(0) printf("str1=%s|\n", str1);
  }
  free(optionstr);
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
