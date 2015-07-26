/* utilities.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99 */

#include <time.h>

#include "sgrid.h"
#include "main.h"

/* for POSIX.1-2001 mkdir, opendir, fork, wait functions */
#include <unistd.h>     /* for fork */
#include <sys/stat.h>
#include <sys/types.h>  /* for pid_t */
#include <sys/wait.h>   /* for wait */
#include <dirent.h>     /* for opendir */


/* global vars for timing */
double time_in_s_at_sgrid_start; /* set in main.c */



/* debug */
void yo(void)     {fprintf(stdout, "yo!\n");fflush(stdout);}
void Yo(double x) {fprintf(stdout, "Yo:%g\n", x);fflush(stdout);}



/* pretty */
void prdivider(int n)
{
  printf("------------------------------------------------------------------------------\n");
  fflush(stdout);
}


/* get current time in seconds */
double getTimeIn_s(void)
{
  double t_in_s = time(NULL) - time_in_s_at_sgrid_start;
  return t_in_s;
}

/* print current time */
void prTimeIn_s(char *comment)
{
  double t_in_s = time(NULL) - time_in_s_at_sgrid_start;
  printf("%s%.2fs\n", comment, t_in_s);
  fflush(stdout);
}

/* get current CPU clock time in seconds */
double getClockTimeIn_s(void)
{
  double t_in_s = clock();
  t_in_s = t_in_s/CLOCKS_PER_SEC;
  return t_in_s;
}

/* print current CPU clock time */
void prClockTimeIn_s(char *comment)
{
  double t_in_s = clock();
  t_in_s = t_in_s/CLOCKS_PER_SEC;
  printf("%s%gs\n", comment, t_in_s);
  fflush(stdout);
}


/* minimum and maximum, also works for integers in most places */
double min2(double x, double y)
{
  return (x < y) ? x : y;
}
double min3(double x, double y, double z)
{
  return min2(min2(x, y), z);
}
double max2(double x, double y)
{
  return (x > y) ? x : y;
}
double max3(double x, double y, double z)
{
  return max2(max2(x, y), z);
}


/* find min in a 1d array f, return min and set imin to index of min */
double min_in_1d_array(double *f, int n, int *imin)
{
  int i;   
  double fmin=f[0];

  *imin=0;
  for(i=1; i<n; i++)
    if(f[i]<fmin) 
    {
      fmin = f[i];
      *imin= i;
    }
  return fmin;
}
/* find max in a 1d array f, return max and set imax to index of max */
double max_in_1d_array(double *f, int n, int *imax)
{
  int i;   
  double fmax=f[0];

  *imax=0;
  for(i=1; i<n; i++)
    if(f[i]>fmax) 
    {
      fmax = f[i];
      *imax= i;
    }
  return fmax;
}
/* find min in two 1d arrays f0, f1, 
   return min, set ai to 0 or 1 depending on which array the min is in,
   set imin to index of min in array ai */
double min2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imin)
{
  int imin0, imin1;   
  double fmin0, fmin1, fmin;
  
  fmin0 = min_in_1d_array(f0, n0, &imin0);
  fmin1 = min_in_1d_array(f1, n1, &imin1);

  if(fmin1<fmin0) { *ai=1;  fmin=fmin1;  *imin=imin1; }
  else            { *ai=0;  fmin=fmin0;  *imin=imin0; }  
  return fmin;
}
/* find max in two 1d arrays f0, f1, 
   return max, set ai to 0 or 1 depending on which array the max is in,
   set imax to index of max in array ai */
double max2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imax)
{
  int imax0, imax1;   
  double fmax0, fmax1, fmax;
  
  fmax0 = max_in_1d_array(f0, n0, &imax0);
  fmax1 = max_in_1d_array(f1, n1, &imax1);

  if(fmax1>fmax0) { *ai=1;  fmax=fmax1;  *imax=imax1; }
  else            { *ai=0;  fmax=fmax0;  *imax=imax0; }  
  return fmax;
}
/* find min in three 1d arrays f0, f1, f2, 
   return min, set ai to 0,1,2 depending on which array the min is in,
   set imin to index of min in array ai */
double min3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imin)
{
  double fmin;

  min2_in_1d_array(f1,n1, f2,n2, ai, imin);
  if(*ai==1) { fmin = min2_in_1d_array(f0,n0, f2,n2, ai, imin); *ai = (*ai)*2; }
  else       { fmin = min2_in_1d_array(f0,n0, f1,n1, ai, imin); }
  return fmin;
}
/* find max in three 1d arrays f0, f1, f2, 
   return max, set ai to 0,1,2 depending on which array the max is in,
   set imax to index of max in array ai */
double max3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imax)
{
  double fmax;

  max2_in_1d_array(f1,n1, f2,n2, ai, imax);
  if(*ai==1) { fmax = max2_in_1d_array(f0,n0, f2,n2, ai, imax); *ai = (*ai)*2; }
  else       { fmax = max2_in_1d_array(f0,n0, f1,n1, ai, imax); }
  return fmax;
}


/* make copy of a file: cp fname newname */
int copy_file(char *fname, char *newname)
{
  FILE *in, *out;
  void *buffer;
  size_t BUFSIZE=16777216; /* 16MiB */
  size_t bcount;

  printf("copy_file(\"%s\", \"%s\");\n", fname, newname);

  /* open source file */
  in = fopen(fname, "rb");
  if(!in) errorexits("failed opening %s", fname);

  /* open destination file */
  out = fopen(newname, "wb");
  if(!out) errorexits("failed opening %s", newname);

  /* copy char by char */
  /* while( (ch=fgetc(in)) != EOF)
       fputc(ch, out);             */
  buffer = calloc(BUFSIZE , sizeof(char));
  if(!buffer)
    errorexiti("copy_file: out of memory for buffer (%d chars)", BUFSIZE);
  do
  {
    bcount = fread(buffer, sizeof(char), BUFSIZE, in);
    fwrite(buffer, sizeof(char), bcount, out);
  } while(bcount==BUFSIZE);
  free(buffer);

  fclose(out);
  fclose(in);
  return 0;
}

/* make a copy of a file in some dir */
int copy_file_into_dir(char *fname, char *dir)
{
  char newname[10000];
  int i;

  /* find / in fname, to determine filename without dirname */
  for(i=strlen(fname)-1; i>=0; i--) if(fname[i]=='/') break;
  snprintf(newname, 9999, "%s/%s", dir, fname+i+1);

  return copy_file(fname, newname);
}

/* use opendir to scan through dir and remove the entire dir */
int remove_dir(char *which_dir)
{
  DIR           *d;
  struct dirent *dir;
  char file[1000];
  struct stat s;
  
  d = opendir(which_dir);
  if(d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      /* exclude . and .. directories */
      if( strcmp( dir->d_name, "." ) == 0 || strcmp( dir->d_name, ".." ) == 0)
        continue;

      snprintf(file, 999, "./%s/%s", which_dir, dir->d_name);
      //printf("*"); //print * for every deleted file

      if(opendir(file)!=NULL) remove_dir(file);
      else
      {
        if(remove(file) != 0)
        {
          //printf("\n%s\n", file);
          //perror("Remove failed");
          closedir(d);
          return -2;
        }
      }
    } /* end of while loop */

    closedir(d);

    /* delete directory */
    if(remove(which_dir) != 0)
    {
      //printf("%s\n", which_dir);
      //perror("Remove failed");
      return -1;
    }
  }
  return 0;
}

/* ugh, but how universal are those built in functions? */
int system2(char *s1, char *s2) 
{
  return system3(s1, s2, "");
}
int system3(char *s1, char *s2, char *s3) 
{
  char command[10000];
  int status = 0;

  /* check for special cases where we can use c-functions */
  if( strcmp(s1,"mv")==0 || strcmp(s1,"mv -f")==0 ) /* use rename */
  {
    sprintf(command, "rename(\"%s\", \"%s\");", s2, s3);
    status = rename(s2, s3);
    printf("ANSI C call: %s\n", command);
  }
  else if( strcmp(s1,"rm -rf")==0 ) /* use remove */
  {
    if(strlen(s2)>0)
    {
      status = remove_dir(s2);
      printf("remove_dir(\"%s\");\n", s2);
    }
    if(strlen(s3)>0)
    {
      status = remove_dir(s3);
      printf("remove_dir(\"%s\");\n", s3);
    }
  }
  else if( strcmp(s1,"mkdir")== 0 ) /* use POSIX.1-2001 mkdir function */
  {
    sprintf(command, "mkdir(\"%s\", S_IRWXU | S_IRWXG);", s2);
    status = mkdir(s2, S_IRWXU | S_IRWXG);
    printf("POSIX.1-2001 call: %s\n", command);
  }
  else /* use system */
  { 
    sprintf(command, "%s %s %s", s1, s2, s3);
    status = system(command);
    printf("System call: %s\n", command);
  }
  
  if(status!=0) printf(" -> WARNING: Return value = %d\n", status);
  return status;
}

/* construct an argv array from a string and return number of args */
int construct_argv(char *str, char ***argv)
{
  char *str1, *token, *saveptr;
  int count;

  *argv = NULL;
  for(count=0, str1=str; ; count++, str1=NULL)
  {
    *argv = (char **) realloc(*argv, sizeof(char *)*(count+1));
    token = strtok_r(str1, " ", &saveptr);
    //printf("token=%p:%s\n", token,token);
    (*argv)[count] = token;
    if(token == NULL) break;
  }
  //printf("saveptr=%p:%s\n", saveptr,saveptr);
  return count;
}

/* run a command, without a shell */
int system_emu(const char *command)
{
  char *com = strdup(command); /* duplicate since construct_argv modifies its args */
  int ret, status;
  printf("system_emu: running command:\n%s\n", command);

  /* Spawn a child to run the program. */
  pid_t cpid = fork();
  if(cpid==0) /* child process */
  {
    char **argv;
    construct_argv(com, &argv);
    ret = execv(argv[0], argv);
    printf("*** WARNING: command not found, (execv returned %d) ***\n", ret);
    exit(127); /* exit child, only if execv fails */
  }
  else /* cpid!=0; parent process */
  {
    waitpid(cpid, &ret, 0); /* wait for child to exit */
    status = WEXITSTATUS(ret);
  }
  if(status!=0) printf(" -> WARNING: Return value = %d\n", status);
  free(com);
  return status;
}



/* preliminary ... */
double *dmalloc(int n)
{
  double *p = (double *) malloc(sizeof(double) * n);
  
  if (!p) errorexiti("out of memory (%d double)", n);
  return p;
}

int *imalloc(int n)
{
  int *p = (int *) malloc(sizeof(int) * n);
  
  if (!p) errorexiti("out of memory (%d int)", n);
  return p;
}

char *cmalloc(int n)
{
  char *p = (char *) malloc(sizeof(char) * n);
  
  if (!p) errorexiti("out of memory (%d char)", n);
  return p;
}

void *pmalloc(int n)
{
  void *p = malloc(sizeof(void *) * n);
  
  if (!p) errorexiti("out of memory (%d void *)", n);
  return p;
}





/* the one function every program should have */
/* note that sgrid_main.h defines a macro so that the user does not have
   to specify __FILE__ and __LINE__ for location where the error occured
*/
#undef errorexit
#undef errorexits
#undef errorexiti

void errorexit(char *file, int line, char *s) 
{
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  //sgridpi_finalize();
  exit(1);
}

void errorexits(char *file, int line, char *s, char *t)
{
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  //sgridpi_finalize();
  exit(1);
}

void errorexiti(char *file, int line, char *s, int i)
{
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, i);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  //sgridpi_finalize();
  exit(1);
}

/* do not write functions beyond this line because the undef/define 
   method for the errorexit functions means that they should go last */
