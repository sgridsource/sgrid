/* utilities.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99 */

#include <time.h>

#include "sgrid.h"
#include "main.h"




/* debug */
void yo(void)     {fprintf(stdout, "yo!\n");fflush(stdout);}
void Yo(double x) {fprintf(stdout, "Yo:%g\n", x);fflush(stdout);}



/* pretty */
void prdivider(int n)
{
  printf("------------------------------------------------------------------------------\n");
  fflush(stdout);
}


/* print current time */
void prTimeIn_s(char *comment)
{
  double t_in_s = clock()/CLOCKS_PER_SEC;
  printf("%s%g\n", comment, t_in_s);
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
    printf("ANSI C call:  %s\n", command);
  }
  else if( strcmp(s1,"rm -rf")==0 && 0 ) /* use remove */
  {
    if(strlen(s2)>0)
    {
      sprintf(command, "remove(\"%s\");", s2);
      status = remove(s2); /* Note: remove fails if dir is not empty */
      printf("ANSI C call:  %s\n", command);
    }
    if(strlen(s3)>0)
    {
      sprintf(command, "remove(\"%s\");", s3);
      status = remove(s3);
      printf("ANSI C call:  %s\n", command);
    }
  }
  else /* use system */
  { 
    sprintf(command, "%s %s %s", s1, s2, s3);
    status = system(command);
    printf("System call:  %s\n", command);
  }
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
