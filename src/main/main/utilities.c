/* utilities.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99 */

#include "sgrid.h"
#include "main.h"




/* debug */
void yo() {fprintf(stdout, "Yo!\n");fflush(stdout);}
void yo1() {fprintf(stdout, "Yo1!\n");fflush(stdout);}
void yo2() {fprintf(stdout, "Yo2!\n");fflush(stdout);}
void yo3() {fprintf(stdout, "Yo3!\n");fflush(stdout);}




/* pretty */
void prdivider(int n)
{
  printf("------------------------------------------------------------------------------\n");
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




/* ugh, but how universal are those built in functions? */
int system2(char *s1, char *s2) 
{
  return system3(s1, s2, "");
}
int system3(char *s1, char *s2, char *s3) 
{
  char command[10000];
  int status = 0;

  if (1) { 
    sprintf(command, "%s %s %s", s1, s2, s3);
    printf("System call:  %s\n", command);
    status = system(command);
  } else {
    sprintf(command, "%s %s %s", s1, s2, s3);
    //if (1) printf("System call for proc %d:  %s\n", sgridpi_rank(), command);
    status = system(command);
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
