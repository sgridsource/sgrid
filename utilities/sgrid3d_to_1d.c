/* sgrid3d_to_1d.c,  (c) Wolfgang Tichy 2006 */
/* dump var along a line: 
  here files have the following format: */
/*
# "time = 0" , box0: n1=30 n2=20 n3=18 , N=10800 numbers , type=float
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# "time = 2" , box0: n1=30 n2=20 n3=18 , N=10800 numbers , type=float
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# "time = 4" , box0: n1=30 n2=20 n3=18 , N=10800 numbers , type=float
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

*/
/* where @@@@... stands for 10800 floats.
e.g. sgrid3d_to_1d -j=5 -k=0 foo
will dump the file foo along the line j=5,k=0 */

#define STRLEN 65535

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <endian.h>

/* my functions */

/* my wolfio functions */
#include "../src/utility/checkpoint/wolfio.h"
#include "../src/utility/checkpoint/wolfio.c"


/* main program */
int main(int argc, char **argv)
{
  char str[STRLEN];
  double time;
  int ijk, i,j,k, l,lm, n1,n2,n3, npoints;
  float *data=NULL;
  FILE *file;
  int mode;

  if(argc != 4)
  {
    printf("Usage: sgrid3d_to_1d [-i=<num>] [-j=<num>] [-k=<num>] filename\n");
    printf("[note: you need exactly 2 option arguments]\n");
    return -10;
  }
  if(argv[1][1]=='i' && argv[2][1]=='j')
  { i=atoi(argv[1]+3); j=atoi(argv[2]+3);  mode=3;}
  else if(argv[1][1]=='i' && argv[2][1]=='k')
  { i=atoi(argv[1]+3); k=atoi(argv[2]+3);  mode=2;}
  else if(argv[1][1]=='j' && argv[2][1]=='k')
  { j=atoi(argv[1]+3); k=atoi(argv[2]+3);  mode=1;}
  else
  { printf("you need exactly 2 different option arguments\n"); return -2; }

  //printf("filename = argv[i] = %s\n", argv[i]);
  /* open file */
  file=fopen(argv[3], "r");
  if(file==NULL)
  {
    printf("%s not found. \n", argv[3]);
    return -1;
  }
  while(fgetparameter(file, "\"time", str)!=EOF)
  {
    str[strlen(str)-1]=0;
    time=atof(str);
    printf("\"time=%f\n", time);
  
    fgetparameter(file, "n1", str);  n1=atoi(str);
    fgetparameter(file, "n2", str);  n2=atoi(str);
    fgetparameter(file, "n3", str);  n3=atoi(str);
    //printf("# n1,n2,n2 = %d,%d,%d\n", n1,n2,n3);
    npoints=n1*n2*n3;
    fscanline(file, str); /* read til end of line */
    
    /* read the actual data */
    if(data==NULL) data = calloc(npoints, sizeof(data));
    fread(data, sizeof(*data), npoints, file);
  
    if(mode==1)      lm=n1;
    else if(mode==2) lm=n2;
    else             lm=n3;
    for(l=0; l<lm; l++)
    {
      if(mode==1)      i=l;
      else if(mode==2) j=l;
      else             k=l;
      if(i<0 || i>=n1) { printf("i=%d is not allowed!\n",i); exit(11); }
      if(j<0 || j>=n2) { printf("j=%d is not allowed!\n",j); exit(11); }
      if(k<0 || k>=n3) { printf("k=%d is not allowed!\n",k); exit(11); }
      ijk = ((i)+n1*((j)+n2*(k)));
      printf("%d  %g\n", l, data[ijk]);
    }
    printf("\n");
  }

  fclose(file);
  free(data);
  return 0;
}
