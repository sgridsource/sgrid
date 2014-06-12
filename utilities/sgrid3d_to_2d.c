/* sgrid3d_to_2d.c,  (c) Wolfgang Tichy 2006 */
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
e.g. sgrid3d_to_2d -j=5 -k=0 foo
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
  int ijk, i,j,k, l1,l1m, l2,l2m, n1,n2,n3, npoints;
  float *data=NULL;
  FILE *file;
  int mode;

  if(argc != 3)
  {
    printf("Usage: sgrid3d_to_2d [-i=<num>] [-j=<num>] [-k=<num>] filename\n");
    printf("[note: you need exactly 1 option argument]\n");
    return -10;
  }
  if(argv[1][1]=='i')
  { i=atoi(argv[1]+3);  mode=3;}
  else if(argv[1][1]=='j')
  { j=atoi(argv[1]+3);  mode=2;}
  else if(argv[1][1]=='k')
  { k=atoi(argv[1]+3);  mode=1;}
  else
  { printf("you need exactly 1 option argument\n"); return -2; }

  //printf("filename = argv[i] = %s\n", argv[i]);
  /* open file */
  file=fopen(argv[2], "r");
  if(file==NULL)
  {
    printf("%s not found. \n", argv[2]);
    return -1;
  }
  while(fgetparameter(file, "\"time", str)!=EOF)
  {
    str[strlen(str)-1]=0;
    time=atof(str);
    printf("# \"time=%f\"\n", time);
  
    fgetparameter(file, "n1", str);  n1=atoi(str);
    fgetparameter(file, "n2", str);  n2=atoi(str);
    fgetparameter(file, "n3", str);  n3=atoi(str);
    //printf("# n1,n2,n2 = %d,%d,%d\n", n1,n2,n3);
    npoints=n1*n2*n3;
    fscanline(file, str); /* read til end of line */
    
    /* read the actual data */
    if(data==NULL) data = calloc(npoints, sizeof(data));
    fread(data, sizeof(*data), npoints, file);
  
    if(mode==1)      { l1m=n1;  l2m=n2; }
    else if(mode==2) { l1m=n1;  l2m=n3; }
    else             { l1m=n2;  l2m=n3; }
    for(l2=0; l2<l2m; l2++)
    {
      for(l1=0; l1<l1m; l1++)
      {
        if(mode==1)      { i=l1;  j=l2; }
        else if(mode==2) { i=l1;  k=l2; }
        else             { j=l1;  k=l2; }
        if(i<0 || i>=n1) { printf("i=%d is not allowed!\n",i); exit(11); }
        if(j<0 || j>=n2) { printf("j=%d is not allowed!\n",j); exit(11); }
        if(k<0 || k>=n3) { printf("k=%d is not allowed!\n",k); exit(11); }
        ijk = ((i)+n1*((j)+n2*(k)));
        printf("%d %d %g\n", l1,l2, data[ijk]);
      }
      printf("\n");
    }
    printf("\n");
  }

  fclose(file);
  free(data);
  return 0;
}
