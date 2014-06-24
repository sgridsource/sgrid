/* extractvar_from_checkpoint.c */
/* transform or modify checkpoint data */
/* (c) Wolfgang Tichy 2014 */


#define STRLEN 262144

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "../src/utility/checkpoint/wolfio.h"
#include "../src/utility/checkpoint/wolfio.c"



/* functions */

/* read a var from checkpoint */
int readvariable(FILE *in, char *infostr, int *ndata, double **data)
{
  char *str;

  fscanline(in, infostr);
  /* return 0 if infostr is empty */
  if(strstr(infostr, "$variable")==NULL) return 0;

  /* find number of doubles in var */
  str=strstr(infostr, "->nnodes = ");
  str += 11; /* go 11 chars forward */
  *ndata = atoi(str);

  /* realloc data to make space */
  *data = (double *) realloc(*data, sizeof(double *) * (*ndata));
    
  /* read var data */
  fread(*data, sizeof(double), *ndata, in);
  fscanline(in, infostr+strlen(infostr)+1); /* read \n after var data */

  return 1;
}

/* write a single var */
void dumpvariable(FILE *out, char *name, char *time, int boxi,
                  int n1, int n2, int n3, double *data)
{
  int i;
  int N= n1*n2*n3;
  float dat;

  fprintf(out,
    "# \"time = %s\" , box%d: n1=%d n2=%d n3=%d , N=%d numbers , type=float\n",
    time, boxi, n1,n2,n3, N);
  for(i=0;i<N;i++)
  {
    dat = data[i];
    fwrite(&dat, sizeof(dat), 1, out);
  }
  fprintf(out, "\n");
}


/* main prog. */
int main(int argc, char *argv[])
{
  FILE *in;
  FILE *out;
  char str[STRLEN];
  char str1[STRLEN];
  char str2[STRLEN];
  char str3[STRLEN];
  char infostr[STRLEN];
  char time[STRLEN];
  char bn1[STRLEN];
  char bn2[STRLEN];
  char bn3[STRLEN];
  double *data=NULL;
  int i, ndata, n1,n2,n3, boxi;
  char *astr;
  char *varname;

  /* mem for var data */
  data = (double *) calloc(10, sizeof(double));
  
  printf("# extractvar_from_checkpoint \n");
  if(argc!=6)  
  {
   printf("# usage: extractvar_from_checkpoint -V vname b  checkpoint.0 v.XYZ\n");
   printf("# options: -V vname b  extract var \"vname\" from box b\n");
   printf("# examples:\n");
   printf("# extractvar_from_checkpoint -V gxx 1  checkpoint.0 gxx.XYZ1\n");
   return -1;
  }

  /* default options */
  varname=NULL;

  /* parse command line options, which start with - */
  for(i=1; (i<argc)&&(argv[i][0] == '-'); i++)
  {
   astr = argv[i];

   if( (strcmp(astr+1,"V")==0) )
   {
     if(i>=argc-2) 
     {
       printf("no varname and box number after -V\n");
       return -1;
     }
     varname=argv[i+1];
     i++;
     boxi=atoi(argv[i+1]);
     i++;
   }
   else 
   {
     printf("extractvar_from_checkpoint: unknown argument %s\n", astr);
     return -1;
   }
  }

  //printf("i=%d argc=%d\n", i, argc);

  /* put par about box dims in strings */
  sprintf(bn1, "box%d_n1", boxi);
  sprintf(bn2, "box%d_n2", boxi);
  sprintf(bn3, "box%d_n3", boxi);

  /* info */
  if(varname!=NULL)
    printf("# reading variable %s from box%d\n", varname, boxi);
  else return -1;

  /* open file in */
  printf("# input file: %s",argv[i]);
  in=fopen(argv[i],"r");
  if(in==NULL)
  {
   printf(" not found.\n");
   return -2;
  }
  printf("\n");

  /* open file out */
  printf("# output file: %s",argv[i+1]);
  out=fopen(argv[i+1],"wb");
  if(out==NULL)
  {
   printf(" could not be opened.\n");
   return -2;
  }
  printf("\n");

  fgetparameter(in, "$time", time);

  /* find $BEGIN_variables:  */
  while(fscanline(in,str)!=EOF)
  {
    //printf("|%s|\n", str);
    /* break if we get to BEGIN_variables: */
    if(strstr(str, "$BEGIN_variables:")!=NULL) break;

    /* read 3 box pars */
    if(strstr(str, bn1)==str)
    {
      astr = str + strlen(bn1)+3; /* let astr point to parval */
      n1=atoi(astr);
    }
    if(strstr(str, bn2)==str)
    {
      astr = str + strlen(bn1)+3; /* let astr point to parval */
      n2=atoi(astr);
    }
    if(strstr(str, bn3)==str)
    {
      astr = str + strlen(bn1)+3; /* let astr point to parval */
      n3=atoi(astr);
    }
  }

  /* read, modify and write a variable */
  while(readvariable(in, infostr, &ndata, &data))
  {
    int bi;
    //printf("%s\nndata=%d data[0]=%g\n", infostr, ndata, data[0]);
    /* get variable name in str1 and box ind in str2 */
    sscanf(infostr, "%s%s%s%s%s%s%s", str3, str3, str1, str3, str3, str3, str2);
    bi = atoi(str2);

    /* dump var if it's varname */
    if(strcmp(str1, varname)==0)
      if(bi=boxi)
      {
        dumpvariable(out, varname, time, boxi, n1,n2,n3, data);
        break;
      }
  }

  /* close files, free memory and return */
  fclose(in);
  fclose(out);
  free(data); /* free data allocated by readvariable */

  return 0; 
}
