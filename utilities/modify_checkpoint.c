/* modify_checkpoint.c */
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
      
  return 1;
}

/* write a single var */
void writevariable(FILE *out, char *infostr, int ndata, double *data)
{
  fprintf(out, "%s\n", infostr);
  fwrite(data, sizeof(double), ndata, out);
  fprintf(out, "\n");
}


/* main prog. */
int main(int argc, char *argv[])
{
  int tcol;
  FILE *in;
  FILE *out;
  char str[STRLEN];
  char str1[STRLEN];
  char str2[STRLEN];
  char str3[STRLEN];
  char infostr[STRLEN];
  double *data=NULL;
  int i;
  char *astr;
  double number, add, mul;
  char *varname;

  /* mem for var data */
  data = (double *) calloc(10, sizeof(double));
  
  printf("# modify_checkpoint \n");
  if(argc<3 || argc>9)  
  {
   printf("# usage: modify_checkpoint [-V name] [-M m] [-A a] in.0 out.0\n");
   printf("# options: -V name  modify variable \"name\"\n");
   printf("#          -M m     multiply by \"m\"\n");
   printf("#          -A a     add \"a\"\n");
   printf("# examples:\n");
   printf("# modify_checkpoint -V gxx -A 25.13274122871834 checkpoint.0 l\n");
   return -1;
  }

  /* default options */
  add=0.0;
  mul=1.0;
  varname=NULL;

  /* parse command line options, which start with - */
  for(i=1; (i<argc)&&(argv[i][0] == '-'); i++)
  {
   astr = argv[i];

   if( (strcmp(astr+1,"V")==0) )
   {
     if(i>=argc-1) 
     {
       printf("no varname after -V\n");
       return -1;
     }
     varname=argv[i+1];
     i++;
   }
   else if( (strcmp(astr+1,"A")==0) )
   {
     if(i>=argc-1) 
     {
       printf("argument a is needed after -A\n");
       return -1;
     }
     add=atof(argv[i+1]);
     i++;
   }
   else if( (strcmp(astr+1,"M")==0) )
   {
     if(i>=argc-1) 
     {
       printf("argument m is needed after -M\n");
       return -1;
     }
     mul=atof(argv[i+1]);
     i++;
   }
   else 
   {
     printf("modify_checkpoint: unknown argument %s\n", astr);
     return -1;
   }
  }

  //printf("i=%d argc=%d\n", i, argc);

  /* info about colmuns */
  if(varname!=NULL) printf("# modifying variable %s\n", varname);

  
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

  /* find $BEGIN_variables:  */
  while(fscanline(in,str)!=EOF)
  {
    /* break if we get to BEGIN_variables: */
    if(strstr(str, "BEGIN_variables:")!=NULL) break;
    
    /* modify pars ??? */
    // if(strstr(str, "parname")==str)  do_someting
  }

  /* read, modify and write a variable */
  while(readvariable(in, infostr, &ndata, &data))
  { 
    /* get variable name in str3 */
    sscanf(infostr, "%s%s%s", str1, str2, str3);
    /* modify var if it's varname */
    if(strcmp(str3, varname)==0)
      for(i=0; i<ndata; i++)
      {
        /* multiply */
        if(mul!=1.0) data[i] *= mul;
        /* add */
        if(add!=0.0) data[i] += add;
      }

    /* write var into file out */
    writevariable(out, infostr, ndata, data);
  }

  /* write $END_variables */
  fprintf(out, "$END_variables\n\n", );

  /* close files, free memory and return */
  fclose(in);
  fclose(out);
  free(data); /* free data allocated by readvariable */

  return 0; 
}
