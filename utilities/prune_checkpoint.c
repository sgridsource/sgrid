/* prune_checkpoint.c */
/* write a pruned checkpint that contains only the vars listed in the
   vars_to_keep.txt file */
/* (c) Wolfgang Tichy 2023 */


#define STRLEN 262144
#define NVARSMAX 10000

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
void writevariable(FILE *out, char *infostr, int ndata, double *data)
{
  fprintf(out, "%s\n", infostr);
  fwrite(data, sizeof(double), ndata, out);
  fprintf(out, "\n");
}


/* main prog. */
int main(int argc, char *argv[])
{
  FILE *vin;
  FILE *in;
  FILE *out;
  char str[STRLEN];
  char str1[STRLEN];
  char str2[STRLEN];
  char str3[STRLEN];
  char infostr[STRLEN];
  double *data=NULL;
  int i, ndata;
  char *astr;
  double add, mul;
  char *varname, *parname, *parval;
  char *var_to_keep[NVARSMAX];
  int nvars_to_keep;

  /* mem for var data */
  data = (double *) calloc(10, sizeof(double));
  
  printf("# prune_checkpoint \n");
  if(argc != 4)  
  {
   printf("# usage: prune_checkpoint vars_to_keep.txt in.0 ou.0\n");
   printf("# examples:\n");
   printf("# prune_checkpoint vars_to_keep.txt checkpoint.0 out.0\n");
   return -1;
  }

  /* open file with vars we keep */
  printf("# var file: %s",argv[1]);
  vin=fopen(argv[1],"r");
  if(vin==NULL)
  {
   printf(" not found.\n");
   return -2;
  }
  printf("\n");

  /* read vars we keep */
  nvars_to_keep = 0;
  while(fscanline(vin,str)!=EOF)
  {
    sscanf(str, "%s", str1);
    var_to_keep[nvars_to_keep] = strdup(str1);
    nvars_to_keep++;
  }
  fclose(vin);
  if(nvars_to_keep<=0)
  {
   printf(" no vars found in file.\n");
   return -3;
  }
  printf("# we keep these vars:\n");
  for(i=0; i<nvars_to_keep; i++)
  {
    printf("# %s\n", var_to_keep[i]);
  }
  printf("\n");

  /* open file in */
  printf("# input file: %s",argv[2]);
  in=fopen(argv[2],"r");
  if(in==NULL)
  {
   printf(" not found.\n");
   return -2;
  }
  printf("\n");

  /* open file out */
  printf("# output file: %s",argv[3]);
  out=fopen(argv[3],"wb");
  if(out==NULL)
  {
   printf(" could not be opened.\n");
   return -2;
  }
  printf("\n");

  /* find $BEGIN_variables:  */
  while(fscanline(in,str)!=EOF)
  {
    /* add line with str to file out */
    fprintf(out, "%s\n", str);

    /* break if we get to BEGIN_variables: */
    if(strstr(str, "$BEGIN_variables:")!=NULL)
    {
      break;
    }
  }

  /* read a variable and write it to output if we keep it */
  while(readvariable(in, infostr, &ndata, &data))
  {
    //printf("%s\nndata=%d data[0]=%g\n", infostr, ndata, data[0]);
    /* get variable name in str3 */
    sscanf(infostr, "%s%s%s", str1, str2, str3);

    /* check if var in str3 is in var_to_keep */
    for(i=0; i<nvars_to_keep; i++)
    {
      if(strcmp(str3, var_to_keep[i])==0)
      {
        /* write var into file out */       
        writevariable(out, infostr, ndata, data);
        break;
      }
    }
  }

  /* write $END_variables */
  fprintf(out, "$END_variables\n\n");

  /* close files, free memory and return */
  fclose(in);
  fclose(out);
  free(data); /* free data allocated by readvariable */

  for(i=0; i<nvars_to_keep; i++) free(var_to_keep[i]);

  return 0; 
}
