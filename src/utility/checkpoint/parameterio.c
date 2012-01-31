/* parameterio.c */
/* Wolfgang Tichy 12/2008 */

/* output current pars and update them from files */

#include "sgrid.h"
#include "checkpoint.h"



/* write all current pars into a file */
int parameterio_write_current_pars(tGrid *g)
{
  char *curparfile = "sgrid_current_parameters.par";
  char *filename;
  FILE *fp;
  int i;

  /* name of file for current pars */
  filename = cmalloc(strlen(Gets("outdir")) + strlen(curparfile) + 20);
  sprintf(filename, "%s/%s", Gets("outdir"), curparfile);

  /* open file for writing */
  fp = fopen(filename, "wb");
  if (!fp) errorexits("failed opening %s", filename);

  /* write parameter database */
  for (i = 0; i < GetnParameters(); i++)
  {
    fprintf(fp, "%s = ", GetnameInd(i));
    fprintf(fp, "%s\n", GetsInd(i));
  }
  fclose(fp);
  free(filename);
  return 0;
}

/* overwrite sgrid pars with the values in file sgrid_update_parameters.par */
int parameterio_update_pars(tGrid *g)
{
  char *updatefile = "sgrid_update_parameters.par";
  char *filename;
  FILE *fp;
  int i, j, ret;
  char name[10000], value[10000], str[10000], *currentvalue;
  int parsread;

  /* name of file with pars we update */
  filename = cmalloc(strlen(Gets("outdir")) + strlen(updatefile) + 20);
  sprintf(filename, "%s/%s", Gets("outdir"), updatefile);

  /* open file for reading */
  fp = fopen(filename, "rb");
  if (!fp) return 0; /* errorexits("failed opening %s", filename); */

  /* read parameters */
  if(1) printf("parameterio_update_pars:\n");
  parsread = 0;
  while( (ret=fscanline(fp, str)) != EOF )
  {
    while(ret==0)  ret = fscanline(fp, str);
    if( ret==EOF ) break;

    if(0) printf("  str: %s\n", str);
    extrstr_before_after_EQ(str, name, value);
 
    /* remove spaces at the end of name and the beginning of value */
    strcpy(str, name);
    sscanf(str, "%s", name);
    strcpy(str, value);
    j = strspn(str, " ");  /* strcpy(value, str+j); */
    if(j>0) strcpy(value, str+1); /* <--we remove at most one space */
    if(0) printf("  found parameter %s = ", name);
    if(0) printf("%s\n", value);
    currentvalue = GetsLax(name);
    if(currentvalue==0)
    {
      printf("  skipping %s = %s\n", name, value);
      continue;
    }
    if(strcmp(value, currentvalue) != 0)
    {
      if(1) printf("  overwriting %s = %s -> %s\n",
		      name, currentvalue, value);
      Sets(name, value);
      parsread++;
    }
  }
  printf("  updated %d parameters\n", parsread);
  fclose(fp);

  /* rename sgrid_update_parameters.par */  
  sprintf(str, "%s_done", filename);
  system3("mv -f", filename, str);

  free(filename);
  return 1;
}
