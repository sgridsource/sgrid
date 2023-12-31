/* parameters.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 12/99 */

/* basic parameter file syntax:
   parameter_name = parameter_value parameter_value ...
   
   special characters:
   #: this and the rest of the line is ignored
   ": counts as white space

   before interpretation the file content is converted to the form
   par=val ... val par=val ... val ...
*/


#include <ctype.h>
#include "sgrid.h"



/* global vars for parameter data base
   (tParameter is defined in src/main/main/sgrid_main.h) */
tParameter *pdb = NULL;
int pdb_iStart = 0, npdb = 0, npdbmax = 10000;

tParameter *findparameter(char *name, int fatal);
void setparameter(char *name, char *value);
void makeparameter(char *name, char *value, char *description);
void printparameter(tParameter *p);
void printparameters(void);
void translatevalue(char **value);
int set_numericalvalue_byIndex(tParameter *pdb1, int ind, int npdb1max);
int set_booleanvalue_byIndex(tParameter *pdb1, int ind, int npdb1max);



/* parse given parameter file */
void parse_parameter_file(char *parfile)
{
  FILE *fp;
  int c, i, j;
  int nbuffer;
  char *buffer;
  char *par, *val;
  int lpar, lval;

  if (1) printf("Reading parameter file \"%s\"\n", parfile);

  /* read file into memory, add one space at end and beginning */
  fp = fopen(parfile, "r");
  if (!fp)
  {
    printf("Could not open parameter file \"%s\"\n", parfile);
    errorexit("");
  }

  buffer = 0;
  for (i = nbuffer = 0;; i++)
  {
    if (i >= nbuffer-2)
    {
	if (nbuffer > 1000000) 
	  errorexit("Sanity forbids parameter files bigger than 1MB");
	buffer = (char *) realloc(buffer, sizeof(char)*(nbuffer += 1000));
	if (!buffer) errorexit("Out of memory while reading parameter file.");
    }
    if (i == 0) buffer[i++] = ' ';
    if ((c = fgetc(fp)) == EOF) break;
    buffer[i] = c;
  }
  fclose(fp);
  buffer[i++] = ' ';
  buffer[i] = '\0';
  nbuffer = strlen(buffer);

  if (0) { printf("%s", buffer); yo(); }

  /* white out comments and quotes */
  for (i = 0; i < nbuffer; i++) { 
    if (buffer[i] == '#') 
      while (i < nbuffer && buffer[i] != '\n')
	buffer[i++] = ' ';
    if (buffer[i] == '"')
      buffer[i] = ' ';
  }

  /* collapse all white space to single space */
  for (i = j = 1; i < nbuffer; i++) {
    if (!isspace(buffer[i])) buffer[j++] = buffer[i];
    else if (!isspace(buffer[j-1])) buffer[j++] = ' ';
  }
  buffer[j] = '\0';
  nbuffer = strlen(buffer);
  if (0) printf("|%s|\n", buffer);

  /* put a # in front of each par name */
  for (i=0; i<nbuffer-2; i++)
  {
    if(buffer[i]==' ') j=i; /* save pos of newest space found in j */
    if(buffer[i+1]=='=' || buffer[i+2]=='=') /* there's a '=' after pos i */
    {
      buffer[j]='#';  /* put a '#' in front of name, i.e. at pos j */
      i++;
    }
  }
  if (0) printf("|%s|\n", buffer);

  /* and remove spaces around = */
  for (i = j = 1; i < nbuffer; i++) {
    if (buffer[i] != ' ' || (buffer[i-1] != '=' && buffer[i+1] != '=')) 
      buffer[j++] = buffer[i];
  }
  buffer[j] = '\0';
  if(buffer[j-1] == ' ') buffer[j-1] = '\0';
  nbuffer = strlen(buffer);
  if (0) printf("|%s|\n", buffer);

  /* now buffer is 
     |#par=val ... val#par=val ... val| 
  */

  /* split parameter names and values by replacing '#' and '=' by zero */
  for(i=0; i<nbuffer; i++)
  {
    if(buffer[i] == '#' || buffer[i] == '=')
      buffer[i] = '\0';
  }

  /* loop over all parameter/value pairs */
  for (i = 1; i < nbuffer; i += lpar + lval + 2) {
    par = buffer+i;
    lpar = strlen(par);
    val = par + lpar + 1;
    lval = strlen(val);
    if (0) printf("%s = |%s|\n", par, val);

    if (!findparameter(par, 0))
      makeparameter(par, val, "not in libs");
    else
      setparameter(par, val);
  }  

  /* print parameters */
  if (0) {
    printf("after reading the parameterfile:\n");
    printparameters();
  }
  free(buffer);
}





/***************************************************************************/
/* parameter data base */

/* make new parameter in parameter data base, merge if already there */
void makeparameter(char *name, char *value, char *description) 
{
  tParameter *p;

  if (0) printf("Makp %s = %s,  %s\n", name, value, description);

  if(!pdb)
  {
    pdb = (tParameter *) calloc(npdbmax, sizeof(tParameter));
    if (!pdb) errorexit("makeparameter: out of memory");
    npdb = 0;
  }

  p = findparameter(name, 0);

  if (!p) {
    p = &pdb[npdb++];
    p->name = (char *) calloc(strlen(name)+1, sizeof(char));
    p->value = (char *) calloc(strlen(value)+1, sizeof(char));
    strcpy(p->name, name);
    strcpy(p->value, value);
    translatevalue(&p->value);
    set_numericalvalue_byIndex(pdb, npdb-1, npdb);
    set_booleanvalue_byIndex(pdb, npdb-1, npdb);
  } else {
    free(p->description);
  }
  p->description = (char *) calloc(strlen(description)+1, sizeof(char));
  strcpy(p->description, description);

  if (npdb >= npdbmax) 
    errorexit("makeparameter: lazy coding, no more space for new parameters");

  if (0) printparameters();
}




/* set parameter */
void setparameter(char *name, char *value)
{
  tParameter *p = findparameter(name, 1);

  free(p->value);
  p->value = strdup(value);
  translatevalue(&p->value);
  set_numericalvalue_byIndex(p, 0, 1);
  set_booleanvalue_byIndex(p, 0, 1);
}




/* find parameter index */
int findparameterindex(char *name, int fatal)
{
  int i, iS = pdb_iStart;

  if(!name) errorexit("findparameterindex: called without parameter name");

  if( (iS < 0) || (iS >= npdb) ) iS = 0; /* make sure first i is in range */

  for(i = iS; i < npdb; i++)
    if(!strcmp(pdb[i].name, name))
      return i;

  for(i = 0; i < iS; i++)
    if(!strcmp(pdb[i].name, name))
      return i;

  if(fatal) errorexits("Could not find parameter \"%s\"", name);

  /* -1 means par was not found */
  return -1;
}

/* find parameter */
tParameter *findparameter(char *name, int fatal)
{
  int i = findparameterindex(name, fatal);

  if(i<0) return 0;

  return &pdb[i];
}




/* translate parameter value 
   should be much more elaborate, say, implement simple arithmetic
*/
void translatevalue(char **value)
{
  double x = 0;

  if (strcmp(*value, "pi")    == 0) x =    PI;
  if (strcmp(*value, "-pi")   == 0) x =   -PI;
  if (strcmp(*value, "pi/2")  == 0) x =    PI/2;
  if (strcmp(*value, "-pi/2") == 0) x =   -PI/2;
  if (strcmp(*value, "2*pi")  == 0) x =  2*PI;
  if (strcmp(*value, "-2*pi") == 0) x = -2*PI;

  if (x) {
    char newvalue[100];
    sprintf(newvalue, "%.18e", x);
    free(*value);
    *value = strdup(newvalue);
  }
}

/* Write the numerical (double) value of the par with index ind into the 
   par cache. 
   Note: we keep the numerical (double) values of each par in a cache */
int set_numericalvalue_byIndex(tParameter *pdb1, int ind, int npdb1max)
{
  if(pdb1!=NULL && ind>=0 && ind<npdb1max)
    pdb1[ind].numericalvalue = atof(pdb1[ind].value);
  else
    errorexit("set_numericalvalue_byIndex: index out of range");

  return 1;
}

/* Write the boolean (int) value of the par with index ind into the 
   par cache. 
   Note: we keep the boolean values of each par in a cache */
int set_booleanvalue_byIndex(tParameter *pdb1, int ind, int npdb1max)
{
  if(pdb1!=NULL && ind>=0 && ind<npdb1max)
  {
    char *par = pdb1[ind].value;
    int boolval=atoi(par); /* try to read integer value */

    if(strstr(par, "yes")) boolval=1;
    if(strstr(par, "Yes")) boolval=1;
    if(strstr(par, "YES")) boolval=1;
    if(strstr(par, "on")) boolval=1;
    if(strstr(par, "On")) boolval=1;
    if(strstr(par, "ON")) boolval=1;
    if(strstr(par, "true")) boolval=1;
    if(strstr(par, "True")) boolval=1;
    if(strstr(par, "TRUE")) boolval=1;
    
    pdb1[ind].booleanvalue = boolval;
  }
  else
    errorexit("set_booleanvalue_byIndex: index out of range");

  return 1;
}


void printparameter(tParameter *p)
{
  printf("%16s = %-16s,  %s\n                 = %g  ->  %d\n",
  p->name, p->value, p->description,  p->numericalvalue, p->booleanvalue);
}

/* print parameters */
void printparameters(void)
{
  int i;

  for (i = 0; i < npdb; i++)
    printf("pdb[%3d]:  %16s = %-16s,  %s\n", 
	   i, pdb[i].name, pdb[i].value, pdb[i].description);
}






/***************************************************************************/
/* functions for external calls */

/* creation functions */
void AddPar(char *name, char *value, char *description)
{
  makeparameter(name, value, description);
  printf("  parameter %-25s  =  %s\n", name, Gets(name));
}

void AddOrModifyPar(char *name, char *value, char *description)
{
  if (!findparameter(name, 0))
    makeparameter(name, value, description);
  else
    setparameter(name, value);
  printf("  parameter %-25s  =  %s\n", name, Gets(name));
}



/* assignment functions */
void Sets(char *name, char *value)
{
  setparameter(name, value);
}

void Seti(char *name, int i)
{
  char value[100];
  sprintf(value, "%d", i);
  setparameter(name, value);
}

void Setd(char *name, double d)
{
  char value[100];
  sprintf(value, "%.20e", d);
  setparameter(name, value);
}

void Appends(char *name, char *value)
{
  if (Getv(name, value)) return;

  {
    char *oldvalue = Gets(name);
    char *newvalue = cmalloc(strlen(oldvalue) + strlen(value) + 2);
    sprintf(newvalue, "%s %s", oldvalue, value);
    setparameter(name, newvalue);
    free(newvalue);
  }
}




/* query functions */
char *Gets(char *name) 
{
  tParameter *p = findparameter(name, 1);
  return p->value; 
}

char *GetsLax(char *name) 
{
  tParameter *p = findparameter(name, 0);
  return p ? p->value : 0; 
}

int Geti(char *name)
{
  tParameter *p = findparameter(name, 1);
  return atoi(p->value); 
}

double Getd(char *name)
{
  tParameter *p = findparameter(name, 1);
  return atof(p->value);
}

/* "get value" returns 1 if value is in the list of values and 0 else 
   (not equivalent to value being a substring of string parameter) 
*/
int GetvFlag(char *name, char *value, int fatal)
{
  tParameter *p = findparameter(name, fatal);
  char *s=NULL;
  char *parval;
  int lv, ls, lp, startok, endok;

  if (!p) return 0;
  parval = p->value;
  while( (s = strstr(parval, value)) )
  {
    lv = strlen(value);
    //printf("ls=%d  value=%s| s=%s| p->value=%s|\n", ls, value, s, p->value); 
    if( s[lv]==' ' || s[lv]==0 ) break;
    parval = s+1;
  }
  if (!s) return 0;
  ls = strlen(s);
  lp = strlen(p->value);
  startok = (s == p->value || *(s-1) == ' ');  /* how robust is this? */
  endok   = (s+ls == p->value+lp || *(s+ls) == ' ');
  return startok && endok ? 1 : 0;
}

int GetvLax(char *name, char *value)
{
  return GetvFlag(name, value, 0);
}

int Getv(char *name, char *value)
{
  return GetvFlag(name, value, 1);
}





/* "for each"
   return next string delimited by SINGLE space (as in parameter data base) 
   0 if end of list, then restart 
   list = 0 or *list = "" restarts, too
*/
/* how universal is strdup? strsep? */
char *NextEntry(char *list) 
{
  static char *copyoflist = 0;
  static int i = 0, l = 0;
  char *s;

  if (!list || !*list || (i && i == l)) {
    i = l = 0;
    return 0;
  }

  if (!i) {
    //FIX: crashes now: free(copyoflist);  /* no op if null */
    free(copyoflist);
    l = strlen(list);
    copyoflist = (char *) malloc(sizeof(char)*(l+1));
    strcpy(copyoflist, list);
  } else 
    i++;

  s = copyoflist + i;
  for (; i < l && copyoflist[i] != ' '; i++);
  if (copyoflist[i] == ' ') copyoflist[i] = 0;
  
  return s;
}




/* direct index access to parameter data base 
   should not be needed except for, say, dumping the database 
   for checkpointing */
char *GetsInd(int i)
{
  if (0 <= i && i < npdb)
    return pdb[i].value;
  return 0;
}

char *GetnameInd(int i)
{
  if (0 <= i && i < npdb)
    return pdb[i].name;
  return 0;
}

tParameter *GetPointerTo_pbd(void)
{
  return pdb;
}

int GetnParameters(void)
{
  return npdb;
}

/* get the cached numerical value by Index */
double GetCachedNumValByParIndex(int i)
{
  return pdb[i].numericalvalue;
}

/* get the cached boolean value by Index */
int GetCachedBoolValByParIndex(int i)
{
  return pdb[i].booleanvalue;
}

/* get the par Index */
int GetParIndex(char *name)
{
  return findparameterindex(name, 1);
}

/* set the global var pdb_iStart to the index of par name */
int Set_pdb_iStart_AtPar(char *name)
{
  pdb_iStart = findparameterindex(name, 1);
  return pdb_iStart;
}


/**************************************************************************/
int iterate_parameters(int next)
{
  static int ncall = 0;
  tParameter *p;
  char *list, *name, *newvalue, *value;
  char iterpar[100] = "iterate_parameter1";
  char newoutdir[10000], *outdirp;
  int i, j, l;

  /* reset ncall to zero if next=0 */
  if(next==0) { ncall=0;  return 0; }

  /* the default is that we don't want to iterate */
  if (!Getv("iterate_parameters", "yes")) {

    /* return 1 for first call, but 0 for second call, which exits sgrid */
    if (ncall < 0) 
      return 0;
    ncall = -1;
    return 1;
  }

  /* so we want to iterate */
  if(ncall<0) ncall=0;
  printf("\n");
  prdivider(0);
  printf("Iterating parameters:\n");

  sprintf(newoutdir, "%s", Gets("parameterfile"));
  *strstr(newoutdir, ".par") = '\0';

  p = findparameter(iterpar, 0);
  if (!p)
    errorexit("none found, specify iterate_parameter1\n");

  j = 1;
  while (p) {
    list = Gets(iterpar);

    name = NextEntry(list);
    if (!findparameter(name, 0))
      errorexit("iterate_parameterN has to start "
		"with name of existing parameter");
 
    value = NextEntry(list);
    if (!value)
      errorexit("iterate_parameterN needs at least one value");

    i = 0;
    newvalue = 0;
    while (value) {
      if (i == ncall) newvalue = value;
      value = NextEntry(list);
      i++;
    }

    if (newvalue) {
      Sets(name, newvalue);
      printf("%38s = %s\n", name, Gets(name));

      l = strlen(newoutdir);
      sprintf(newoutdir+l, "_%s", newvalue);
      Sets("outdir", newoutdir);
    }

    j++;
    l = strlen("iterate_parameter");
    sprintf(iterpar+l, "%d", j);
    p = findparameter(iterpar, 0);
  }

  if (newvalue) {
    printf("Starting iteration %d\n", ncall);
    if (1) printf("  outdir = %s\n", Gets("outdir"));
    outdirp = (char *) calloc(strlen(newoutdir)+40, sizeof(char));
    strcpy(outdirp, newoutdir);
    strcat(outdirp, "_previous");
    system2("rm -rf", outdirp);
    system3("mv", newoutdir, outdirp); 
    system2("mkdir", newoutdir);
    free(outdirp);
    ncall++;
    return 1;
  }
  printf("  done with iterations.\n");
  prdivider(0);
  return 0;
}

/* print pars with index i1 to i2 in pdb */
void print_pdb_i1_i2(tParameter *pdb, int i1, int i2, int pr_ind, int pr_cache)
{
  int i;

  printf("print_pdb_i1_i2: pdb=%p i1=%d i2=%d\n", pdb, i1,i2);
  for(i=i1; i<=i2; i++)
  {
    if(pr_ind) printf("#%d# ", i);
    if(pr_cache) printf("|%g|%d| ", pdb[i].numericalvalue, pdb[i].booleanvalue);
    printf("%s:\n", pdb[i].description);
    printf("%s = %s\n", pdb[i].name, pdb[i].value);
  }
}
/* print entire parameter database */
void print_parameter_database(void)
{
  print_pdb_i1_i2(pdb, 0, npdb-1, 0,0);
}

/**********************************************/
/* functions to copy the entire par data base */
/**********************************************/

/* Create a copy of parameter database pdb1 in pdb2. 
   This allocates all memory needed for pdb2.
   The caller has to free pdb2 later on its own, e.g. with free_pdb */
/* if we have: tParameter *pdb2;
   call like:  create_copy_of_pdb1_in_pdb2(pdb,npdb,npdbmax, &pdb2); */
void create_copy_of_pdb1_in_pdb2(tParameter *pdb1, int npdb1, int npdb1max,
                                 tParameter **pdb2)
{
  int i;
  tParameter *p2;

  /* allocate array for p2 */
  p2 = (tParameter *) calloc(npdb1max, sizeof(tParameter));
  if(!p2) errorexit("create_copy_of_pdb1_in_pdb2: out of memory");

  /* go over pars in pdb1 and use strdup to create copies in p2 */
  for(i=0; i<npdb1; i++)
  {
    p2[i].name  = strdup(pdb1[i].name);
    p2[i].value = strdup(pdb1[i].value);
    p2[i].description = strdup(pdb1[i].description);
    p2[i].numericalvalue = pdb1[i].numericalvalue;
    p2[i].booleanvalue = pdb1[i].booleanvalue;
    if(!p2[i].name || !p2[i].value || !p2[i].description)
      errorexit("create_copy_of_pdb1_in_pdb2: out of memory");
  }
  /* set pdb2 */
  *pdb2=p2;
}

/* Make an empty par database of max length npdb1max.
   Usage: pdb2 = make_empty_pdb(npdbmax); */
tParameter *make_empty_pdb(int npdb1max)
{
  int i;
  tParameter *pdb1;
  /* allocate array for pdb1 */
  pdb1 = (tParameter *) calloc(npdb1max, sizeof(tParameter));
  if(!pdb1) errorexit("make_empty_pdb: out of memory");

  /* set all entries to NULL */
  for(i=0; i<npdb1max; i++)
  {
    pdb1[i].name  = NULL;
    pdb1[i].value = NULL;
    pdb1[i].description = NULL;
    pdb1[i].numericalvalue = 0.0;
    pdb1[i].booleanvalue = 0;
  }

  return pdb1;
}

/* copy pdb1 into pdb2 */
void copy_pdb(tParameter *pdb1, int npdb1, tParameter *pdb2)
{
  int i;

  /* Go over pars in pdb1. Use realloc and strcpy to create copies in pdb2 */
  for(i=0; i<npdb1; i++)
  {
    pdb2[i].name = 
     (char *) realloc(pdb2[i].name, 
                      sizeof(char)*(strlen(pdb1[i].name)+1));
    pdb2[i].value = 
     (char *) realloc(pdb2[i].value, 
                      sizeof(char)*(strlen(pdb1[i].value)+1));
    pdb2[i].description = 
     (char *) realloc(pdb2[i].description, 
                      sizeof(char)*(strlen(pdb1[i].description)+1));
    if(!pdb2[i].name || !pdb2[i].value || !pdb2[i].description)
      errorexit("copy_pdb: out of memory");
    strcpy(pdb2[i].name, pdb1[i].name);
    strcpy(pdb2[i].value, pdb1[i].value);
    strcpy(pdb2[i].description, pdb1[i].description);
    pdb2[i].numericalvalue = pdb1[i].numericalvalue;
    pdb2[i].booleanvalue = pdb1[i].booleanvalue;
  }
}

/* free the parameter database content in pdb1 */
void free_pdb_contents(tParameter *pdb1, int npdb1)
{
  int i;

  /* go over pars in pdb1 and free them */
  for(i=0; i<npdb1; i++)
  {
    free(pdb1[i].name);
    free(pdb1[i].value);
    free(pdb1[i].description);
  }
}

/* free the parameter database pdb1 */
void free_pdb(tParameter *pdb1, int npdb1)
{
  free_pdb_contents(pdb1, npdb1);
  free(pdb1);
}

/* free the global parameter database pdb and set npdb=0 */
void free_global_parameter_database_contents(void)
{
  free_pdb_contents(pdb, npdb);
  npdb = 0;
  pdb_iStart = 0;
}

/* free the global parameter database pdb and set npdb=0 */
void free_global_parameter_database(void)
{
  free_global_parameter_database_contents();
  free(pdb);
  pdb = NULL;
}
