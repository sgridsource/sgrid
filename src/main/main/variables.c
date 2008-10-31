/* variables.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 12/99 */


#include "sgrid.h"


typedef struct {
  double iotime[4];
  int ioiter[4];
  int ioflag[4];
  int ioatall;
} tIO;


typedef struct tVAR {
  char *name;
  char *tensorindices;
  char *description;
  int index;
  int ncomponents;
  int component;
  tIO *io;
  double farlimit;
  double falloff;
  double propspeed;
  int sym[3];
  int constant;
} tVar;

tVar *vdb = 0;
int nvdb = 0, globalnvariables = 0;




/* add variable to data base */
void AddVar(char *name, char *tensorindices, char *description)
{
  tVar *new;
  int i, j;
  char fullname[100];
  int nilist;
  char *ilist[NINDEXLIST]; 
  int sym[3*NINDEXLIST];
  char *symsigns[3] = {"-", "0", "+"}, **ss = symsigns+1;

  if (0) printf("AddVar: name %s, tensorindices %s\n", name, tensorindices);

  /* construct list of tensor indices */
  tensorindexlist(tensorindices, &nilist, ilist, sym);

  /* for each tensor index */
  for (j = 0; j < nilist; j++) {

    /* construct name of variable */
    snprintf(fullname, 100, "%s%s", name, ilist[j]);
    printf("  variable  %s\n", fullname);
    if (0) printf("%13s%s%s\n", ss[sym[3*j]], ss[sym[3*j+1]], ss[sym[3*j+2]]);

    /* make sure that this variable does not exist yet */
    for (i = 0; i < nvdb; i++) 
      if (!strcmp(vdb[i].name, fullname))
	errorexits("AddVar: variable \"%s\" already exists\n", fullname);

    /* variable does not exist, so add a new element to data base */
    vdb = (tVar *) realloc(vdb, sizeof(tVar)*(nvdb+1));
    new = &vdb[nvdb];

    /* initialize and fill in structure */
    memset(new, 0, sizeof(tVar));
    new->name          = strdup(fullname);
    new->tensorindices = strdup(tensorindices);
    new->description   = strdup(description);
    new->index         = nvdb;
    new->ncomponents   = nilist;
    new->component     = j;
    new->io            = NULL;
    new->sym[0]        = sym[3*j];
    new->sym[1]        = sym[3*j+1];
    new->sym[2]        = sym[3*j+2];

    nvdb++;
    globalnvariables = nvdb;
  }
}




/* add constant variable to data base */
void AddConstantVar(char *name, char *tensorindices, char *description)
{
  AddVar(name, tensorindices, description);
  VarNameSetConstantFlag(vdb[nvdb-1].name);
}




/* return index of variable */
int Ind(char *name) 
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (!strcmp(vdb[i].name, name)) {
      if (0) printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }
  errorexits("Ind: variable \"%s\" does not exist\n", name);
  return 0;
}

/* return index of variable or -1 if it was not found */
int IndLax(char *name) 
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (!strcmp(vdb[i].name, name)) {
      if (0) printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }
  return -1;
}

/* return index of variable given pointer */
int IndFromPtr(tBox *box, double *p) 
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (box->v[i] == p)
      return vdb[i].index;
  return -1;
}

/* return name given index */
char *VarName(int i) 
{
  if (i < 0 || i >= nvdb)
    errorexit("VarName: index out of range");

  return vdb[i].name;
}

/* return number of components */
int VarNComponents(int i)
{
  if (i < 0 || i >= nvdb)
    errorexit("VarNComponents: index out of range");
  if (vdb[i].component != 0)
  {
    /* errorexit("VarNComponents: you have to use index of zeroth component"); */
    i = IndComponent0(i);
  }
  return vdb[i].ncomponents;
}

/* return component */
int VarComponent(int i)
{
  return vdb[i].component;
}

/* return index of component 0 */
int IndComponent0(int i)
{
  return i - vdb[i].component;
}

/* return name of component 0 for a given name */
char *VarNameComponent0(char *name)
{
  return VarName(IndComponent0(Ind(name)));
}

/* return string with tensor indices */
char *VarTensorIndices(int i)
{
  return vdb[i].tensorindices;
}

/* set information on how variable behaves at Boundary*/
void VarNameSetBoundaryInfo(char *name, 
			    double farlimit, double falloff, double propspeed)
{
  int i = Ind(name);

  vdb[i].farlimit = farlimit;
  vdb[i].falloff = falloff;
  vdb[i].propspeed = propspeed;
}

/* set information on how variable behaves at Boundary*/
void VarNameSetConstantFlag(char *name)
{
  int i, i0 = IndComponent0(Ind(name));
  int n  = VarNComponents(i0);

  for (i = 0; i < n; i++) {
    vdb[i+i0].constant = 1;
    if (0) printf("  setting %s constant\n", vdb[i+i0].name);
  }
}

/* return various boundary information */
double VarFallOff(int i) {return vdb[i].falloff;}
double VarFarLimit(int i) {return vdb[i].farlimit;}
double VarPropSpeed(int i) {return vdb[i].propspeed;}
int VarSymmetry(int i, int dir) {return vdb[i].sym[dir];}
int VarConstantFlag(int i) {return vdb[i].constant;}





/************************************************************************/
/* utility functions for variable lists */


/* print variable list */
void prvarlist(tVarList *v)
{
  int i, j;

  printf("VarList=%p  grid=%p  time=%g  n=%d\n", v, v->grid, v->time, v->n);
  for (i = 0; i < v->n; i++)
  {
    j = v->index[i];
    printf(" %d  VarIndex=%d  %s\n", i, j, VarName(j));
  }
}

/* print variable list in one box */
void prvarlist_inbox(tBox *box, tVarList *v)
{
  int i, j;

  printf("box%d=%p: VarList=%p  grid=%p  time=%g  n=%d\n",
          box->b, box, v, v->grid, v->time, v->n);
  for (i = 0; i < v->n; i++)
  {
    j = v->index[i];
    printf(" %d  VarIndex=%d  box->v[%d]=%p  %s\n",
           i, j, j, box->v[j], VarName(j));
  }
}


/* allocate an empty variable list */
tVarList *vlalloc(tGrid *grid)
{
  tVarList *u;

  u = calloc(1, sizeof(tVarList));
  u->grid = grid;
  if(grid) u->time = grid->time;
  return u;
} 




/* free a variable list */
void vlfree(tVarList *u)
{
  if (u) {
    if (u->index) free(u->index);
    free(u);
  }
} 




/* add a variable (one component) to a variable list */
void vlpushone(tVarList *v, int vi)
{
  v->n += 1;
  v->index = realloc(v->index, sizeof(int) * v->n); 
  v->index[v->n-1] = vi;
}




/* add a variable with all its components to a variable list */
void vlpush(tVarList *v, int vi)
{
  int i, n = VarNComponents(vi);

  if(IndComponent0(vi)!=vi) 
    errorexit("vlpush (23.8.2008): vi needs to be index of component 0. "
              "Consider using vlpushone.");
  v->n += n;
  v->index = realloc(v->index, sizeof(int) * v->n); 
  for (i = 0; i < n; i++)
    v->index[v->n-n+i] = vi + i;
}




/* add a variable list to a variable list */
void vlpushvl(tVarList *v, tVarList *u)
{
  int i;

  if (!v || !u) return;
  v->n += u->n;
  v->index = realloc(v->index, sizeof(int) * v->n); 
  for (i = 0; i < u->n; i++)
    v->index[v->n - u->n + i] = u->index[i];
}




/* drop a variable (one component) from a variable list */
void vldropone(tVarList *v, int vi)
{
  int i;
  
  for (i = 0; i < v->n; i++)
    if (v->index[i] == vi) {
      v->n -= 1;
      for (; i < v->n; i++)
	v->index[i] = v->index[i+1];
      break;
    }
}

/* drop a variable with all its components from a variable list */
void vldrop(tVarList *v, int vi)
{
  int i, n = VarNComponents(vi);

  for (i = 0; i < n; i++)
    vldropone(v, vi+i);
}

/* drop last n variables from a variable list */
void vldropn(tVarList *v, int n)
{
  if (n <= 0) 
    return;
  if (n >= v->n) 
    v->n = 0;
  else
    v->n -= n;
}




/* duplicate variable list */
tVarList *vlduplicate(tVarList *v)
{
  int i;
  tVarList *u = vlalloc(v->grid);

  u->time = v->time;

  for (i = 0; i < v->n; i++) 
    vlpushone(u, v->index[i]);

  return u;
}




/* enable all variables in a variable list */
void vlenable(tVarList *v)
{
  enablevarlist(v);
}

void vlenablegrid(tGrid *grid, tVarList *v)
{
  v->grid = grid;
  enablevarlist(v);
}




/* disable all variables in a variable list */
void vldisable(tVarList *v)
{
  disablevarlist(v);
}




/* create, enable, return pointer for a 1 variable VarList */
tVarList *VLPtrEnable1(tGrid *grid, char *varname)
{
  tVarList *vl = vlalloc(grid);
  int i = Ind(varname);
  
  enablevar(grid, i);
  vlpush(vl, i);
  return vl;
}




/* disable variables in a VarList and free VarList */
void VLDisableFree(tVarList *vl)
{
  disablevarlist(vl);
  vlfree(vl);
}




/* add variables based on an existing variable list and a postfix 
   note that we add each component as a scalar but fix it later because
   we want gxx_p, gxy_p, ...  and not gxx_pxx, gxx_pxy ...
*/
tVarList *AddDuplicate(tVarList *vl, char *postfix) 
{
  char name[1000];
  int i, j;
  int nadded = 0;
  tVarList *newvl;
  tVar *var, *newvar;

  /* new variable list with same number of indices */
  newvl = vlduplicate(vl);

  /* for all scalar variables in list */
  for (i = 0; i < vl->n; i++) {

    /* construct new name */
    var = &vdb[vl->index[i]];
    snprintf(name, 1000, "%s%s", var->name, postfix);

    /* if variable already exists, don't add it again */
    /* note that we nevertheless return a corresponding variable list */
    if ((j = IndLax(name)) >= 0) {
      newvl->index[i] = j;
      continue;
    }

    /* add scalar variable with new name to variable database */
    AddVar(name, "", var->description);
    nadded++;

    /* get index of new variable and overwrite index in duplicate */
    newvl->index[i] = Ind(name);

    /* get pointer to old variable again since AddVar reallocates vdb */
    var = &vdb[vl->index[i]];

    /* set structure in variable data base */
    newvar = &vdb[newvl->index[i]];
    free(newvar->tensorindices);
    newvar->tensorindices = strdup(var->tensorindices);
    newvar->component     = var->component;
    newvar->ncomponents   = var->ncomponents;
    newvar->farlimit      = var->farlimit;
    newvar->falloff       = var->falloff;
    newvar->propspeed     = var->propspeed;
    newvar->constant      = var->constant;
    for (j = 0; j < 3; j++)
      newvar->sym[j] = var->sym[j];
  }

  /* create storage for as many variables as have been actually added 
     do it on all boxs so that nvariables remains the same on all boxs
  */
  if (vl->grid && nadded) {
    tGrid *g = vl->grid;
    int n = g->nvariables + nadded;
    realloc_gridvariables(g, n);
  }
  if (0) printf("nvdb is now %d\n", nvdb);
  
  return newvl;
}




/* add duplicate and enable variables */
tVarList *AddDuplicateEnable(tVarList *vl, char *postfix) 
{
  tVarList *newvl;

  newvl = AddDuplicate(vl, postfix);
  enablevarlist(newvl);
  return newvl;
}




/* set: u = c */   
void vlsetconstant(tVarList *u, const double c)
{
  tGrid *grid = u->grid;
  double *pu;
  int i, n;
  int b;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int nnodes = box->nnodes;

    for (n = 0; n < u->n; n++)
    {
      pu = box->v[u->index[n]];

      #pragma omp parallel for
      for (i = 0; i < nnodes; i++)
        pu[i] = c;
    }
  }
}




/* copy: v = u */   
void vlcopy(tVarList *v, tVarList *u)
{
  tGrid *grid = v->grid;
  double *pu, *pv;
  int i, n;
  int b;

  /* copy time */
  v->time = u->time;

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    int nnodes = box->nnodes;

    for (n = 0; n < v->n; n++)
    {
      pu = box->v[u->index[n]];
      pv = box->v[v->index[n]];

      #pragma omp parallel for
      for (i = 0; i < nnodes; i++)
        pv[i] = pu[i];
    }
  }
}

void vlcopygrid(tGrid *grid, tVarList *v, tVarList *u)
{
  if (!grid || !v || !u) return;
  v->grid = u->grid = grid;
  vlcopy(v, u);
}

/* wrapper for single variable: v = u (iv/u is index of v/u) */
void varcopy(tGrid *grid, int iv, int iu)
{
  tVarList *v = vlalloc(grid);
  tVarList *u = vlalloc(grid);
  vlpushone(v, iv);
  vlpushone(u, iu);
  vlcopy(v, u);
  vlfree(u);
  vlfree(v);
}
            

/* average: r=(a+b)/2 */   
void vlaverage(tVarList *r, tVarList *a, tVarList *b)
{
  tGrid *grid = r->grid;
  double *pr, *pa, *pb;
  double c = 0.5;
  int i, n;
  int bi;

  for (bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;

    for (n = 0; n < r->n; n++)
    {
      pr = box->v[r->index[n]];
      pa = box->v[a->index[n]];
      pb = box->v[b->index[n]];

      #pragma omp parallel for
      for (i = 0; i < nnodes; i++)
        pr[i] = c * (pa[i] + pb[i]);
    }
  }
  /* average times as well */
  r->time = c * (a->time + b->time);
}




/* subtract two var lists: r = a - b
   can be called as vlsubtract(r,a,b); or vlsubtract(a,a,b); */   
void vlsubtract(tVarList *r, tVarList *a, tVarList *b)
{
  tGrid *grid = r->grid;
  double *pr, *pa, *pb;
  int i, n;
  int bi;

  for (bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;

    for (n = 0; n < r->n; n++)
    {
      pr = box->v[r->index[n]];
      pa = box->v[a->index[n]];
      pb = box->v[b->index[n]];

      #pragma omp parallel for
      for (i = 0; i < nnodes; i++)
        pr[i] = pa[i] - pb[i];
    }
  }
  /* subtract times as well */
  r->time = a->time - b->time;
}




/* linear combination of two var lists: r = ca*a + cb*b
   should change function name
   one function can catch several special cases like cb == 0 (unfinished)
   important: if coefficient is zero we guarantee that memory is not accessed
*/
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList *b) 
{
  tGrid *grid = r->grid;
  double *pr, *pa, *pb;
  int i, n;
  int bi;

  for (bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;

    for (n = 0; n < r->n; n++)
    {
      pr = box->v[r->index[n]];
      if(ca!=0)  pa = box->v[a->index[n]];
      if(cb!=0)  pb = box->v[b->index[n]];

      if (ca == 0 && cb == 0) {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] = 0; }

      else if (ca == 0) {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] = cb * pb[i]; }

      else if (cb == 0) {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] = ca * pa[i]; }

      else {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] = ca * pa[i] + cb * pb[i]; }
    }
  }
  /* add times as well */
  if (ca == 0 && cb == 0) r->time = 0.0;
  else if (ca == 0)	  r->time = cb * b->time;
  else if (cb == 0)	  r->time = ca * a->time;
  else			  r->time = ca * a->time + cb * b->time;
}

/* wrapper for single variable: r = ca*a + cb*b (ia/b/r is index of a/b/r) */
void varadd(tGrid *grid, int ir, double ca, int ia, double cb, int ib)
{
  tVarList *a = vlalloc(grid);
  tVarList *b = vlalloc(grid);
  tVarList *r = vlalloc(grid);
  vlpushone(a, ia);
  vlpushone(b, ib);
  vlpushone(r, ir);
  vladd(r, ca,a, cb,b);
  vlfree(a);
  vlfree(b);
  vlfree(r);
}

/* add second var list to first: r += ca*a
   special treatment for ca = 1 and ca = -1 */
void vladdto(tVarList *r, const double ca, tVarList *a) 
{
  tGrid *grid = r->grid;
  int bi;

  if (ca == 0) return;
  
  for (bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int nnodes = box->nnodes;
    double *pr, *pa;
    int i, n;

    for (n = 0; n < r->n; n++)
    {
      pr = box->v[r->index[n]];
      pa = box->v[a->index[n]];
  
      if (ca == 1) {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] += pa[i]; }
      else if (ca == -1) {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] -= pa[i]; }
      else {
        #pragma omp parallel for
        for (i = 0; i < nnodes; i++) pr[i] += ca * pa[i]; }
    }
  }
  /* add times as well */
  r->time += ca * a->time;
}
