/* storage.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "MemoryMan.h"

#define PR 0

int storage_verbose = 0;



/**************************************************************************/
/* basic memory management */

/* allocate grid */
tGrid *alloc_grid(int nboxes, int nvariables)
{
  int b;
  tGrid *grid;

  grid = (tGrid *) calloc( 1, sizeof(tGrid) );
  if (!grid) errorexit("alloc_grid:  out of memory");

  grid->nboxes = nboxes;
  grid->nvariables = nvariables;
  grid->iteration = 0;
  grid->time = 0;
  
  /* allocate mem for array of tBox * */
  grid->box = (tBox **) calloc( nboxes, sizeof(tBox *) );
  if (!grid->box) errorexit("alloc_grid:  out of memory 2");
  
  /* allocate mem for each box */
  for(b=0; b<nboxes ; b++)
  {
    grid->box[b] = (tBox *) calloc( 1, sizeof(tBox) );
    /* printf("b=%d grid->box[b]=%p\n", b, grid->box[b]); */
  }
  
  storage_verbose = Getv("storage_verbose", "yes");
  return grid;
}


/* allocate box */
tBox *alloc_box(tGrid *g, int b, int n1, int n2, int n3) 
{
  tBox *box = g->box[b];
  int nnodes = n1*n2*n3;
  int i,j,k;

  /* printf("b=%d box=%p g->box=%p\n", b, box, g->box); */

  /* allowed? */
  if (b >= g->nboxes || b<0 )
    errorexit("alloc_box: index out of range");

  box->grid = g;
  box->b = b;
  box->n1 = n1;
  box->n2 = n2;
  box->n3 = n3;
  box->nnodes = nnodes;
    
  /* allocate storage for nodes  */
  box->node = alloc_nodes(box, nnodes);

  /* allocate storage for data pointers, they default to NULL */
  box->v = NULL;
  realloc_boxvariables(box, g->nvariables);

  /* get mem. for diff. matrices */
  box->D1 = box->D2 = box->D3 = box->DD1 = box->DD2 = box->DD3 = NULL;
  box->D1 = (double *) calloc( n1*n1, sizeof(double) );
  box->D2 = (double *) calloc( n2*n2, sizeof(double) );
  box->D3 = (double *) calloc( n3*n3, sizeof(double) );
  box->DD1 = (double *) calloc( n1*n1, sizeof(double) );
  box->DD2 = (double *) calloc( n2*n2, sizeof(double) );
  box->DD3 = (double *) calloc( n3*n3, sizeof(double) );
  if( !(box->D1 && box->D2 && box->D3 && 
        box->DD1 && box->DD2 && box->DD3) )
    errorexit("alloc_box: out of memory for diff. matrices");

  /* let the coord trafo funcs point to NULL */
  for(i=0;i<4;i++)
  {
    box->x_of_X[i] = NULL;
    box->Sing_d_dx[i] = NULL;
    for(j=0;j<4;j++)
    {
      box->dX_dx[i][j] = NULL;
      box->dx_dX[i][j] = NULL;
      for(k=0;k<4;k++)
        box->ddX_dxdx[i][j][k] = NULL;
    }
  }

  /* get mem. for filter matrices */
  box->F1 = box->F2 = box->F3 = NULL;
  box->F1 = (double *) calloc( n1*n1, sizeof(double) );
  box->F2 = (double *) calloc( n2*n2, sizeof(double) );
  box->F3 = (double *) calloc( n3*n3, sizeof(double) );
  if( !(box->F1 && box->F2 && box->F3) )
    errorexit("alloc_box: out of memory for filter matrices");

  /* get mem. for matrices to obtain coeffs */
  box->Mcoeffs1 = box->Mcoeffs2 = box->Mcoeffs3 = NULL;
  box->Mcoeffs1 = (double *) calloc( n1*n1, sizeof(double) );
  box->Mcoeffs2 = (double *) calloc( n2*n2, sizeof(double) );
  box->Mcoeffs3 = (double *) calloc( n3*n3, sizeof(double) );
  if( !(box->Mcoeffs1 && box->Mcoeffs2 && box->Mcoeffs3) )
    errorexit("alloc_box: out of memory for matrices to obtain coeffs");

  /* get mem. for matrices for evaluation on points */
  box->Meval1 = box->Meval2 = box->Meval3 = NULL;
  box->Meval1 = (double *) calloc( n1*n1, sizeof(double) );
  box->Meval2 = (double *) calloc( n2*n2, sizeof(double) );
  box->Meval3 = (double *) calloc( n3*n3, sizeof(double) );
  if( !(box->Meval1 && box->Meval2 && box->Meval3) )
    errorexit("alloc_box: out of memory for matrices for evaluation on points");

  /* get mem. for matrices for integration */
  box->Int1 = box->Int2 = box->Int3 = NULL;
  box->Int1 = (double *) calloc( n1*n1, sizeof(double) );
  box->Int2 = (double *) calloc( n2*n2, sizeof(double) );
  box->Int3 = (double *) calloc( n3*n3, sizeof(double) );
  if( !(box->Int1 && box->Int2 && box->Int3) )
    errorexit("alloc_box: out of memory for matrices for integration");

  return box;
} 



/* allocate memory for several nodes */
tNode *alloc_nodes(tBox *box, int nnodes) 
{
  tNode *new;
  int i;

  /* get memory for node structure */
  new = (tNode *) calloc(nnodes, sizeof(tNode));
  if (!new) errorexit("alloc_nodes:  out of memory");

  /* fill in available data */
  for (i = 0; i < nnodes; i++) {
    new[i].i = i;
  }

  return new;
}  


/* add a bface to a box, f denotes the box face (0 to 5),
   return index in bface list */
int add_empty_bface(tBox *box, int f)
{
  int fi = box->nbfaces; /* add bface in this pos. in bface list */
  void *ret;

  /* increase size of bface list */
  ret = realloc( box->bface, (sizeof( *(box->bface) ))*(fi+1) );
  if(ret==NULL)  errorexit("add_bface: not enough memory for box->bface");
  box->bface = ret;

  /* mem for new bface */
  ret = calloc( 1, sizeof( *(box->bface[fi]) ) );
  if(ret==NULL) errorexit("add_bface: not enough memory for box->bface[n]");
  /* add new bface */
  box->bface[fi] = ret;
  box->nbfaces = fi+1;

  /* set some bface info */
  box->bface[fi]->grid = box->grid;
  box->bface[fi]->b    = box->b;
  box->bface[fi]->f    = f;
  box->bface[fi]->fi   = fi;
  box->bface[fi]->ob   = -1; /* other box is not known yet */
  return fi;
}

/* add a point ijk on face f to a bface with index fi, returns fi */
/* if called with fi<0, it first calls add_empty_bface and returns the new fi */
int add_point_to_bface_inbox(tBox *box, int fi, int ijk, int f)
{
  tBface *bface;
  /* make new bface when needed */
  if(fi<0 || fi>=box->nbfaces) fi = add_empty_bface(box, f);
  bface = box->bface[fi];
  /* make PointList id needed */
  if(bface->fpts==NULL)
    bface->fpts = AllocatePointList(box->grid);
  /* add point ijk */
  AddToPointList(bface->fpts, box->b, ijk);
  /* check if we have more than one face on this bface */
  if(bface->f!=f) bface->f = -1;
  return fi;
}

/* duplicate bface without pointlist fpts */
tBface *duplicate_bface_without_fpts_for_grid(tBface *bface0, tGrid *grid)
{
  tBface *bface;
  void *ret;
  if(bface0==NULL) return NULL;

  /* mem for new bface */
  ret = calloc( 1, sizeof( *(bface) ) );
  if(ret==NULL)
    errorexit("duplicate_bface_without_fpts_for_grid: not enough memory");
  bface = ret;

  /* make a shallow copy of the struct */
  *bface = *bface0;
  /* now set grid pointer */
  bface->grid = grid;
  /* remove pointer to bface0->fpts */
  bface->fpts = NULL;
  
  return bface;
}

/* duplicate bface with pointlist fpts */
tBface *duplicate_bface_for_grid(tBface *bface0, tGrid *grid)
{
  tBface *bface = duplicate_bface_without_fpts_for_grid(bface0, grid);
  /* copy fpts */
  bface->fpts = DuplicatePointList_for_grid(bface0->fpts, grid);
  return bface;
}

/* free a bface */
void free_bface(tBface *bface)
{
  if(bface!=NULL)
  {
    /* free the point lists of the faces */ 
    FreePointList(bface->fpts);
    free(bface);
  }
}

/* remove a bface with index fi from a box, return number of bfaces removed */
int remove_bface(tBox *box, int fi)
{
  void *ret;
  int nbfaces = box->nbfaces;
  int i;

  /* return 0 if bface does not exist */
  if(fi<0 || fi>=nbfaces) return 0;

  /* free the bface */
  free_bface(box->bface[fi]);

  /* shift bfaces behind fi one position to the front */
  for(i=fi; i<nbfaces-1; i++)
  {
    box->bface[i] = box->bface[i+1];
    box->bface[i]->fi = i;
  }

  /* reduce size of bface list */
  nbfaces--;
  ret = realloc( box->bface, (sizeof( *(box->bface) ))*(nbfaces) );
  if(ret==NULL && nbfaces!=0)
    errorexit("remove_bface: not enough memory for box->bface");
  box->bface = ret;
  box->nbfaces = nbfaces;

  return 1;
}

/* look for empty bfaces and remove them */
int remove_bfaces_with_NULL_ftps(tBox *box)
{
  int n, fi;
  n=0;
  for(fi=0; fi<box->nbfaces; fi++)
    if(box->bface[fi]->fpts == NULL)
    {
      int r = remove_bface(box, fi); /* this decreases box->nbfaces */
      fi=fi-r;  /* go back by one in fi so that for-loop covers fi again */
      n++;      /* count number of bfaces removed */
    }
  return n;
}

/* free box, currently leaves grid untouched */
void free_box(tBox *box) 
{
  int i;

  if (!box) return;

  for (i = 0; i < box->grid->nvariables; i++)
    disablevarcomp_inbox(box, i);
  free(box->v);
  free(box->node);
  free(box->D1);
  free(box->D2);
  free(box->D3);
  free(box->DD1);
  free(box->DD2);
  free(box->DD3);
  free(box->F1);
  free(box->F2);
  free(box->F3);
  free(box->Mcoeffs1);
  free(box->Mcoeffs2);
  free(box->Mcoeffs3);
  free(box->Meval1);
  free(box->Meval2);
  free(box->Meval3);
  free(box->Int1);
  free(box->Int2);
  free(box->Int3);

  for(i=0; i < box->nbfaces; i++)  free_bface(box->bface[i]);
  free(box->bface);

  free(box);
}


/* free grid without freeing boxes */
void free_grid_only(tGrid *g)
{
  if (!g) return;
  if (g->box) free(g->box);
  free(g);
}


/* free grid */
void free_grid(tGrid *g)
{
  int i;

  if (!g) return;
  for (i = 0; i < g->nboxes; i++)
    free_box(g->box[i]);
  free_grid_only(g);
}




/* several utility functions to enable and disable storage for variables */
int nenabled = 0;


/* enable one component of variable on one box */
void enablevarcomp_inbox(tBox *box, int i)
{
  if (!box->v[i]) {
    box->v[i] = calloc( box->nnodes, sizeof(double) );
    if (!box->v[i]) {
      printf("enable storage for variable %d = %s: out of memory\n",
	     i, VarName(i));
      errorexit("");
    }
    nenabled++;
    if (storage_verbose) {
      printf("enabling variable %d = %s in b=%d, nenabled=%d => ",
	     i, VarName(i), box->b, nenabled);
      printf("%d*%d = %d bytes\n", 
	     box->nnodes, (int) sizeof(double), 
	     box->nnodes * ((int) (sizeof(double))) );
    } 
  }
}


/* disable one component of variable */
void disablevarcomp_inbox(tBox *box, int i) 
{
  if (box->v[i]) {
    free(box->v[i]);
    box->v[i] = 0;
    nenabled--;
    if (storage_verbose) {
      printf("disabling variable %d = %s in b=%d, nenabled=%d => ",
             i, VarName(i), box->b, nenabled);
      printf("%d*%d = %d bytes\n", 
	     box->nnodes, (int) sizeof(double), 
	     box->nnodes * ((int) (sizeof(double))) );
    } 
  }
}


/* enable all components of a variable on one box */
void enablevar_inbox(tBox *box, int i)
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    enablevarcomp_inbox(box, i+j);
}


/* disable all components of a variable on one box */
void disablevar_inbox(tBox *box, int i)
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    disablevarcomp_inbox(box, i+j);
}




/* enable one component of variable on whole grid */
void enablevarcomp(tGrid *grid, int i)
{
  int b;
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];

    if(storage_verbose)
      printf("enablevarcomp: box=%p, box->b=%d:\n", box, box->b);
    enablevarcomp_inbox(box, i);
  }
}



/* disable one component of variable  on whole grid */
void disablevarcomp(tGrid *grid, int i) 
{
  int b;
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];

    if(storage_verbose)
      printf("disablevarcomp: box=%p, box->b=%d:\n ", box, box->b);
    disablevarcomp_inbox(box, i);
  }
}




/* enable all components of a variable */
void enablevar(tGrid *grid, int i) 
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    enablevarcomp(grid, i+j);
}




/* disable all components of a variable */
void disablevar(tGrid *grid, int i) 
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    disablevarcomp(grid, i+j);
}




/* enable variable list */
void enablevarlist(tVarList *v) 
{
  int i;

  if (v) 
    for (i = 0; i < v->n; i++)
      enablevarcomp(v->grid, v->index[i]);
}




/* disable variable list */
void disablevarlist(tVarList *v) 
{
  int i;

  if (v)
    for (i = 0; i < v->n; i++)
      disablevarcomp(v->grid, v->index[i]);
}




/* enable all variables found on given grid */
void enablesamevars(tGrid *grid, tGrid *newgrid)
{
  int i, b;

  if (grid->nvariables != newgrid->nvariables)
  {
    printf("nvariables: old %d, new %d\n", 
	   grid->nvariables, newgrid->nvariables);
    errorexit("enablesamevars: need same number of variables");
  }
  if (grid->nboxes != newgrid->nboxes)
  {
    printf("nboxes: old %d, new %d\n", 
	   grid->nboxes, newgrid->nboxes);
    errorexit("enablesamevars: need same number of boxes");
  }

  for (i = 0; i < newgrid->nvariables; i++)
    for (b = 0; b < newgrid->nboxes; b++)
    {
      tBox *newbox = newgrid->box[b];
      if(grid->box[b]->v[i])
        enablevarcomp_inbox(newbox, i);
    }
}



/* create room for more variables in one box */
void realloc_boxvariables(tBox *box, int nvariables)
{
  tGrid *grid = box->grid;
  int i;
  int box_nvariables = grid->nvariables;

  if (PR) printf(" realloc_boxvariables in box%d from %d to %d\n", 
		 box->b, grid->nvariables, nvariables);

  if(box->v == NULL) box_nvariables=0;
  box->v = (double **) realloc(box->v, sizeof(double *) * nvariables);
  if(nvariables > box_nvariables) 
    for(i=box_nvariables; i<nvariables; i++)  box->v[i] = 0;
}

/* create room for more variables */
void realloc_gridvariables(tGrid *grid, int nvariables)
{
  int b;

  if (PR) printf("realloc_gridvariables from %d to %d\n", 
		 grid->nvariables, nvariables);

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    realloc_boxvariables(box, nvariables);
  }
  grid->nvariables = nvariables;
}

