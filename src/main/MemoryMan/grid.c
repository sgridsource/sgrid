/* grid.c */
/* Wolfgang Tichy, April 2005 */


#include "sgrid.h"
#include "MemoryMan.h"

#define PR 0




/* initialize grid 
   called in main()
   here we sort out the various parameter options before calling
   for example make_grid_box
*/
tGrid *make_grid(int pr)
{
  tGrid *g;
  tBox *box;
  int b;
  int n1, n2, n3;
  int nboxes = Geti("nboxes");

  /* print info */
  if (pr) {
    prdivider(0);
    printf("Creating grid\n");
  }

  /* make grid structure with nboxes=1, nvariables */
  g = alloc_grid(nboxes, globalnvariables);
  printf("g->nboxes=%d g->box=%p\n", g->nboxes, g->box);

  /* make box structure with so many nodes and add it to g */
  for(b=0; b<nboxes ; b++)
  {  
    char str[1000];
  
    snprintf(str, 999, "box%d_n1", b);
    n1=Geti(str);
    snprintf(str, 999, "box%d_n2", b);
    n2=Geti(str);
    snprintf(str, 999, "box%d_n3", b);
    n3=Geti(str);

    box = alloc_box(g, b, n1, n2, n3);
  }
  
  /* time step */
  g->dt = Getd("dt");

  /* set all pointers in each tBox struct of grid from Pars */
  set_BoxStructures_fromPars(g, 0);
  if(pr) printgrid(g);

  /* return pointer to newly created grid */
  return g;
}


/* set all pointers (except var contents) in each tBox struct 
   of grid from Pars */
int set_BoxStructures_fromPars(tGrid *g, int pr)
{
  tBox *box;
  tNode *node;
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int b, i, j, k, ijk;
  int n1, n2, n3;
  int nboxes = Geti("nboxes");
  int basis;

  /* print info */
  if (pr) {
    prdivider(0);
    printf("set_BoxStructures_fromPars: setting box structures on grid=%p\n",
           g);
  }

  /* enable storage of the coordinate variables */
  enablevar(g, var_X);
  enablevar(g, var_Y);
  enablevar(g, var_Z);

  /* fill in connectivity and coordinates of nodes */
  for(b=0; b<nboxes ; b++)
  {
    int dir;
    char str[1000];
  
    box = g->box[b];

    /* get values of bounding box */
    snprintf(str, 999, "box%d_min1", b);
    box->bbox[0] = Getd(str);
    snprintf(str, 999, "box%d_max1", b);
    box->bbox[1] = Getd(str);

    snprintf(str, 999, "box%d_min2", b);
    box->bbox[2] = Getd(str);
    snprintf(str, 999, "box%d_max2", b);
    box->bbox[3] = Getd(str);

    snprintf(str, 999, "box%d_min3", b);
    box->bbox[4] = Getd(str);
    snprintf(str, 999, "box%d_max3", b);
    box->bbox[5] = Getd(str);

    /* check for periodic dirs */
    for(dir=1; dir<4; dir++)
    {
      snprintf(str, 999, "box%d_basis%d", b, dir);
      if( Getv(str, "Fourier") || Getv(str, "fd2_periodic") )
        box->periodic[dir]=1;
      else
        box->periodic[dir]=0;
    }

    n1 = box->n1;
    n2 = box->n2;
    n3 = box->n3;

    node = box->node;
    ijk = 0;
    for (k = 0; k < n3; k++)
    { 
      for (j = 0; j < n2; j++)
      {  
        for (i = 0; i < n1; i++, ijk++)
        {
	node[ijk].i = ijk;

	/* node type */	
	node[ijk].type = 1;

	/* compute coordinate values */
        for(basis=1; basis<=3 ; basis++)
        {
          int nb, bbi, varb, m;

          nb=n1; bbi=0; varb=var_X; m=i; /* <--for basis=1 */
          if(basis==2)
          {
            nb=n2; bbi=2; varb=var_Y; m=j;
          }
          else if(basis==3)
          {
            nb=n3; bbi=4; varb=var_Z; m=k;
          }

          /* compute coordinate values */
	  snprintf(str, 999, "box%d_basis%d", b, basis);
	  if( Getv(str, "ChebExtrema") )
	  {
	    if(m==0)          box->v[varb][ijk] = box->bbox[bbi];
	    else if(m==nb-1)  box->v[varb][ijk] = box->bbox[bbi+1];
	    else
	      box->v[varb][ijk] 
	        = 0.5*( (box->bbox[bbi] - box->bbox[bbi+1])*cos(m*PI/(nb-1)) 
	               +(box->bbox[bbi] + box->bbox[bbi+1]));
	  }
	  else if( Getv(str, "ChebZeros") )
	    box->v[varb][ijk] 
	      = 0.5*( (box->bbox[bbi] - box->bbox[bbi+1])
	              *cos( ((double) m + 0.5)*PI/nb) 
	             +(box->bbox[bbi] + box->bbox[bbi+1]));
          else if( Getv(str, "Fourier") || Getv(str, "fd2_periodic") )
	    box->v[varb][ijk] 
	      = ( (box->bbox[bbi+1] - box->bbox[bbi])* ((double) m)/nb 
	             +box->bbox[bbi]);
          else if( Getv(str, "fd2_onesided") )
	    box->v[varb][ijk] 
	      = ( (box->bbox[bbi+1] - box->bbox[bbi])* ((double) m)/(nb-1)
	             +box->bbox[bbi]);
        }/* end for basis */
        }
      }
    }  /* end for (k = 0; k < n3; k++) */

    /* set flags which determine if we use FFTs */
    set_TransformType_flags_inbox(box);

    /* initialize funcs in box struct */
    init_spec_functionpointers(box);

    /* initialize diff. and filter matrices */
    {
       void (*get_coeffs)(double *,double *, int)=NULL;
       void (*coeffs_of_deriv)(void *, double, double, double *,double *, int)=NULL;
       void (*coeffs_of_2ndderiv)(void *,double, double, double *,double *, int)=NULL;
       void (*coeffs_of_int)(double, double, double *,double *, int)=NULL;
       void (*eval_onPoints)(double *,double *, int)=NULL;
       void (*filter_coeffs)(double *, int, int)=NULL;
       double (*basisfunc)(void *aux, double a, double b, int k, int N, double X)=NULL;
       int filt1, filt2, filt3;

       snprintf(str, 999, "box%d_filter1", b);
       filt1 = Geti(str);

       snprintf(str, 999, "box%d_filter2", b);
       filt2 = Geti(str);

       snprintf(str, 999, "box%d_filter3", b);
       filt3 = Geti(str);

       /* direction 1 */
       get_spec_functionpointers(box, 1, &get_coeffs, &coeffs_of_deriv,
                                 &coeffs_of_2ndderiv, &coeffs_of_int,
                                 &eval_onPoints, &filter_coeffs, &basisfunc);
       initdiffmatrix(box, 1, box->D1, box->DD1, n1,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F1, n1+1-filt1, n1, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box, 1, box->DD1, n1,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       initMatrix_ForCoeffs(box->Mcoeffs1, n1, get_coeffs);
       initMatrix_ToEvaluate(box->Meval1,  n1, eval_onPoints);
       initIntegrationMatrix(box->bbox[0], box->bbox[1], box->Int1, n1,
                             get_coeffs, coeffs_of_int, eval_onPoints);

       /* direction 2 */
       get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                                 &coeffs_of_2ndderiv, &coeffs_of_int,
                                 &eval_onPoints, &filter_coeffs, &basisfunc);
       initdiffmatrix(box, 2, box->D2, box->DD2, n2,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F2, n2+1-filt2, n2, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box, 2, box->DD2, n2,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       initMatrix_ForCoeffs(box->Mcoeffs2, n2, get_coeffs);
       initMatrix_ToEvaluate(box->Meval2,  n2, eval_onPoints);
       initIntegrationMatrix(box->bbox[2], box->bbox[3], box->Int2, n2,
                             get_coeffs, coeffs_of_int, eval_onPoints);

       /* direction 3 */
       get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                                 &coeffs_of_2ndderiv, &coeffs_of_int,
                                 &eval_onPoints, &filter_coeffs, &basisfunc);
       initdiffmatrix(box, 3, box->D3, box->DD3, n3,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F3, n3+1-filt3, n3, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box, 3, box->DD3, n3,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       initMatrix_ForCoeffs(box->Mcoeffs3, n3, get_coeffs);
       initMatrix_ToEvaluate(box->Meval3,  n3, eval_onPoints);
       initIntegrationMatrix(box->bbox[4], box->bbox[5], box->Int3, n3,
                             get_coeffs, coeffs_of_int, eval_onPoints);
    }
  } /* end box loop */

  /* set coord. info CI inside each box from box_CI pars */
  set_box_CI_struct_from_pars(g);

  if(pr) printgrid(g);
  return 0;
}


/* make an empty grid, into which we can e.g. copy the contents of
   an existing grid */
tGrid *make_empty_grid(int nvariables, int pr)
{
  tGrid *g;
  tBox *box;
  int b;
  int n1, n2, n3;
  int nboxes = Geti("nboxes");

  /* print info */
  if(pr) 
  {
    prdivider(0);
    printf("make_empty_grid:  ");
  }

  /* make grid structure with nboxes, nvariables */
  g = alloc_grid(nboxes, nvariables);
  if(pr) printf("g->nboxes=%d  g->box=%p  nvariables=%d\n",
                g->nboxes, g->box, nvariables);

  /* make box structure with so many nodes and add it to g */
  for(b=0; b<nboxes ; b++)
  {  
    char str[1000];
  
    snprintf(str, 999, "box%d_n1", b);
    n1=Geti(str);
    snprintf(str, 999, "box%d_n2", b);
    n2=Geti(str);
    snprintf(str, 999, "box%d_n3", b);
    n3=Geti(str);

    box = alloc_box(g, b, n1, n2, n3);
  }
  
  /* time step */
  g->dt = Getd("dt");

  if(pr) printgrid(g);

  /* return pointer to newly created grid */
  return g;
}


/* copy all the contents, except the vars, from g_old into g_new */
int copy_grid_withoutvars(tGrid *g_old, tGrid *g_new, int pr)
{
  int b, i, j, k, ijk;
  int n1, n2, n3;

  /* print info */
  if(pr)
  {
    prdivider(0);
    printf("copy_grid_withoutvars:\n");
  }

  /* copy struct tGRID */
  g_new->nboxes     = g_old->nboxes;
  g_new->ninactive  = g_old->ninactive;
  g_new->nvariables = g_old->nvariables;
  g_new->iteration  = g_old->iteration;
  g_new->time	    = g_old->time;
  g_new->dt         = g_old->dt;

  /* copy box contents */
  for(b=0; b<g_new->nboxes ; b++)
  {
    tBox *b_new = g_new->box[b];
    tBox *b_old = g_old->box[b];

    /* copy values of bounding boxes */
    for(i=0; i<6; i++)
    {
      b_new->bbox[i]  = b_old->bbox[i];
      b_new->ibbox[i] = b_old->ibbox[i];
    }

    n1 = b_new->n1;
    n2 = b_new->n2;
    n3 = b_new->n3;

    /* copy nodes */
    ijk = 0;
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++, ijk++)
        {
          /* copy nodes */
          b_new->node[ijk].i    = b_old->node[ijk].i;	
          b_new->node[ijk].type = b_old->node[ijk].type;
        }

    /* copy function pointers */
    for (i = 0; i < 4; i++, ijk++)
    {
      b_new->x_of_X[i]    = b_old->x_of_X[i];
      b_new->XYZ_Of_xyz   = b_old->XYZ_Of_xyz;
      b_new->Sing_d_dx[i] = b_old->Sing_d_dx[i];
      for (j = 0; j < 4; j++)
      {
        b_new->dX_dx[i][j] = b_old->dX_dx[i][j];
        b_new->dx_dX[i][j] = b_old->dx_dX[i][j];
        for (k = 0; k < 4; k++)
          b_new->ddX_dxdx[i][j][k] = b_old->ddX_dxdx[i][j][k];
      }
    }
    b_new->isSing = b_old->isSing;
    b_new->basis1 = b_old->basis1;
    b_new->basis2 = b_old->basis2;
    b_new->basis3 = b_old->basis3;

    /* copy FT func pointers */
    b_new->get_coeffs1 = b_old->get_coeffs1;
    b_new->get_coeffs2 = b_old->get_coeffs2;
    b_new->get_coeffs3 = b_old->get_coeffs3;
    b_new->coeffs_of_deriv1 = b_old->coeffs_of_deriv1;
    b_new->coeffs_of_deriv2 = b_old->coeffs_of_deriv2;
    b_new->coeffs_of_deriv3 = b_old->coeffs_of_deriv3;
    b_new->coeffs_of_2ndderiv1 = b_old->coeffs_of_2ndderiv1;
    b_new->coeffs_of_2ndderiv2 = b_old->coeffs_of_2ndderiv2;
    b_new->coeffs_of_2ndderiv3 = b_old->coeffs_of_2ndderiv3;
    b_new->coeffs_of_int1 = b_old->coeffs_of_int1;
    b_new->coeffs_of_int2 = b_old->coeffs_of_int2;
    b_new->coeffs_of_int3 = b_old->coeffs_of_int3;
    b_new->eval_onPoints1 = b_old->eval_onPoints1;
    b_new->eval_onPoints2 = b_old->eval_onPoints2;
    b_new->eval_onPoints3 = b_old->eval_onPoints3;
    b_new->filter_coeffs1 = b_old->filter_coeffs1;
    b_new->filter_coeffs2 = b_old->filter_coeffs2;
    b_new->filter_coeffs3 = b_old->filter_coeffs3;
    /* copy type of FT */
    b_new->TransformType1 = b_old->TransformType1;
    b_new->TransformType2 = b_old->TransformType2;
    b_new->TransformType3 = b_old->TransformType3;
    /* copy periodic info */
    for(i = 0; i < 4; i++)
       b_new->periodic[i] = b_old->periodic[i];

    /* copy attributes of boxes */
    b_new->inactive = b_old->inactive;
    for(i = 0; i < NATTRIBS; i++)
      b_new->Attrib[i] = b_old->Attrib[i];

    /* copy bfaces: free all in b_new and then duplicate b_old */
    free_all_bfaces(b_new);
    b_new->nbfaces = b_old->nbfaces;
    b_new->bface = calloc(b_new->nbfaces, sizeof( *(b_new->bface) ));
    if(b_new->bface==NULL)
      errorexit("copy_grid_withoutvars: not enough memory for box->bface");
    for(i = 0; i < b_old->nbfaces; i++)
      b_new->bface[i] = duplicate_bface_for_grid(b_old->bface[i], g_new);

    /* copy diff., filter matrices, ..., and all other arrays */
    for (i = 0; i < n1*n1; i++, ijk++)
    {
      b_new->D1[i]  = b_old->D1[i];
      b_new->DD1[i] = b_old->DD1[i];
      b_new->F1[i]  = b_old->F1[i];
      b_new->Mcoeffs1[i] = b_old->Mcoeffs1[i];
      b_new->Meval1[i]   = b_old->Meval1[i];
      b_new->Int1[i]     = b_old->Int1[i];
    }
    for (j = 0; j < n2*n2; j++)
    {
      b_new->D2[j]  = b_old->D2[j];
      b_new->DD2[j] = b_old->DD2[j];
      b_new->F2[j]  = b_old->F2[j];
      b_new->Mcoeffs2[j] = b_old->Mcoeffs2[j];
      b_new->Meval2[j]   = b_old->Meval2[j];
      b_new->Int2[j]     = b_old->Int2[j];
    }
    for (k = 0; k < n3*n3; k++)
    {
      b_new->D3[k]  = b_old->D3[k];
      b_new->DD3[k] = b_old->DD3[k];
      b_new->F3[k]  = b_old->F3[k];
      b_new->Mcoeffs3[k] = b_old->Mcoeffs3[k];
      b_new->Meval3[k]   = b_old->Meval3[k];
      b_new->Int3[k]     = b_old->Int3[k];
    }

    /* shallow copy of coord. info in CI */
    b_new->CI[0] = b_old->CI[0];
  }
  if(pr) printgrid(g_new);

  return 1;
}


/* copy all the contents, including the vars, from g_old into g_new */
int copy_grid(tGrid *g_old, tGrid *g_new, int pr)
{
  int b, ind, ijk;

  /* print info */
  if(pr)
  {
    prdivider(0);
    printf("copy_grid:\n");
  }
  
  /* copy all the contents, except the vars, from g_old into g_new */
  copy_grid_withoutvars(g_old, g_new, 0);

  /* enable all vars needed io g_old */
  enablesamevars(g_old, g_new);

  /* copy all vars in all boxes */
  forallboxes(g_new, b)
  {
    tBox *box_new = g_new->box[b];
    tBox *box_old = g_old->box[b];

    for(ind=0; ind<g_new->nvariables; ind++)
      if(box_new->v[ind] && box_old->v[ind])
        for(ijk=0; ijk < min2(box_new->nnodes, box_old->nnodes); ijk++)
          box_new->v[ind][ijk] = box_old->v[ind][ijk];
  }
  if(pr) printgrid(g_new);

  return 1;
}

/* copy var with index vi from g_old into g_new */
int copy_gridvar(int vi, tGrid *g_old, tGrid *g_new)
{
  int b, ijk;

  /* copy var with index vi in all boxes */
  forallboxes(g_new, b)
  {
    tBox *box_new = g_new->box[b];
    tBox *box_old = g_old->box[b];

    if(box_new->v[vi] && box_old->v[vi])
      for(ijk=0; ijk < min2(box_new->nnodes, box_old->nnodes); ijk++)
        box_new->v[vi][ijk] = box_old->v[vi][ijk];
  }
  return 1;
}


/* let var pointers in g_new point to the same memory as in g_old */
int point_grid_tosamevars(tGrid *g_old, tGrid *g_new, int pr)
{
  int b, ind;

  /* print info */
  if(pr)
  {
    prdivider(0);
    printf("point_grid_tosamevars:\n");
  }
  
  /* copy all var pointers in all boxes */
  forallboxes(g_new, b)
  {
    tBox *box_new = g_new->box[b];
    tBox *box_old = g_old->box[b];

    for(ind=0; ind<g_new->nvariables; ind++)
        box_new->v[ind] = box_old->v[ind];
  }
  if(pr) printgrid(g_new);

  return 1;
}


/* let var pointers in g_new point to NULL */
/* CAUTION: If we use free_grid after set_gridvars_toNULL the vars 
            do not get freed, and if their pointers are not saved in 
            another grid we get a memory leak! */
int set_gridvars_toNULL(tGrid *g_new, int pr)
{
  int b, ind;

  /* print info */
  if(pr)
  {
    prdivider(0);
    printf("set_gridvars_toNULL:\n");
  }

  /* set var pointers in all boxes to NULL */
  forallboxes(g_new, b)
  {
    tBox *box_new = g_new->box[b];

    for(ind=0; ind<g_new->nvariables; ind++)
        box_new->v[ind] = NULL;
  }
  if(pr) printgrid(g_new);

  return 1;
}


/**********************************************/
/* functions to deactivate and activate boxes */
/**********************************************/

/* set grid->ninactive to number of boxes with with box->inactive==0 */
void set_ninactive(tGrid *grid)
{
  int b, n=0;
  forallboxes(grid,b) n++;
  grid->ninactive = grid->nboxes - n;
}

/* set box->inactive = fl for a list */
void set_box_inactiveflags(tGrid *grid, intList *ilist, int fl)
{
  int i;
  for(i=0; i<ilist->n; i++) grid->box[ilist->e[i]]->inactive = fl;
  set_ninactive(grid);
}
/* deactivate boxes in ilist */
void deactivate_boxes(tGrid *grid, intList *ilist)
{
  set_box_inactiveflags(grid, ilist, 1);
}
/* activate boxes in ilist */
void activate_boxes(tGrid *grid, intList *ilist)
{
  set_box_inactiveflags(grid, ilist, 0);
}
/* deactivate all boxes */
void deactivate_allboxes(tGrid *grid)
{
  int b;
  for(b=0; b<grid->nboxes; b++) grid->box[b]->inactive = 1;
  set_ninactive(grid);
}
/* activate all boxes */
void activate_allboxes(tGrid *grid)
{
  int b;
  for(b=0; b<grid->nboxes; b++) grid->box[b]->inactive = 0;
  set_ninactive(grid);
}

/* find index b of ith box that has box->inactive==0
   e.g. if box[0]->inactive=1 and box[1]->inactive=0, 
   IndexOf_ith_NonSkipBox(grid, 0) returns 1  */
int IndexOf_ith_ActiveBox(tGrid *grid, int i)
{
  int b, nn=0;

  if(grid->ninactive == 0) 
  {
    if(i>=0 && i < grid->nboxes) return i;
    else                         return -1;
  }
  forallboxes(grid,b)
  {
    if(nn==i) break;
    nn++;
  }
  if(b < grid->nboxes) return b;
  else                 return -1;
}
