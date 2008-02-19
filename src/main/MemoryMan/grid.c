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

  /* provide a hook for special grid preparation */
  RunFun(PRE_GRID, 0);

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

  /* provide a hook for special treatmet after grid creation */
  RunFun(POST_GRID, g);

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
	    box->v[varb][ijk] 
	      = 0.5*( (box->bbox[bbi] - box->bbox[bbi+1])*cos(m*PI/(nb-1)) 
	             +(box->bbox[bbi] + box->bbox[bbi+1]));
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

    /* initialize diff. and filter matrices */
    {
       void (*get_coeffs)(double *,double *, int)=NULL;
       void (*coeffs_of_deriv)(double, double, double *,double *, int)=NULL;
       void (*coeffs_of_2ndderiv)(double, double, double *,double *, int)=NULL;
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
       initdiffmatrix(box->bbox[0], box->bbox[1], box->D1, box->DD1, n1,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F1, n1+1-filt1, n1, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box->bbox[0], box->bbox[1], box->DD1, n1,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       box->basis1=basisfunc;
       initMatrix_ForCoeffs(box->Mcoeffs1, n1, get_coeffs);
       initMatrix_ToEvaluate(box->Meval1,  n1, eval_onPoints);
       initIntegrationMatrix(box->bbox[0], box->bbox[1], box->Int1, n1,
                             get_coeffs, coeffs_of_int, eval_onPoints);

       /* direction 2 */
       get_spec_functionpointers(box, 2, &get_coeffs, &coeffs_of_deriv,
                                 &coeffs_of_2ndderiv, &coeffs_of_int,
                                 &eval_onPoints, &filter_coeffs, &basisfunc);
       initdiffmatrix(box->bbox[2], box->bbox[3], box->D2, box->DD2, n2,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F2, n2+1-filt2, n2, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box->bbox[2], box->bbox[3], box->DD2, n2,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       box->basis2=basisfunc;
       initMatrix_ForCoeffs(box->Mcoeffs2, n2, get_coeffs);
       initMatrix_ToEvaluate(box->Meval2,  n2, eval_onPoints);
       initIntegrationMatrix(box->bbox[2], box->bbox[3], box->Int2, n2,
                             get_coeffs, coeffs_of_int, eval_onPoints);

       /* direction 3 */
       get_spec_functionpointers(box, 3, &get_coeffs, &coeffs_of_deriv,
                                 &coeffs_of_2ndderiv, &coeffs_of_int,
                                 &eval_onPoints, &filter_coeffs, &basisfunc);
       initdiffmatrix(box->bbox[4], box->bbox[5], box->D3, box->DD3, n3,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F3, n3+1-filt3, n3, 
                        get_coeffs, filter_coeffs, eval_onPoints);
       if(coeffs_of_2ndderiv!=NULL)
         initdiffmatrix2(box->bbox[4], box->bbox[5], box->DD3, n3,
                         get_coeffs, coeffs_of_2ndderiv, eval_onPoints);
       box->basis3=basisfunc;
       initMatrix_ForCoeffs(box->Mcoeffs3, n3, get_coeffs);
       initMatrix_ToEvaluate(box->Meval3,  n3, eval_onPoints);
       initIntegrationMatrix(box->bbox[4], box->bbox[5], box->Int3, n3,
                             get_coeffs, coeffs_of_int, eval_onPoints);
    }
  } /* end box loop */
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
  g_new->nvariables = g_old->nvariables;
  g_new->iteration  = g_old->iteration;
  g_new->time	    = g_old->time;
  g_new->dt         = g_old->dt;

  /* copy box contents */
  for(b=0; b<g_new->nboxes ; b++)
  {
    /* copy values of bounding boxes */
    for(i=0; i<6; i++)
    {
      g_new->box[b]->bbox[i]  = g_old->box[b]->bbox[i];
      g_new->box[b]->ibbox[i] = g_old->box[b]->ibbox[i];
    }

    n1 = g_new->box[b]->n1;
    n2 = g_new->box[b]->n2;
    n3 = g_new->box[b]->n3;

    /* copy nodes */
    ijk = 0;
    for (k = 0; k < n3; k++)
      for (j = 0; j < n2; j++)
        for (i = 0; i < n1; i++, ijk++)
        {
          /* copy nodes */
          g_new->box[b]->node[ijk].i    = g_old->box[b]->node[ijk].i;	
          g_new->box[b]->node[ijk].type = g_old->box[b]->node[ijk].type;
        }

    /* copy function pointers */
    for (i = 0; i < 4; i++, ijk++)
    {
      g_new->box[b]->x_of_X[i]    = g_old->box[b]->x_of_X[i];
      g_new->box[b]->Sing_d_dx[i] = g_old->box[b]->Sing_d_dx[i];
      for (j = 0; j < 4; j++)
      {
        g_new->box[b]->dX_dx[i][j] = g_old->box[b]->dX_dx[i][j];
        for (k = 0; k < 4; k++)
          g_new->box[b]->ddX_dxdx[i][j][k] = g_old->box[b]->ddX_dxdx[i][j][k];
      }
    }
    g_new->box[b]->basis1 = g_old->box[b]->basis1;
    g_new->box[b]->basis2 = g_old->box[b]->basis2;
    g_new->box[b]->basis3 = g_old->box[b]->basis3;

    /* copy diff., filter matrices, ..., and all other arrays */
    for (i = 0; i < n1*n1; i++, ijk++)
    {
      g_new->box[b]->D1[i]  = g_old->box[b]->D1[i];
      g_new->box[b]->DD1[i] = g_old->box[b]->DD1[i];
      g_new->box[b]->F1[i]  = g_old->box[b]->F1[i];
      g_new->box[b]->Mcoeffs1[i] = g_old->box[b]->Mcoeffs1[i];
      g_new->box[b]->Meval1[i]   = g_old->box[b]->Meval1[i];
      g_new->box[b]->Int1[i]     = g_old->box[b]->Int1[i];
    }
    for (j = 0; j < n2*n2; j++)
    {
      g_new->box[b]->D2[j]  = g_old->box[b]->D2[j];
      g_new->box[b]->DD2[j] = g_old->box[b]->DD2[j];
      g_new->box[b]->F2[j]  = g_old->box[b]->F2[j];
      g_new->box[b]->Mcoeffs2[j] = g_old->box[b]->Mcoeffs2[j];
      g_new->box[b]->Meval2[j]   = g_old->box[b]->Meval2[j];
      g_new->box[b]->Int2[j]     = g_old->box[b]->Int2[j];
    }
    for (k = 0; k < n3*n3; k++)
    {
      g_new->box[b]->D3[k]  = g_old->box[b]->D3[k];
      g_new->box[b]->DD3[k] = g_old->box[b]->DD3[k];
      g_new->box[b]->F3[k]  = g_old->box[b]->F3[k];
      g_new->box[b]->Mcoeffs3[k] = g_old->box[b]->Mcoeffs3[k];
      g_new->box[b]->Meval3[k]   = g_old->box[b]->Meval3[k];
      g_new->box[b]->Int3[k]     = g_old->box[b]->Int3[k];
    }
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
      for(ijk=0; ijk < min2(box_new->nnodes, box_old->nnodes); ijk++)
        if(box_new->v[ind])
          box_new->v[ind][ijk] = box_old->v[ind][ijk];
  }
  if(pr) printgrid(g_new);

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
