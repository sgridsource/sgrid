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
  tNode *node;
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int b, i, j, k, ijk;
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
	snprintf(str, 999, "box%d_basis1", b);
	if( Getv(str, "ChebExtrema") )
	  box->v[var_X][ijk] 
	    = 0.5*( (box->bbox[0] - box->bbox[1])*cos(i*PI/(n1-1)) 
	           +(box->bbox[0] + box->bbox[1]));
        else if( Getv(str, "Fourier") )
	  box->v[var_X][ijk] 
	    = ( (box->bbox[1] - box->bbox[0])* ((double) i)/n1 
	           +box->bbox[0]);

	snprintf(str, 999, "box%d_basis2", b);
	if( Getv(str, "ChebExtrema") )
	  box->v[var_Y][ijk] 
	    = 0.5*( (box->bbox[2] - box->bbox[3])*cos(j*PI/(n2-1)) 
	           +(box->bbox[2] + box->bbox[3])); 
        else if( Getv(str, "Fourier") )
	  box->v[var_Y][ijk] 
	    = ( (box->bbox[3] - box->bbox[2])* ((double) j)/n2 
	           +box->bbox[2]);

	snprintf(str, 999, "box%d_basis3", b);
	if( Getv(str, "ChebExtrema") )
	  box->v[var_Z][ijk] 
	    = 0.5*( (box->bbox[4] - box->bbox[5])*cos(k*PI/(n3-1)) 
	           +(box->bbox[4] + box->bbox[5])); 
        else if( Getv(str, "Fourier") )
	  box->v[var_Z][ijk] 
	    = ( (box->bbox[5] - box->bbox[4])* ((double) k)/n3 
	           +box->bbox[4]);
        }
      }
    }  /* end for (k = 0; k < n3; k++) */

    /* initialize diff. and filter matrices */
    {
       void (*get_coeffs)(double *,double *, int);
       void (*coeffs_of_deriv)(double, double, double *,double *, int);
       void (*eval_onPoints)(double *,double *, int);
       void (*filter_coeffs)(double *, int, int);
       int filt1, filt2, filt3;

       snprintf(str, 999, "box%d_filter1", b);
       filt1 = Geti(str);

       snprintf(str, 999, "box%d_filter2", b);
       filt2 = Geti(str);

       snprintf(str, 999, "box%d_filter3", b);
       filt3 = Geti(str);
  
       snprintf(str, 999, "box%d_basis1", b);
       if( Getv(str, "ChebExtrema") )
       {
         get_coeffs = cheb_coeffs_fromExtrema;
         coeffs_of_deriv = cheb_deriv;
         eval_onPoints = cheb_eval_onExtrema;
         filter_coeffs = cheb_filter;
       }
       else if( Getv(str, "Fourier") )
       {
         get_coeffs = four_coeffs;
         coeffs_of_deriv = four_deriv;
         eval_onPoints = four_eval;
         filter_coeffs = four_filter;
       }
       initdiffmatrix(box->bbox[0], box->bbox[1], box->D1, box->DD1, n1,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F1, n1+1-filt1, n1, 
                        get_coeffs, filter_coeffs, eval_onPoints);

        
       snprintf(str, 999, "box%d_basis2", b);
       if( Getv(str, "ChebExtrema") )
       {
         get_coeffs = cheb_coeffs_fromExtrema;
         coeffs_of_deriv = cheb_deriv;
         eval_onPoints = cheb_eval_onExtrema;
         filter_coeffs = cheb_filter;
       }
       else if( Getv(str, "Fourier") )
       {
         get_coeffs = four_coeffs;
         coeffs_of_deriv = four_deriv;
         eval_onPoints = four_eval;
         filter_coeffs = four_filter;
       }
       initdiffmatrix(box->bbox[2], box->bbox[3], box->D2, box->DD2, n2,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F2, n2+1-filt2, n2, 
                        get_coeffs, filter_coeffs, eval_onPoints);


       snprintf(str, 999, "box%d_basis3", b);
       if( Getv(str, "ChebExtrema") )
       {
         get_coeffs = cheb_coeffs_fromExtrema;
         coeffs_of_deriv = cheb_deriv;
         eval_onPoints = cheb_eval_onExtrema;
         filter_coeffs = cheb_filter;
       }
       else if( Getv(str, "Fourier") )
       {
         get_coeffs = four_coeffs;
         coeffs_of_deriv = four_deriv;
         eval_onPoints = four_eval;
         filter_coeffs = four_filter;
       }
       initdiffmatrix(box->bbox[4], box->bbox[5], box->D3, box->DD3, n3,
                      get_coeffs, coeffs_of_deriv, eval_onPoints);
       initfiltermatrix(box->F3, n3+1-filt3, n3, 
                        get_coeffs, filter_coeffs, eval_onPoints);
    }
  }

  if(1) printgrid(g);

  /* return pointer to newly created grid */
  return g;
}

