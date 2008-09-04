/* YlmFilter.c */
/* Wolfgang Tich 8/2007 */
/* use s2kit to use the naive Plm trafos to filter unew */

#include "sgrid.h"
#include "s2kit.h"

#include "s2kit_pmls.h"
#include "s2kit_makeweights.h"
#include "s2kit_naive_synthesis.h"


/* use s2kit's Plm to project onto Ylm */
/* here lmax = bw-1 + lmshift, where bw = n2/4 */
/* currently we need: lmshift<=0 */
void Naive_YlmFilter_lmshift(tVarList *unew, int lmshift)
{
  tGrid *grid = unew->grid;
  int b;

  if(lmshift>0) errorexit("YlmFilter: we need lmshift<=0.");

  /* loop over boxes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, m, vi;
    int l, bw;
    double *samples, *coeffs;
    double **plm, *weights, *workspace;

    /* we need an even n3 and an n2 divisible by 4 */
    if(n2 % 4)
      errorexit("YlmFilter: box->n2 has to be divisible by 4.");
    if(n3 % 2)
      errorexit("YlmFilter: box->n3 has to be divisible by 2.");

    /* set sampling rate B=bw */
    bw = n2/4;

    /* space for samples and coefficients */
    samples   = (double *) malloc(sizeof(double) * n2);
    coeffs    = (double *) malloc(sizeof(double) * bw);

    /* space for precomputed Plms */
    plm = (double **) malloc(sizeof(double *) * bw );
    for(m=0; m<bw; m++)
      plm[m] = (double *) malloc(sizeof(double) * 2 * bw * bw );

    /* for weights */
    weights = (double *) malloc(sizeof(double) * 4 * bw);

    /* workspace space */
    workspace = (double *) malloc(sizeof(double) * 18 * bw);

    /* precompute the Plms */
    for(m=0; m<bw; m++)
      PmlTableGen(bw, m, plm[m], workspace);

    /* make the weights */
    makeweights( bw, weights );


    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *vc  = box->v[Ind("temp1")];
//printf("Naive_YlmFilter: %s\n", VarName(unew->index[vi]));

      /* compute coeffs of var after Fourier trafo in phi*/
      spec_analysis1(box, 3, box->Mcoeffs3, var, vc);

      /* loop over all phi-coeffs, i.e. all k or m */
      for(k=0; k<n3; k++)
      {
        m = (k+1)/2; /* set m */
        if(m<bw + lmshift)
        {
          /* loop over all radii */
          for(i=0; i<n1; i++)
          {
            /* get sample for a particular radius and m (i.e. k) */
            get_memline(vc, samples, 2, i, k, n1, n2, n3);
            /* Note: samples need to be taken at the points:
               theta_j = PI*(2j + 1)/(4bw), where j = 0,1, ... ,2bw-1
               these are the points of SphericalDF if n2%4=0 */

            /* now do forward naive transform */
            Naive_AnalysisX(samples, bw, m, weights, coeffs, plm[m], workspace);
            /* Note l_max=bw-1=n2/4-1, so we filter out half of all l modes!!! */

            /* set coeffs with l = l_max + 1 + lmshift, ..., lmax to zero */
            for(l=bw-m + lmshift; l<bw-m; l++) coeffs[l]=0.0;
            /* coeffs is an array of size (bw-m). First coefficient 
               is coefficient for Pmm */

            /* do inverse naive transform */
            Naive_SynthesizeX(coeffs, bw, m, samples, plm[m]);
            /* Note: this writes in samples only for theta<PI !!!*/

            /* put filtered sample back */
            put_memline(vc, samples, 2, i, k, n1, n2, n3);
          }
        }
        else /* if(m>=bw + lmshift) */
        {
          /* filter all k modes with m>=bw=l+1 */
          for(j=0; j<n2; j++)
            for(i=0; i<n1; i++)
              vc[Index(i,j,k)]=0.0;
        }
      }
      /* use modified coeffs to change var */
      spec_synthesis1(box, 3, box->Meval3, var, vc);
      /* Note: vc was wrong for all theta>PI !!! */

      /* copy var into double covered regions */
      for(k = 0;    k < n3/2; k++)
      for(j = n2/2; j < n2; j++)
      for(i = 0;    i < n1; i++)
        var[Index(i,j,k)] = var[Index(i,n2-j-1,k+n3/2)];
      for(k = n3/2; k < n3; k++)
      for(j = n2/2; j < n2; j++)
      for(i = 0;    i < n1; i++)
        var[Index(i,j,k)] = var[Index(i,n2-j-1,k-n3/2)];
    }

    /* now free up memory for naive trafo */
    for(m=0; m<bw; m++) free(plm[m]);
    free(plm);
    free(workspace);
    free(weights);
    free(coeffs);
    free(samples);
  }
}    

/* wrapper that sets lmshift=0, i.e. lmax = bw-1 = n2/4 - 1 */
void Naive_YlmFilter(tVarList *unew)
{
  Naive_YlmFilter_lmshift(unew, 0);
}
