/* s2kit_test.c */
/* Wolfgang Tichy, 2008 */
/* test Naive_YlmFilter(vlfil); */
  

#include "sgrid.h"
#include "s2kit.h"



/* test Naive_YlmFilter(vlfil) */
int s2kit_test_Naive_YlmFilter(tGrid *grid) 
{
  tVarList *vlfil;
  tVarList *vlvar;
  tVarList *vldif;
  int b;
  int lmshift = Geti("s2kit_test_lmshift");

  /* I need to enable my vars */
  //enablevar(grid, Ind("temp1"));
  enablevar(grid, Ind("s2kit_test_var1"));
  enablevar(grid, Ind("s2kit_test_var2"));
  enablevar(grid, Ind("s2kit_test_fil1"));    
  enablevar(grid, Ind("s2kit_test_fil2"));
  enablevar(grid, Ind("s2kit_test_dif1"));    
  enablevar(grid, Ind("s2kit_test_dif2"));

  /* set var1 and var2 */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i;
    char str[1000];
    int n2 = box->n2;
    double *v1 = box->v[Ind("s2kit_test_var1")];
    double *v2 = box->v[Ind("s2kit_test_var2")];
    double *f1 = box->v[Ind("s2kit_test_fil1")];
    double *f2 = box->v[Ind("s2kit_test_fil2")];
    double *X = box->v[Ind("X")];
    double *Y = box->v[Ind("Y")];
    double *Z = box->v[Ind("Z")];

    snprintf(str, 999, "box%d_Coordinates", b);
    if(!Getv(str, "SphericalDF")) errorexit("use SphericalDF");
    if(n2%4) errorexit("use n2 divisible by 4.");

    printf("s2kit_test_Naive_YlmFilter: Set up var1 to contain l=2 only.\n");
    printf("s2kit_test_Naive_YlmFilter: Set up var2 to contain l=2 and l=3.\n");
    forallpoints(box,i)
    {
      double r     = X[i];
      double theta = Y[i] + PI/((1+n2%2)*n2);
      double phi   = Z[i];
      double m20 = 3*cos(theta)*cos(theta)-1.0;      /* l=2 m=0 mode */
      double m21 = sin(theta)*cos(theta)*sin(phi);   /* l=2 m=1 mode */
      double m22 = sin(theta)*sin(theta)*sin(2*phi); /* l=2 m=2 mode */
      double m30 = 5*cos(theta)*cos(theta)*cos(theta)-3*cos(theta);   /* l=3 m=0 mode */
      double m31 = (5*cos(theta)*cos(theta)-1.0)*sin(theta)*sin(phi); /* l=3 m=1 mode */
      double m32 = sin(theta)*sin(theta)*cos(theta)*sin(2*phi); /* l=3 m=2 mode */
      double m33 = sin(theta)*sin(theta)*sin(theta)*sin(3*phi); /* l=3 m=3 mode */
      
      /* set fil1=var1 and fil2=var2 */
      f1[i] = v1[i] = r*(m20 + m21 + m22); /* a l=2 mode */
      //f1[i] = v1[i] = sqrt(fabs(theta*phi));
      f2[i] = v2[i] = r*(m20 + m21 + m22 + 
                         m30 + m31 + m32 + m33); /* both l=2 and 3 */
    }
  }

  /* make var list to filter */
  vlfil = vlalloc(grid);
  vlpush(vlfil, Ind("s2kit_test_fil1"));
  vlpush(vlfil, Ind("s2kit_test_fil2"));

  /* now filter vlfil */
  printf("Naive_YlmFilter_lmshift: filtering var1 and var2 with lmshift=%d, "
         "and writing results in fil1 and fil2\n", lmshift);
  Naive_YlmFilter_lmshift(vlfil, lmshift);
  
  /* compute difference */
  vlvar = vlalloc(grid);
  vlpush(vlvar, Ind("s2kit_test_var1"));
  vlpush(vlvar, Ind("s2kit_test_var2"));
  vldif = vlalloc(grid);
  vlpush(vldif, Ind("s2kit_test_dif1"));
  vlpush(vldif, Ind("s2kit_test_dif2"));
  /* subtract: vldif = vlvar-vlfil */
  printf("computing dif1/2 = var1/2 - fil1/2.\n");
  vladd(vldif, 1.0,vlvar, -1.0,vlfil);

  vlfree(vlfil);
  vlfree(vlvar);
  vlfree(vldif);
  return 0;
}
