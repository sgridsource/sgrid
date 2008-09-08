/* ParManipulator.c */
/* Wolfgang Tichy 9/2008 */

#include "sgrid.h"
#include "ParManipulator.h"


/* set some pars automatically in PRE_GRID */
int ParManipulator_setPars(tGrid *grid)
{
  int nboxes = Geti("nboxes");
  int b;

  for(b=0; b<nboxes; b++)
  {
    char str[1000];
    char val[1000];

    /* set the n1,n2,n3 of in each box */
    snprintf(str, 999, "box%d_n1", b);
    if(Getv(str, "n1")) 
    {
      snprintf(val, 999, "%s n1", Gets("n1"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_n2", b);
    if(Getv(str, "n2")) 
    {
      snprintf(val, 999, "%s n2", Gets("n2"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_n3", b);
    if(Getv(str, "n3")) 
    {
      snprintf(val, 999, "%s n3", Gets("n3"));
      Sets(str, val);
    }

    /* set the min1,min2,min3 of in each box */
    snprintf(str, 999, "box%d_min1", b);
    if(Getv(str, "min1")) 
    {
      snprintf(val, 999, "%s min1", Gets("min1"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_min2", b);
    if(Getv(str, "min2")) 
    {
      snprintf(val, 999, "%s min2", Gets("min2"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_min3", b);
    if(Getv(str, "min3")) 
    {
      snprintf(val, 999, "%s min3", Gets("min3"));
      Sets(str, val);
    }

    /* set the max1,max2,max3 of in each box */
    snprintf(str, 999, "box%d_max1", b);
    if(Getv(str, "max1")) 
    {
      snprintf(val, 999, "%s max1", Gets("max1"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_max2", b);
    if(Getv(str, "max2")) 
    {
      snprintf(val, 999, "%s max2", Gets("max2"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_max3", b);
    if(Getv(str, "max3")) 
    {
      snprintf(val, 999, "%s max3", Gets("max3"));
      Sets(str, val);
    }

    /* set the basis1,basis2,basis3 of in each box */
    snprintf(str, 999, "box%d_basis1", b);
    if(Getv(str, "basis1")) 
    {
      snprintf(val, 999, "%s basis1", Gets("basis1"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_basis2", b);
    if(Getv(str, "basis2")) 
    {
      snprintf(val, 999, "%s basis2", Gets("basis2"));
      Sets(str, val);
    }
    snprintf(str, 999, "box%d_basis3", b);
    if(Getv(str, "basis3")) 
    {
      snprintf(val, 999, "%s basis3", Gets("basis3"));
      Sets(str, val);
    }

    /* set dt */
    if(Getv("ParManipulator_dt_scaling", "1/n1^2"))
    {
      int n1 = Geti("n1");
      double coeff = Getd("ParManipulator_ds_min_coeff");
      Setd("ParManipulator_ds_min", coeff/(n1*n1));
      Setd("dt", Getd("ParManipulator_dtfac")*Getd("ParManipulator_ds_min"));
    }
  }

  return 0;
}

