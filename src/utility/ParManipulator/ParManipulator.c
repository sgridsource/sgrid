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

    /* set the n1,n2,n3 of in each box */
    snprintf(str, 999, "box%d_n1", b);
    if(Getv(str, "n1")) Sets(str, Gets("n1"));
    snprintf(str, 999, "box%d_n2", b);
    if(Getv(str, "n2")) Sets(str, Gets("n2"));
    snprintf(str, 999, "box%d_n3", b);
    if(Getv(str, "n3")) Sets(str, Gets("n3"));

    /* set the min1,min2,min3 of in each box */
    snprintf(str, 999, "box%d_min1", b);
    if(Getv(str, "min1")) Sets(str, Gets("min1"));
    snprintf(str, 999, "box%d_min2", b);
    if(Getv(str, "min2")) Sets(str, Gets("min2"));
    snprintf(str, 999, "box%d_min3", b);
    if(Getv(str, "min3")) Sets(str, Gets("min3"));

    /* set the max1,max2,max3 of in each box */
    snprintf(str, 999, "box%d_max1", b);
    if(Getv(str, "max1")) Sets(str, Gets("max1"));
    snprintf(str, 999, "box%d_max2", b);
    if(Getv(str, "max2")) Sets(str, Gets("max2"));
    snprintf(str, 999, "box%d_max3", b);
    if(Getv(str, "max3")) Sets(str, Gets("max3"));

    /* set the basis1,basis2,basis3 of in each box */
    snprintf(str, 999, "box%d_basis1", b);
    if(Getv(str, "basis1")) Sets(str, Gets("basis1"));
    snprintf(str, 999, "box%d_basis2", b);
    if(Getv(str, "basis2")) Sets(str, Gets("basis2"));
    snprintf(str, 999, "box%d_basis3", b);
    if(Getv(str, "basis3")) Sets(str, Gets("basis3"));

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

