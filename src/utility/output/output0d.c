/* output0d.c */
/* Wolfgang Tichy, August 2005 */

#include "sgrid.h"
#include "output.h"



/* output one value */
void output0d_value(char *filename, double time, double val)
{
  FILE *fp;
  
  /* open file */
  fp = fopen(filename, "a");
  if (!fp) errorexits("failed opening %s", filename);

  /* write value */
  fprintf(fp, "%.16g %.16g\n", time, val);

  /* close file */
  fclose(fp);
}


/* do 0d output with one var */
void output0d_boxvar(tBox *box, char *name)
{
  char filename[1000];
  char str[1000];
  double max, min, rms, meanAbs, mean, VolInt;
  double Vol, Int;
  double (*VolIntergral)(tBox *box, double *u, double *U)=NULL;
  double *temp1  = box->v[Ind("temp1")];
  double *temp2  = box->v[Ind("temp2")];
  double *VolJac = box->v[Ind("temp3")]; /* Jacobian or volume element */
  double *var    = box->v[Ind(name)];
  int i;
  double time=box->grid->time;

  /* find max and min of var */
  max = min = var[0];
  forallpoints(box ,i)
  {
    if(var[i]>max) max=var[i];
    if(var[i]<min) min=var[i];
  }

  /* determine function used to compute volume integrals */
  snprintf(str, 999, "box%d_Coordinates", box->b);
  if( Getv(str, "Cartesian") )
  {
    VolIntergral = spec_3dIntegral;
    forallpoints(box ,i) VolJac[i] = 1.0;
  }
  else if( Getv(str, "SphericalDF") )
  {
    VolIntergral = spec_sphericalDF3dIntegral;
    forallpoints(box ,i) VolJac[i] = 1.0; /* not r^2 sin(theta), since spec_sphericalDF3dIntegral is special */
  }
  else
  {
    prdivider(0);
    printf("WARNING!!!\n");
    printf("output0d_boxvar: I don't know how to do volume integrals in\n"
           "%s coordinates...\n", Gets(str));
    VolIntergral = spec_3dIntegral;
    forallpoints(box ,i) VolJac[i] = 1.0;
    printf("Defaulting to same method as for Cartesian coordinates.\n");
    printf("WARNING: 0doutput of %s may be meaningless!!!\n", name);
    prdivider(0);
  }
          
  /* compute volume */
  Vol=VolIntergral(box, VolJac, temp2);

  /* integrate var^2 and compute rms */
  forallpoints(box ,i) temp1[i] = var[i]*var[i]*VolJac[i];
  Int = VolIntergral(box, temp1, temp2);
  rms = sqrt(Int/Vol);

  /* integrate |var| and compute meanAbs */
  forallpoints(box ,i) temp1[i] = fabs(var[i])*VolJac[i];
  Int = VolIntergral(box, temp1, temp2);
  meanAbs = Int/Vol;

  /* integrate var and compute VolInt and mean */
  forallpoints(box ,i) temp1[i] = var[i]*VolJac[i];
  VolInt = VolIntergral(box, temp1, temp2);
  mean = VolInt/Vol;
  
  /* output max, min, rms, meanAbs, mean, VolInt */
  snprintf(filename, 999, "%s/%s_max.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, max);

  snprintf(filename, 999, "%s/%s_min.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, min);

  snprintf(filename, 999, "%s/%s_rms.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, rms);

  snprintf(filename, 999, "%s/%s_meanAbs.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, meanAbs);

  snprintf(filename, 999, "%s/%s_mean.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, mean);

  snprintf(filename, 999, "%s/%s_VolInt.%d", Gets("outdir"), name, box->b);
  output0d_value(filename, time, VolInt);
}
