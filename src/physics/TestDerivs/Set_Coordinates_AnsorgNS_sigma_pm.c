/* TestDerivs.c */
/* Wolfgang Tichy 2005 */


#include "sgrid.h"
#include "TestDerivs.h"


/* initialize Coordinates_AnsorgNS_sigma_pm */
int Set_Test_Coordinates_AnsorgNS_sigma_pm(tGrid *grid)
{
  int isigma;
  int b;

  if(!Getv("Coordinates_AnsorgNS_sigma_pm_vars", "yes"))
    return 0;

  isigma = Ind("Coordinates_AnsorgNS_sigma_pm");

  /* set Coordinates_AnsorgNS_sigma_pm in boxes where we use AnsorgNS */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *sigma;
    double sig1, c1;
    char str[1000];
    int i;

    snprintf(str, 999, "box%d_Coordinates", b);
    if( strstr(Gets(str), "AnsorgNS")==NULL ) continue;
    if( strstr(Gets(str), "AnsorgNS2") || strstr(Gets(str), "AnsorgNS3") )
    { sig1 = -1.0;  c1=0.30; }
    else
    { sig1 = +1.0;  c1=0.40; }
    printf("Setting %s in box%d:  sig1=%g c1=%g\n",
           VarName(isigma), b, sig1, c1);
    enablevar_inbox(box, isigma);
    sigma = box->v[isigma];

    forallpoints(box,i)
    {
      double B  = pY[i];
      double phi= pZ[i];

      /* set Coordinates_AnsorgNS_sigma_pm */
      sigma[i] = sig1 + c1*B*(1.0-B)*cos(phi);
    }
  }

  return 0;
}
