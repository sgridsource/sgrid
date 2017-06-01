/* Singularities.c */
/* Wolfgang Tichy 5/2017 */

#include "sgrid.h"
#include "Coordinates.h"


/* print tSingInfo */
void prSingInfo(tSingInfo *si)
{
  int i,j;

  if(si==NULL) { printf("si = NULL\n"); return; }

  printf("si->f[j] = ");
  for(j=0; j<6; j++) printf("%d ", si->f[j]);
  printf("\n");
  printf("si->x[j] = ");
  for(j=1; j<=3; j++) printf("%c ", si->x[j]);
  printf("\n");
  printf("si->dx_dX[i][j] =\n");
  for(i=1; i<=3; i++)
  {
    for(j=1; j<=3; j++) printf("%c ", si->dx_dX[i][j]);
    printf("\n");
  }
}

/* set SingInfo for domains 1 and 2 of AnsorgNS */
int isSing_AnsorgNS12(void *aux, double X, double Y, double Z,
                      int update, tSingInfo *si)
{
  int sing=0;
  int i;

  if( dequal(Y, 0.0) || dequal(Y, 1.0) ) /* on axis */
  {
    if(update)
    {
      for(i=0; i<6; i++) si->f[i] = 0;
      si->f[2+dequal(Y, 1.0)] = 1;

      si->x[1] = '.';
      si->x[2] = si->x[3] = '0';

      si->dx_dX[1][1] = '.';
      si->dx_dX[2][1] = si->dx_dX[3][1] = '0';
      for(i=1; i<=3; i++) si->dx_dX[i][2] = '.';
      for(i=1; i<=3; i++) si->dx_dX[i][3] = '0';
      if( dequal(X, 1.0) && dequal(Y, 0.0) ) /* at infinity */
      {
        si->x[1] = '0';
        si->x[2] = si->x[3] = 'i';
      }
    }
    sing=1;
  }
  return sing;
}

/* set SingInfo for domains 3 and 4 of AnsorgNS */
int isSing_AnsorgNS03(void *aux, double X, double Y, double Z,
                      int update, tSingInfo *si)
{
  int sing=0;
  int i;

  /* X=0 is deactivated, because we don't use grids where X=0 */
  if( 0 && dequal(X, 0.0) )
  {
    if(update)
    {
      for(i=1; i<6; i++) si->f[i] = 0;
      si->f[0] = 1;
      for(i=1; i<=3; i++)   si->dx_dX[i][1] = '.';
      for(i=1; i<=3; i++) { si->dx_dX[i][2] = si->dx_dX[i][3] = '0'; }
    }
    sing=1;
  }
  if( dequal(Y, 0.0) || dequal(Y, 1.0) )
  {
    if(update)
    {
      for(i=0; i<6; i++) si->f[i] = 0;
      si->f[2+dequal(Y, 1.0)] = 1;

      si->dx_dX[1][1] = '.';
      si->dx_dX[2][1] = si->dx_dX[3][1] = '0';
      if(!sing)
      {
        for(i=1; i<=3; i++) si->dx_dX[i][2] = '.';
        for(i=1; i<=3; i++) si->dx_dX[i][3] = '0';
      }
    }
    sing=1;
  }
  return sing;
}
