/* find_ijk.c */
/* Wolfgang Tichy 9/2018 */

#include "sgrid.h"
#include "Coordinates.h"


/* find index of X-coord, returns -1 if not found */
int find_i_Of_X(tBox *box, double X0)
{
  double *X = box->v[Ind("X")];
  int i;
  for(i=0; i<box->n1; i++)
    if(dequal(X[i], X0)) return i;
  return -1;
}

/* find index of Y-coord, returns -1 if not found */
int find_j_Of_Y(tBox *box, double Y0)
{
  double *Y = box->v[Ind("Y")];
  int j;
  for(j=0; j<box->n2; j++)
    if(dequal(Y[j], Y0)) return j;
  return -1;
}

/* find index of Z-coord, returns -1 if not found */
int find_k_Of_Z(tBox *box, double Z0)
{
  double *Z = box->v[Ind("Z")];
  int k;
  for(k=0; k<box->n3; k++)
    if(dequal(Z[k], Z0)) return k;
  return -1;
}
