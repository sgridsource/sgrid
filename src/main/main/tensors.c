/* tensors.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 6/02, 11/02 */

#include "sgrid.h"
#include <ctype.h>



/* helper function: 
   decode tensor index string
   return list of indices into variable list 
   return sign under reflections for symmetry boundaries

   should be made automatic, but for now this is simpler 
*/
void tensorindexlist(char *t, int *nilist, char **ilist, int *sym)
{
  /* name of coordinates, could be made variable */
  char *coord[3] = {"x", "y", "z"};
  int i, j, k;
  int n = 0;
  char *tensorindices = strdup(t);

  /* convert local copy to lower case since we ignore co/contra-variance */
  for (t = tensorindices; *t; t++)
    *t = tolower(*t);

  /* initialize symmetries */
  for (i = 0; i < 3*NINDEXLIST; i++) 
    sym[i] = 1;


  /* now treat each case separately */
  if (strcmp(tensorindices, "") == 0) {
    ilist[n] = calloc(sizeof(char), 8);
    sprintf(ilist[n++], "%s", "");
  }
  
  if (strcmp(tensorindices, "i") == 0) {
    for (i = 0; i < 3; i++) {
      sym[3*n+i] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord[i]);
    }
  }
  
  if (strcmp(tensorindices, "ij") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }

  if (strcmp(tensorindices, "ij+ji") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }
  
  if (strcmp(tensorindices, "ij-ji") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i+1; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }
  
  if (strcmp(tensorindices, "ijk") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk+ikj") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) 
    for (k = j; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk+jik") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk-jik") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i+1; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (n == 0) {
    printf("Error in index string %s.\n", tensorindices);
    printf("Legal combinations besides the empty string are\n");
    printf("i, ij, ij+ji, ij-ji, ijk, ijk+ikj, ijk+jik\n");
    printf("Anything else can be easily added to main/tensors.c.\n");
    errorexit("");
  }

  *nilist = n;
  free(tensorindices);
}

