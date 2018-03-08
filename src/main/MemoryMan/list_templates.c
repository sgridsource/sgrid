/* main/main/list_templates.h */
/* Wolfgang Tichy, May 2017 */

#include "list_templates.h"
/* the list template is defined in list_templates.h */

/* To compile list functions with lists of type int, we need to use:
#define TYP int
#include "list_templates.c"
#undef TYP

To use them in another file we need to add
#define TYP int
#include "list_templates.h"
#undef TYP
to this file. */


/************************************************************************/
/* utility functions for these lists */


/* allocate an empty list */
LIST(TYP) *FN(alloc,LIST(TYP))(void)
{
  LIST(TYP) *u;

  u = calloc(1, sizeof(LIST(TYP)));
  return u;
} 

/* free a list */
void FN(free,LIST(TYP))(LIST(TYP) *u)
{
  if (u)
  {
    if (u->e) free(u->e);
    free(u);
    u = NULL;
  }
} 

/* clear a list */
void FN(clear,LIST(TYP))(LIST(TYP) *u)
{
  if (u)
  {
    if(u->e) free(u->e);
    u->e = NULL;
    u->n = 0;
  }
}

/* print list u */
void FN(pr,LIST(TYP))(LIST(TYP) *u)
{
  int i;
  double ui;
  printf("n=%d  e =", u->n);
  for(i=0; i<u->n; i++) { ui=u->e[i]; printf(" %g" , ui); }
  printf("\n");
}

/* add an entry to a list */
void FN(push,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  v->n += 1;
  v->e = realloc(v->e, sizeof(TYP) * v->n); 
  v->e[v->n-1] = vi;
}

/* append a list to a list */
void FN(pushlist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;

  if (!v || !u) return;
  v->n += u->n;
  v->e = realloc(v->e, sizeof(TYP) * v->n); 
  for (i = 0; i < u->n; i++)
    v->e[v->n - u->n + i] = u->e[i];
}

/* add to list (if not already in it) */
void FN(unionpush,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int addvi=1;
  /* add vi only if it is not already in blist */
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { addvi=0; break; }
  if(addvi) FN(push,LIST(TYP))(v, vi);
}

/* v = union(v, u): add all of u to list v (if not already in v) */
void FN(unionpushlist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;
  for(i=0; i<u->n; i++) FN(unionpush,LIST(TYP))(v, u->e[i]);
}

/* drop an entry from a list */
void FN(dropindex,LIST(TYP))(LIST(TYP) *v, int ind)
{
  int i;
  if(ind<0 || ind >= v->n) return;
  for(i = ind; i < v->n-1; i++)  v->e[i] = v->e[i+1];
  v->n -= 1;
}

/* drop an entry from a list */
void FN(drop,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  
  for (i = 0; i < v->n; i++)
    if (v->e[i] == vi)
    {
      v->n -= 1;
      for (; i < v->n; i++)
        v->e[i] = v->e[i+1];
      break;
    }
}

/* drop last n entries from a list */
void FN(droplastn,LIST(TYP))(LIST(TYP) *v, int n)
{
  if (n <= 0) 
    return;
  if (n >= v->n) 
    v->n = 0;
  else
    v->n -= n;
}

/* drop all in u from v */
void FN(droplist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;
  for(i=0; i<u->n; i++) FN(drop,LIST(TYP))(v, u->e[i]);
}

/* duplicate a list */
LIST(TYP) *FN(duplicate,LIST(TYP))(LIST(TYP) *v)
{
  int i;
  LIST(TYP) *u = FN(alloc,LIST(TYP))();

  for (i = 0; i < v->n; i++)
    FN(push,LIST(TYP))(u, v->e[i]);

  return u;
}

/* return 1 if vi is in list v */
int FN(in,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int in=0;
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { in=1; break; }
  return in;
}

/* return index of first element vi in list v, returns -1 if not in list */
int FN(index,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int in=-1; /* is not in list */
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { in=i; break; }
  return in;
}

/* return index of first element in list v that has prop returning 1,
   we start checking with element i0, returns -1 if prop returns 0 for all in v */
/* the function prop could be as simple as:
   int prop(int vi, void *p)
   {
     int *pi = (int *) p;
     return (vi == *pi);
   }
*/
int FN(index_prop,LIST(TYP))(LIST(TYP) *v, int i0,
                             int (*prop)(TYP vi, void *p), void *pars)
{
  int i;
  int in=-1; /* is not in list */
  if(i0<0) i0=0;
  for(i=i0; i<v->n; i++) if(prop(v->e[i], pars)) { in=i; break; }
  return in;
}
