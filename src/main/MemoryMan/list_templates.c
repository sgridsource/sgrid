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
void FN(unionof,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int addbi=1;
  /* add bi only if it is not already in blist */
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { addbi=0; break; }
  FN(push,LIST(TYP))(v, vi);
}

/* drop an entry from a variable list */
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

/* duplicate a list */
LIST(TYP) *FN(duplicate,LIST(TYP))(LIST(TYP) *v)
{
  int i;
  LIST(TYP) *u = FN(alloc,LIST(TYP))();

  for (i = 0; i < v->n; i++)
    FN(push,LIST(TYP))(u, v->e[i]);

  return u;
}
