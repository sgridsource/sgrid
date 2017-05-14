/* boxlists.c */
/* Wolfgang Tichy 5/2017 */


#include "sgrid.h"
#include "MemoryMan.h"


/* ************************************************************************* */
/* WT: Here are some utility routines to use box lists of type intList       */
/* ************************************************************************* */

void pr_boxlist(intList *bl)
{
  int i;
  printf("%d boxes:", bl->n);
  for(i=0; i<bl->n; i++) printf(" %d", bl->e[i]);
  printf("\n");
}

/* return number of boxes n and list bl with boxes that have an
   attribute Attrib[iAttr] with value AttrVal. */
int boxlist_fromAttrib(tGrid *grid, int iAttr, int AttrVal, intList *bl)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] == AttrVal)  push_intList(bl, b);
  }
  return bl->n;
}

/* return number of boxes n and list bl2 with boxes that have an
   attribute Attrib[iAttr] with value AttrVal and are in bl1. */
int boxlist2_ifAttrib(tGrid *grid, intList *bl1,
                      int iAttr, int AttrVal, intList *bl2)
{
  int i;
  for(i=0; i<bl1->n; i++)
  {
    int b = bl1->e[i];
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] == AttrVal)  push_intList(bl2, b);
  }
  return bl2->n;
}

/* make list bl of all neighbors of box, return list len */
int boxlist_neighbors(tBox *box, intList *bl)
{
  int i, i2;
  for(i=0; i<box->nbfaces; i++)
  {
    int nb = box->bface[i]->ob; // get neighbor index nb
    int addnb = 1;
    // add nb only if it is not already in blist
    for(i2=0; i2<bl->n; i2++) if(bl->e[i2]==nb) { addnb=0; break; }
    if(addnb) push_intList(bl, nb);
  }
  return bl->n;
}                      
