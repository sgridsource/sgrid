/* boxlists.c */
/* Wolfgang Tichy 5/2017 */


#include "sgrid.h"
#include "MemoryMan.h"


/* ****************************************************************************** */
/* WT: Here are some utility routines to use box lists (bl) of type intList       */
/* ****************************************************************************** */

void pr_boxlist(intList *bl)
{
  int i;
  printf("%d boxes:", bl->n);
  for(i=0; i<bl->n; i++) printf(" %d", bl->e[i]);
  printf("\n");
}

/* return number of boxes n and list bl with boxes that have an
   attribute Attrib[iAttr] with value AttrVal. */
int bladd_ifAttrib(tGrid *grid, int iAttr, int AttrVal, intList *bl)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] == AttrVal)  push_intList(bl, b);
  }
  return bl->n;
}

/* return number of boxes n and list bl with boxes that do not have an
   attribute Attrib[iAttr] with value AttrVal. */
int bladd_ifnotAttrib(tGrid *grid, int iAttr, int AttrVal, intList *bl)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] != AttrVal)  push_intList(bl, b);
  }
  return bl->n;
}

/* return number of boxes and list bl2 with boxes that have an
   attribute Attrib[iAttr] with value AttrVal and are in bl1. */
int bl2_ifAttrib(tGrid *grid, intList *bl1, int iAttr, int AttrVal, intList *bl2)
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

/* return number of boxes and list bl2 with boxes that do not have an
   attribute Attrib[iAttr] with value AttrVal and are in bl1. */
int bl2_ifnotAttrib(tGrid *grid, intList *bl1, int iAttr, int AttrVal, intList *bl2)
{
  int i;
  for(i=0; i<bl1->n; i++)
  {
    int b = bl1->e[i];
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] != AttrVal)  push_intList(bl2, b);
  }
  return bl2->n;
}

/* make list bl of all neighbors of box, return list len */
int bladd_neighbors(tBox *box, intList *bl)
{
  int i, i2;
  for(i=0; i<box->nbfaces; i++)
  {
    int nb, addnb;

    nb = box->bface[i]->ob; /* get neighbor index nb */
    if(nb<0) continue;      /* do nothing if there is no neighbor */

    /* add nb only if it is not already in blist */
    addnb = 1;
    for(i2=0; i2<bl->n; i2++) if(bl->e[i2]==nb) { addnb=0; break; }
    if(addnb) push_intList(bl, nb);
  }
  return bl->n;
}                      

/* make list bl of all neighbors of box, return list len */
int bladd_neighbors_of_neighbors(tBox *box, intList *bl)
{
  int i;

  intList *bln  = alloc_intList();
  intList *bln2 = alloc_intList();
  intList *blnn = alloc_intList();

  /* neigbors of box */
  bladd_neighbors(box, bln);

  /* go over neighbors in bln */
  for(i=0; i<bln->n; i++)
  {
    tBox *box2 = box->grid->box[bln->e[i]]; /* a neighbor box */

    /* put neighbors of box2 into bln2 */
    clear_intList(bln2);
    bladd_neighbors(box2, bln2);

    /* add all in bln2 to blnn */
    unionpushlist_intList(blnn, bln2);
  }

  /* remove all neigbors in bln from blnn */
  droplist_intList(blnn, bln);

  /* drop box from list blnn */
  drop_intList(blnn, box->b);

  /* add all next nearest neighbors in blnn to bl */
  pushlist_intList(bl, blnn);

  free_intList(blnn);
  free_intList(bln2);
  free_intList(bln);

  return bl->n;
}
