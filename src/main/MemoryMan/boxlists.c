/* boxlists.c */
/* Wolfgang Tichy 5/2017 */


#include "sgrid.h"
#include "MemoryMan.h"


/* ******************************************************************** */
/* WT: Here are some utility routines to use box primitive lists        */
/* We could make this more sphisticated like the point lists in point.c */
/* Note: the caller must make surte there is memory for these blists!   */
/* ******************************************************************** */

void pr_boxlist(int *blist, int n)
{
  int i;
  printf("%d boxes:", n);
  for(i=0; i<n; i++) printf(" %d", blist[i]);
  printf("\n");
}

/* add bi to boxlist (if not already in it), returns n+1 if bi is new */
int addto_boxlist(int bi, int **blist, int n)
{
  int i;
  int addbi=1;
  /* add bi only if it is not already in blist */
  for(i=0; i<n; i++) if((*blist)[i]==bi) { addbi=0; break; }
  //printf(" n=%d\n",n);
  if(addbi) { (*blist)[n] = bi; n++; }
  return n;
}                      

/* return number of boxes n and list blist with boxes that have an 
   attribute Attrib[iAttr] with value AttrVal. */
int boxlist_fromAttrib(tGrid *grid, int iAttr, int AttrVal, int **blist)
{
  int b, n=0;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] == AttrVal) { (*blist)[n] = b; n++; }
  }
  return n;
}

/* return number of boxes n2 and list blist2 with boxes that have an 
   attribute Attrib[iAttr] with value AttrVal and are in blist1. */
int boxlist2_ifAttrib(tGrid *grid, int *blist1, int n1,
                      int iAttr, int AttrVal, int **blist2)
{
  int i, n2=0;
  for(i=0; i<n1; i++)
  {
    int b = blist1[i];
    tBox *box = grid->box[b];
    if(box->Attrib[iAttr] == AttrVal) { (*blist2)[n2] = b; n2++; }
  }
  return n2;
}

/*
// make blist of all neighbors of box, return list len in n
int boxlist_boxneighbors(tBox *box, int **blist)
{
  int i, i2, n=0;
  for(i=0; i<box->nbfaces; i++)
  {
    int nb = box->bface[i]->obox->b; // get neighbor index nb
    int addnb = 1;
    // add nb only if it is not already in blist
    for(i2=0; i2<n; i2++) if((*blist)[i2]==nb) { addnb=0; break; }
    if(addnb) { (*blist)[n] = b; n++; }
  }
  return n;
}                      
*/
