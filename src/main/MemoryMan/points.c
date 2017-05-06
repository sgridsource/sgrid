/* points.c */
/* Wolfgang Tichy 5/2005 */


#include "sgrid.h"
#include "MemoryMan.h"




/* ******************************************************************** */
/* WT: Here are some utility routines to use point lists of the         */
/*     type tPointList *, defined in sgrid_MemoryMan.h                  */
/* ******************************************************************** */

/* Allocate memory for a PointList */
tPointList *AllocatePointList(tGrid *grid)
{
 int boxindex;
 tPointList *PL;
 
 PL=calloc(1, sizeof(*PL) );

 PL->grid=grid;
 PL->npoints=calloc(grid->nboxes, sizeof( *(PL->npoints) ) );
 PL->point  =calloc(grid->nboxes, sizeof( *(PL->point) ) );
 PL->blist = NULL;
 PL->nblist= 0;
 forallboxes(grid,boxindex)
 {
   PL->npoints[boxindex] = 0;
   PL->point[boxindex] = NULL;
 }

 // printf("AllocatePointList: PL=%p\n",PL); 
 return PL;
}


/* add one point to a PointList on box boxindex */
void AddToPointList(tPointList *PL, int boxindex, int newpoint)
{
 void *ret;

 ret=realloc(PL->point[boxindex], 
             (sizeof( *(PL->point[boxindex]) ))*(PL->npoints[boxindex]+2) );
 if(ret==NULL) 
   errorexit("AddToPointList: not enough memory for PL->point[boxindex]");
 
 PL->point[boxindex]=ret;
 PL->point[boxindex][PL->npoints[boxindex]]=newpoint;
 PL->npoints[boxindex]++;

 ret=realloc(PL->blist, ( sizeof(*(PL->blist)) )*(PL->nblist+2));
 if(ret==NULL) 
   errorexit("AddToPointList: not enough memory for PL->blist");
 PL->blist  = ret;
 //printf("PL->blist=%p\n", PL->blist);
 //printf("ret=%p\n", ret);
 //printf("%d\n",( sizeof(*(PL->blist)) )*(PL->nblist+255));
 PL->nblist = addto_boxlist(boxindex, &(PL->blist), PL->nblist);
}


/* free a PointList */
void FreePointList(tPointList *PL)
{
 if(PL!=NULL)
 {
   int boxindex;

   forallboxes(PL->grid,boxindex)
      free(PL->point[boxindex]);

   free(PL->npoints);
   free(PL->point);
   free(PL->blist);
   free(PL);
 }
}


/* print a PointList */
void prPointList(tPointList *PL)
{
  int i;
  int boxindex;
  
  printf("PointList=%p  PointList->grid=%p\n", PL, PL->grid);
  pr_boxlist(PL->blist, PL->nblist);  
  forallboxes(PL->grid,boxindex)
  {
     printf("PointList->npoints[%d]=%d  ", boxindex, PL->npoints[boxindex]);
     printf("PointList->point[%d]=%p\n", boxindex, PL->point[boxindex]);
  
     printf("*(PointList->point[%d]) = ", boxindex);
     for(i=0; i<PL->npoints[boxindex]; i++)
       printf("%d ",PL->point[boxindex][i]);
     printf("\n");
  }
}

/* return a duplicate of PL0 but put grid into it */
tPointList *DuplicatePointList_for_grid(tPointList *PL0, tGrid *grid)
{
  tPointList *PL;
  int b, i;
  void *ret;

  if(PL0==NULL) return NULL;
  PL = AllocatePointList(grid);

  /* copy all point indices */
  forallboxes(grid, b)
  {
    ret=realloc(PL->point[b], (sizeof( *(PL->point[b]) ))*(PL0->npoints[b]+2));
    if(ret==NULL) 
      errorexit("DuplicatePointList: not enough memory for PL->point[b]");
    PL->point[b] = ret;
    PL->npoints[b] = PL0->npoints[b];
    for(i=0; i<PL->npoints[b]; i++) PL->point[b][i] = PL0->point[b][i];
  }

  /* copy the blist */
  ret=realloc(PL->blist, ( sizeof(*(PL->blist)) )*(PL->nblist+2));
  if(ret==NULL) 
    errorexit("DuplicatePointList: not enough memory for PL->blist");
  PL->blist  = ret;
  PL->nblist = PL0->nblist;
  for(i=0; i<PL->nblist; i++) PL->blist[i] = PL0->blist[i];

  return PL;
}

/* return a duplicate of PL0 */
tPointList *DuplicatePointList(tPointList *PL0)
{
  tPointList *PL;
  if(PL0==NULL) return NULL;
  PL = DuplicatePointList_for_grid(PL0, PL0->grid);
  return PL;
}
