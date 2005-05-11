/* points.c */
/* Wolfgang Tichy 5/2005 */


#include "sgrid.h"
#include "MemoryMan.h"




/* ******************************************************************** */
/* WT: Here are some utility routines to use point lists of the         */
/*     type tPointList *, defined in sgrid_MemoryMan_loops.h            */
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
   free(PL);
 }
}


/* print a PointList */
void prPointList(tPointList *PL)
{
  int i;
  int boxindex;
  
  printf("PointList=%p  PointList->grid=%p\n", PL, PL->grid);
  
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
