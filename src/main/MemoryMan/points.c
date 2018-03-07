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
  void *ret;
  int boxindex;
  tPointList *PL;

  PL=calloc(1, sizeof(*PL) );

  PL->grid=grid;

  ret = calloc(grid->nboxes, sizeof( *(PL->npoints) ) );
  if(ret==NULL) 
    errorexit("AllocatePointList: not enough memory for PL->npoints");
  PL->npoints = ret;

  ret = calloc(grid->nboxes, sizeof( *(PL->point) ) );
  if(ret==NULL) 
    errorexit("AllocatePointList: not enough memory for PL->point");
  PL->point = ret;
//printf("sizeof( *(PL->npoints) )=%d\n", sizeof( *(PL->npoints) ));
//printf("sizeof( *(PL->point) )=%d\n", sizeof( *(PL->point) ));

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

/* add bi to blist (if not already in it), returns n+1 if bi is new */
int addto_blist(int bi, int *blist, int n)
{
  int i;
  int addbi=1;
  /* add bi only if it is not already in blist */
  for(i=0; i<n; i++) if(blist[i]==bi) { addbi=0; break; }
  //printf(" n=%d\n",n);
  if(addbi) { blist[n] = bi; n++; }
  return n;
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
  PL->nblist = addto_blist(boxindex, PL->blist, PL->nblist);
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
  if(PL!=NULL)
  {
    printf("PointList=%p  PointList->grid=%p  PointList->nblist=%d\n",
           PL, PL->grid, PL->nblist);
    forallboxes(PL->grid,boxindex)
    {
      //if(PL->point[boxindex]==NULL) continue;
      printf("PointList->npoints[%d]=%d  ", boxindex, PL->npoints[boxindex]);
      printf("PointList->point[%d]=%p\n", boxindex, PL->point[boxindex]);

      printf("*(PointList->point[%d]) = ", boxindex);
      for(i=0; i<PL->npoints[boxindex]; i++)
        printf("%d ",PL->point[boxindex][i]);
      printf("\n");
    }
  }
  else printf("PointList=%p\n", PL);
}

/* print part of PointList in one box b */
void prPointList_inbox(tPointList *PL, int b)
{
  int i;
  if(PL!=NULL)
  {
    printf("PointList->grid=%p  PointList->nblist=%d  "
           "PointList->npoints[%d]=%d\n",
           PL->grid, PL->nblist, b, PL->npoints[b]);
    printf("*(PointList->point[%d]) = ", b);
    for(i=0; i<PL->npoints[b]; i++)
      printf("%d ",PL->point[b][i]);
    printf("\n");
  }
  else printf("PointList=%p\n", PL);
}

/* print part of PointList in one box b as i,j,k */
void prPointList_ijk_inbox(tPointList *PL, int b)
{
  int n;
  int n1 = PL->grid->box[b]->n1;
  int n2 = PL->grid->box[b]->n2;

  if(PL!=NULL)
  {
    //printf("PointList->grid=%p  PointList->nblist=%d  "
    //       "PointList->npoints[%d]=%d\n",
    //       PL->grid, PL->nblist, b, PL->npoints[b]);
    printf("PointList->nblist=%d  PointList->npoints[%d]=%d\n",
           PL->nblist, b, PL->npoints[b]);
    printf("*(PointList->point[%d]) = ", b);
    for(n=0; n<PL->npoints[b]; n++)
    {
      int ijk = PL->point[b][n];
      int k = kOfInd_n1n2(ijk,n1,n2);
      int j = jOfInd_n1n2_k(ijk,n1,n2,k);
      int i = iOfInd_n1n2_jk(ijk,n1,n2,j,k);
      printf("%d,%d,%d ", i,j,k);
    }
    printf("\n");
  }
  else printf("PointList=%p\n", PL);
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
