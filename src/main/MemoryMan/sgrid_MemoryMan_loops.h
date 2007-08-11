/* sgrid_MemoryMan_loops.h */
/* Wolfgang Tichy, April 2005 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details.
*/

/****************************************************************************/
/* loop over all points in a box */
#define forallpoints(box,i) \
  for (i = 0; i < box->nnodes; i++)

/* loop over all boxes */
#define forallboxes(grid,boxindex) \
  for(boxindex=0; boxindex < grid->nboxes; boxindex++)

/* loop over all points in all boxes */
#define forallgridpoints(grid,curbox, boxindex,ijk) \
  forallboxes(grid,boxindex) \
    for(curbox=grid->box[boxindex], ijk=0; ijk < curbox->nnodes; ijk++)

/* loop over planes e.g. i=p plane */
#define  forplane1(i,j,k, n1,n2,n3, p) \
  for(i=(p), k = 0; k < (n3); k++) \
    for(     j = 0; j < (n2); j++)

#define  forplane2(i,j,k, n1,n2,n3, p) \
  for(j=(p), k = 0; k < (n3); k++) \
    for(     i = 0; i < (n1); i++)

#define  forplane3(i,j,k, n1,n2,n3, p) \
  for(k=(p), j = 0; j < (n2); j++) \
    for(     i = 0; i < (n1); i++)



/**************************************************************************/
/* point lists */
typedef struct tPL
{
  tGrid  *grid;  /* grid to which points belong */
  int *npoints;  /* npoints[b] = number of points in list in box b */
  int **point;   /* point[b] = array containing indices of all the points 
                    in the list in box b */ 
} tPointList;


/* loop over all points ijk in all boxes in list pointlist */
#define forPointList(pointlist, boxindex, pi , ijk) \
  forallboxes(pointlist->grid,boxindex) \
    for(pi = 0 , ijk = pointlist->point[boxindex][pi]; \
        pi < pointlist->npoints[boxindex]; \
        ijk = pointlist->point[boxindex][pi+1], pi++ )

/* loop over all points ijk in all boxes in list pointlist */
#define forPointList_inbox(pointlist, box, pi , ijk) \
  for(pi = 0 , ijk = pointlist->point[box->b][pi]; \
        pi < pointlist->npoints[box->b]; \
        ijk = pointlist->point[box->b][pi+1], pi++ )
