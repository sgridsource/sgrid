/* sgrid_MemoryMan_loops.h */
/* Wolfgang Tichy, April 2005 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details.
*/

/****************************************************************************/
/* loop over all points in a box */
#define forallpoints(box,i) \
  for (i = 0; i < box->nnodes; i++)

/* loop over all boxes no matter if they are active or not */
#define forallActiveAndInactiveboxes(grid,boxindex) \
  for(boxindex=0; boxindex < grid->nboxes; boxindex++)

/* number of active boxes */
#define activeboxes(grid) ((grid)->nboxes - (grid)->ninactive)

/* loop over all boxes that are active, same as: forallActiveboxes */
#define forallboxes(grid,boxindex) \
  for(boxindex=0; boxindex < grid->nboxes; boxindex++) \
    if(!(grid->box[boxindex]->inactive))

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

/* same as forplane1/2/3, but we can specify the plane number N */
#define  forplaneN(N, i,j,k, n1,n2,n3, p) \
  for(k=(p)*((N)==3); ( k<(n3) ) && ( ( k==(p) ) || (N)!=3 ); k++) \
  for(j=(p)*((N)==2); ( j<(n2) ) && ( ( j==(p) ) || (N)!=2 ); j++) \
  for(i=(p)*((N)==1); ( i<(n1) ) && ( ( i==(p) ) || (N)!=1 ); i++)

/* same as forplaneN, but omit edges */
#define  forinnerplaneN(N, i,j,k, n1,n2,n3, p) \
  for(k=1+((p)-1)*((N)==3); ( k<(n3)-((N)!=3) ) && ( ( k==(p) ) || (N)!=3 ); k++) \
  for(j=1+((p)-1)*((N)==2); ( j<(n2)-((N)!=2) ) && ( ( j==(p) ) || (N)!=2 ); j++) \
  for(i=1+((p)-1)*((N)==1); ( i<(n1)-((N)!=1) ) && ( ( i==(p) ) || (N)!=1 ); i++)

/* loop over planes smoothly without any jumping in i,j or k */
#define  forplane1_nojump(i,j,k, n1,n2,n3, p) \
  for(i=(p), k = 0; k < (n3); k++) \
    for(j =((n2)-1)*(k%2); j < (n2) && j >= 0; j=j+1-2*(k%2))

#define  forplane2_nojump(i,j,k, n1,n2,n3, p) \
  for(j=(p), k = 0; k < (n3); k++) \
    for(i =((n1)-1)*(k%2); i < (n1) && i >= 0; i=i+1-2*(k%2))

#define  forplane3_nojump(i,j,k, n1,n2,n3, p) \
  for(k=(p), j = 0; j < (n2); j++) \
    for(i =((n1)-1)*(j%2); i < (n1) && i >= 0; i=i+1-2*(j%2))

/* loop over bfaces in a box */
#define forallbfaces(box,fi) \
  for (fi = 0; fi < box->nbfaces; fi++)


/**************************************************************************/
/* loops over point lists */

/* loop over all points ijk in all boxes in list pointlist */
#define forPointList(pointlist, boxindx, pi , ijk) \
  forallboxes(pointlist->grid,boxindx) \
    for(pi = 0 , ijk = pointlist->point[boxindx][pi]; \
        pi < pointlist->npoints[boxindx]; \
        pi++, ijk = pointlist->point[boxindx][pi-(pi==pointlist->npoints[boxindx])])

/* loop over all points ijk in all boxes in list pointlist */
#define forPointList_inbox(pointlist, box, pi , ijk) \
  for(pi = 0 , ijk = pointlist->point[box->b][pi]; \
        pi < pointlist->npoints[box->b]; \
        pi++, ijk = pointlist->point[box->b][pi-(pi==pointlist->npoints[box->b])])
