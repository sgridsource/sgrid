/* sgrid_MemoryMan_loops.h */
/* Wolfgang Tichy, April 2005 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details.
*/

/****************************************************************************/
/* loop over all points in a box */
#define forallpoints(box,i) \
  for (i = 0; i < box->nnodes; i++)

/* loop over all points in all boxes */
#define forallgridpoints(grid,box,boxindex,ijk) \
  for(box=grid->box[boxindex], boxindex=0; boxindex<grid->nboxes; boxindex++) \
    for(ijk=0; ijk<box->nnodes; ijk++)

/* loop over all boxes */
#define forallboxes(grid,boxindex) \
  for(boxindex=0; boxindex < grid->nboxes; boxindex++)
