/* sgrid_MemoryMan_loops.h */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99, 10/02 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details.

*/



/****************************************************************************/
/* for all points without stencil:
   for all is simple because this is how we store them 
   the name of the index is passed to the macro 
*/
#define forallpoints(box,i) \
  for (i = 0; i < box->nnodes; i++)

