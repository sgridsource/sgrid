/* MemoryMan.h */
/* Wolfgang Tichy, April 2005 */


/* storage.c */
tNode *alloc_nodes(tBox *box, int nnodes);
tBox *alloc_box(tGrid *g, int b, int n1, int n2, int n3); 
void realloc_boxvariables(tBox *box, int nvariables);
void realloc_gridvariables(tGrid *grid, int nvariables);
