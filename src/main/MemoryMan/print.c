/* print.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann, 12/99, 3/03 */

#include "sgrid.h"
#include "MemoryMan.h"


/* print information */

void printgrid(tGrid *g) 
{
  int i;

  printf("grid %p: nboxes=%d, nvariables=%d\n",
	 g, g->nboxes, g->nvariables);
  for (i = 0; i < g->nboxes; i++)
  {
    printbox(g->box[i]);
    printbbox(g->box[i], g->box[i]->bbox, g->box[i]->ibbox);
  }
}




void printbox(tBox *box) 
{
  printf(" box %p: grid=%p, b=%d, n1=%d n2=%d n3=%d,\n"
         " nnodes=%d, node=%p, v=%p,\n"
         " D1 =%p D2 =%p D3 =%p,\n"
         " DD1=%p DD2=%p DD3=%p,\n"
         " F1 =%p F2 =%p F3 =%p\n",
	 box, box->grid, box->b, box->n1, box->n2, box->n3, box->nnodes, 
	 box->node, box->v, box->D1, box->D2, box->D3,
	 box->DD1, box->DD2, box->DD3,   box->F1, box->F2, box->F3);
}




void printnode(tNode *n) 
{
  printf(" node %p: i %d", n, n->i);
  printf("\n");
}




void printbbox(tBox *box, double *bbox, int *ibbox)
{
  if (dless(bbox[1], bbox[0]) ||
      dless(bbox[3], bbox[2]) ||
      dless(bbox[5], bbox[4])) {
    printf("b=%d\n",box->b);
    return;
  }

  printf(" b=%d, [%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f] ", 
	 box->b, 
	 bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
  if (ibbox)
    printf("  %dx%dx%d = %d", 
	   ibbox[1]+1, ibbox[3]+1, ibbox[5]+1, 
	   (ibbox[1]+1)*(ibbox[3]+1)*(ibbox[5]+1));
  printf("\n");
}

/* print a variable in a box */
void printvar_inbox(tBox *box, char *name)
{
  double *v = box->v[Ind(name)];
  int i,j,k;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;

  printf("%s, %p: \n", name, v);
  if(v)
  {
    forallijk(i,j,k)
    {
      printf("%g ", v[Index(i,j,k)]);
    }
  }    
  printf("%s ends here.\n", name);
}

/* print a variable */
void printvar(tGrid *grid, char *name)
{
  int b;
  
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    printf("box=%p box->b=%d: ", box, box->b);
    printvar_inbox(box, name);      
  }
}

/* print a variable list */
void printVarList(tVarList *vl)
{
  int i;
  tGrid *grid = vl->grid;

  for(i = 0; i < vl->n; i++)
    printvar(grid, VarName(vl->index[i]));      
}

