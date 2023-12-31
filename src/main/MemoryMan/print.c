/* print.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "MemoryMan.h"


/* print information */

void printgrid(tGrid *g) 
{
  int i;

  printf("grid=%p: nboxes=%d, nvariables=%d, dt=%g\n",
	 g, g->nboxes, g->nvariables, g->dt);
  for (i = 0; i < g->nboxes; i++)
    printbox(g->box[i]);
}

void printbox(tBox *box) 
{
  printf("b=%d: box=%p, grid=%p, n1=%d n2=%d n3=%d,"
         " nnodes=%d, node=%p, v=%p\n",
	 box->b, box, box->grid, box->n1, box->n2, box->n3, 
	 box->nnodes, 
	 box->node, box->v);
  printbbox(box, box->bbox, box->ibbox);
  /*
  printf(" D1 =%p D2 =%p D3 =%p,\n"
         " DD1=%p DD2=%p DD3=%p,\n"
         " F1 =%p F2 =%p F3 =%p,\n"
         " Mcoeffs1=%p Mcoeffs2=%p Mcoeffs3=%p,\n"
         " Meval1  =%p Meval2  =%p Meval3  =%p,\n"
         " Int1    =%p Int2    =%p Int3    =%p\n",
	 box->D1, box->D2, box->D3,
	 box->DD1, box->DD2, box->DD3,  box->F1, box->F2, box->F3,
	 box->Mcoeffs1, box->Mcoeffs2, box->Mcoeffs3,
	 box->Meval1, box->Meval2, box->Meval3,
	 box->Int1, box->Int2, box->Int3);
  */
  printf(" TransformType1=%d TransformType2=%d TransformType3=%d\n",
         box->TransformType1, box->TransformType2, box->TransformType3);
}

void printnode(tNode *n) 
{
  printf(" node %p: i %d", n, n->i);
  printf("\n");
}

void printCI(tBox *box)
{
  tCoordInfo *CI = box->CI;
  int i, j, p;

  printf("b%d: CI->dom=%d CI->type=%d CI->useF=%d\n",
         box->b, CI->dom, CI->type, CI->useF);

  /* s, xc */
  printf("CI->s[] = [ ");
  for(i=0;i<6;i++) printf("%g ", CI->s[i]);
  printf("]\n");
  printf("CI->xc[] = [ ");
  for(i=0;i<4;i++) printf("%g ", CI->xc[i]);
  printf("]\n");

  /* iFS, iSurf */
  for(i=p=0; i<6; i++) if(CI->iFS[i]) p++;
  if(p)
  {
    printf("CI->iFS[] = [ ");
    for(i=0;i<6;i++) printf("%d ", CI->iFS[i]);
    printf("]\n");
    printf("CI->iSurf[] = [ ");
    for(i=0;i<6;i++) printf("%d ", CI->iSurf[i]);
    printf("]\n");
    printf("idSurfdX[][] = [");
    for(i=0;i<6;i++)
    {
      for(j=1;j<4;j++) printf(" %d", CI->idSurfdX[i][j]);
      printf(";");
    }
    printf("]\n");
  }

  /* FSurf */
  if(CI->useF)
  {
    printf("CI->FSurf[] = [ ");
    for(i=0;i<6;i++) printf("%p ", CI->FSurf[i]);
    printf("]\n");
    for(i=p=0; i<6; i++) for(j=1;j<4;j++) if(CI->dFSurfdX[i][j]) p++;
    if(p)
    {

      printf("dFSurfdX[][] = [");
      for(i=0;i<6;i++)
      {
        for(j=1;j<4;j++) printf(" %p", CI->dFSurfdX[i][j]);
        printf(";");
      }
      printf("]\n");
    }
  }
}

void printbface(tBface *bface)
{
  int np;
  printf("b%d f%d fi%d  ob%d ofi%d  oXi=%d oYi=%d oZi=%d  ",
          bface->b,bface->f,bface->fi, bface->ob,bface->ofi,
          bface->oXi,bface->oYi,bface->oZi);
  printf("bits=%u%u%u%u%u%u,%u,%u%u",
         bface->touch,
         bface->sameX, bface->sameY, bface->sameZ,  bface->same_fpts,
         bface->fpts_off_face,
         bface->setnormalderiv, bface->innerbound,  bface->outerbound);
  if(bface->fpts==NULL) np = 0;
  else                  np = bface->fpts->npoints[bface->b];
  printf("  (%d points)\n", np);
  //prPointList(bface->fpts);
  prPointList_ijk_inbox(bface->fpts, bface->b);
}

void printbfaces(tBox *box)
{
  int i;
  printf("box->b=%d\n", box->b);
  for(i=0; i<box->nbfaces; i++)
  {
    tBface *bface = box->bface[i];
    printf("bface[%d]: ", i);
    printbface(bface);
  }
}

void printmatrix(double *M, int n)
{
  int i,j;
  
  for(i=0; i<n; i++)
  {
    printf("M_{%d j}= ", i);
    for(j=0; j<n; j++) printf("%g ", M[n*i + j]);
    printf("\n");
  }
}

void printbbox(tBox *box, double *bbox, int *ibbox)
{
  if (dless(bbox[1], bbox[0]) ||
      dless(bbox[3], bbox[2]) ||
      dless(bbox[5], bbox[4])) {
    printf("b=%d\n",box->b);
    return;
  }

  printf(" X in [%.16g,%.16g],   b=%d\n", bbox[0], bbox[1], box->b);
  printf(" Y in [%.16g,%.16g]\n", bbox[2], bbox[3]);
  printf(" Z in [%.16g,%.16g]", bbox[4], bbox[5]);
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

  printf("%s, %p, Ind=%d:\n", name, v, Ind(name));
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

