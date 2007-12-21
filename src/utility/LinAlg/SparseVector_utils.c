/* hypre_sgrid_utils.c */
/* Wolfgang Tichy 8/2003 */


#include "sgrid.h"
#include "LinAlg.h"

/* ******************************************************************** */
/* WT: Here are some utility routines to use SparseVectors of           */
/*     type tSparseVector *, defined in sgrid_LinAlg.h                  */
/* ******************************************************************** */

/* Allocate memory for a SparseVector */
tSparseVector *AllocateSparseVector(void)
{
 tSparseVector *SV;
 
 SV=calloc(1, sizeof(*SV) );

 SV->entries=0;
 SV->pos=NULL;
 SV->val=NULL;
 
 return SV;
}

/* add one pos to a SparseVector */
void AddToSparseVector(tSparseVector *SV, int newpos, double newval)
{
 void *ret;
 
 ret=realloc(SV->pos, (sizeof( *(SV->pos) ))*(SV->entries+1) );
 if(ret==NULL) return;
 SV->pos=ret;
 SV->pos[SV->entries]=newpos;

 ret=realloc(SV->val, (sizeof( *(SV->val) ))*(SV->entries+1) );
 if(ret==NULL) return;
 SV->val=ret;
 SV->val[SV->entries]=newval;
  
 SV->entries++;
}

/* free a SparseVector */
void FreeSparseVector(tSparseVector *SV)
{
 if(SV!=NULL)
 {
   free(SV->pos);
   free(SV->val);
   free(SV);
 }
}

/* print a SparseVector */
void prSparseVector(tSparseVector *SV)
{
  int i;
  
  printf("entry\tpos\tval\tSparseVector %p  entries=%d\n",
         SV , SV->entries);
  for(i=0; i<SV->entries; i++)
    printf("%d\t%d\t%.16e\n",i, SV->pos[i], SV->val[i]);
}

/* get component comp  */
double GetSparseVectorComponent(tSparseVector *SV, int comp)
{
  int i;
  
  for(i=0; i<SV->entries; i++)
    if(SV->pos[i]==comp) return SV->val[i];
  return 0.0;
}

/* matrix times vector for a matrix Aline made up of lines of tSparseVector
   vectors. Initialize like Aline this: 
    tSparseVector **Aline;
    Aline = calloc(nlines, sizeof(*Aline));
    AddToSparseVector(Aline[line], col, value);  */
void SparseMatrixLines_times_vector(tSparseVector **Aline, int nlines,
                                    double *x, double *f)
{
  int i,j,ent;
  double Aij;

  for(i=0; i<nlines; i++)
  {
    f[i] = 0.0; 
    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       j   = Aline[i]->pos[ent];
       Aij = Aline[i]->val[ent];
       f[i] += Aij * x[j];
    }
  }
}

/* transpose vector times a matrix Aline made up of lines of tSparseVector
   vectors. Initialize like Aline this: 
    tSparseVector **Aline;
    Aline = calloc(nlines, sizeof(*Aline));
    AddToSparseVector(Aline[line], col, value);  */
void vector_times_SparseMatrixLines(double *x, tSparseVector **Aline, 
                                    int nlines, double *f)
{
  int j,k,ent;
  double Akj;

  for(k=0; k<nlines; k++) f[k] = 0.0;
  for(k=0; k<nlines; k++)
    for(ent = 0; ent < Aline[k]->entries; ent++)
    {
       j   = Aline[k]->pos[ent];
       Akj = Aline[k]->val[ent];
       f[j] += x[k] * Akj;
    }
}

/* matrix times vector for the transpose of Aline.  Aline is made up 
   of lines of tSparseVector vectors. Initialize like Aline this: 
    tSparseVector **Aline;
    Aline = calloc(nlines, sizeof(*Aline));
    AddToSparseVector(Aline[line], col, value);  */
void SparseMatrixLinesTranspose_times_vector(tSparseVector **Aline, int nlines,
                                             double *x, double *f)
{
  int i,j,ent;
  double Aij;
 
  for(i=0; i<nlines; i++) f[i] = 0.0;
  for(i=0; i<nlines; i++)
    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       j   = Aline[i]->pos[ent];
       Aij = Aline[i]->val[ent];
       f[j] += Aij * x[i];
    }
}
