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
/* same as AddToSparseVector, but let only one OpenMP thread in */
void AddToSparseVector_critical(tSparseVector *SV, int newpos, double newval)
{
  #pragma omp critical (AddToSparseVector_c)
  { AddToSparseVector(SV, newpos, newval); }
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

/* Allocate memory for a SparseVector array of size n */
tSparseVector **AllocateSparseVectorArray(int n)
{
  int i;
  tSparseVector **A;

  A = calloc(n, sizeof(*A));
  if(A!=NULL)
    for(i=0; i<n; i++)  A[i]=AllocateSparseVector();
  
  return A;
}
/* free the above SparseVector array of size n */
void FreeSparseVectorArray(tSparseVector **A, int n)
{
  if(A!=NULL)
  {
    int i;
    for(i=0; i<n; i++)  FreeSparseVector(A[i]);
    free(A);
  }
}

/* print a SparseVector array */
void prSparseVectorArray(tSparseVector **A, int n)
{
  int i;

  printf("SparseVectorArray containing %d sparse vectors:\n", n);
  for(i=0; i<n; i++) 
  {
    printf("sparse vector %d:\n", i);
    prSparseVector(A[i]);
  }
}

/* write a SparseVector array in Matrix Market Exchange Format:
   if isAcol=1 we assume matrix is made up of sparse column vectors
   otherwise we assume matrix is made up of sparse row vectors */
int write_SparseVectorArray_inMatrixMarketFormat(char *filename,
                                                 tSparseVector **A, int n,
                                                 int isAcol)
{
  FILE *fp;
  int i,j, nz, maxpos, ncols, nrows;

  /* open file */
  fp = fopen(filename, "a");
  if(!fp) { printf("failed opening %s\n", filename);  return -1; }

  /* count number of entries, and find max pos in sparse matrix */
  nz=0; maxpos=0;
  for(j=0; j<n; j++)
  {
    int entries=A[j]->entries;

    nz+=entries;
    for(i=0; i<entries; i++)
      if(A[j]->pos[i] > maxpos)  maxpos=A[j]->pos[i];
  }
  if(isAcol) { ncols=n; nrows=maxpos+1; }
  else       { nrows=n; ncols=maxpos+1; }

  /* write header */
  fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(fp, "%% Written by write_SparseVectorArray_inMatrixMarketFormat\n");
  fprintf(fp, "%% filename=%s\n", filename);
  fprintf(fp, "%% A=%p\n", A);
  fprintf(fp, "%% n=%d\n", n);
  fprintf(fp, "%% isAcol=%d\n", isAcol);
  fprintf(fp, "%% This is a %dx%d matrix with %d nonzeros\n", nrows,ncols, nz);
  fprintf(fp, "%% BEGIN_matrix:\n");
  fprintf(fp, "%d %d %d\n", nrows,ncols, nz);

  /* loop over vectors */
  for(j=0; j<n; j++)
  {
    int row,col;

    for(i=0; i<A[j]->entries; i++)
    {
      if(isAcol) { row=A[j]->pos[i]+1;  col=j+1; }
      else       { row=j+1;             col=A[j]->pos[i]+1; }
      fprintf(fp, "%d %d %.15g\n", row,col, A[j]->val[i]);
    }
  }
    
  /* close file */
  fclose(fp);
  
  return 0;
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
/* Note: CAN be used to compute A^T x if A is stored as columns in e.g. Acol */
void SparseMatrixLines_times_vector(tSparseVector **Aline, int nlines,
                                    double *x, double *f)
{
  int i;

  SGRID_LEVEL4_Pragma(omp parallel for)
  for(i=0; i<nlines; i++)
  {
    int j,ent;
    double Aij;
    int *Aline_i_pos = Aline[i]->pos;
    double *Aline_i_val = Aline[i]->val;
 
    f[i] = 0.0; 
    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       j   = Aline_i_pos[ent];
       Aij = Aline_i_val[ent];
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
  int k;

  for(k=0; k<nlines; k++) f[k] = 0.0;

  SGRID_LEVEL4_Pragma(omp parallel for)
  for(k=0; k<nlines; k++)
  {
    int j,ent;
    double Akj;
    int *Aline_k_pos = Aline[k]->pos;
    double *Aline_k_val = Aline[k]->val;

    for(ent = 0; ent < Aline[k]->entries; ent++)
    {
       j   = Aline_k_pos[ent];
       Akj = Aline_k_val[ent];
       f[j] += x[k] * Akj;
    }
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
  int i;
 
  for(i=0; i<nlines; i++) f[i] = 0.0;

  SGRID_LEVEL4_Pragma(omp parallel for)
  for(i=0; i<nlines; i++)
  {
    int j,ent;
    double Aij;
    int *Aline_i_pos = Aline[i]->pos;
    double *Aline_i_val = Aline[i]->val;

    for(ent = 0; ent < Aline[i]->entries; ent++)
    {
       j   = Aline_i_pos[ent];
       Aij = Aline_i_val[ent];
       f[j] += Aij * x[i];
    }
  }
}
