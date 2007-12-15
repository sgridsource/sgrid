/* ReadMatrixFrom_Fx.c */
/* Wolfgang Tichy 8/2007 */


#include "sgrid.h"
#include "LinAlg.h"


/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* for serveral vars in varlist vlx                                 */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (u0)   (r0)                                        */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (u1)   (r1)                                        */
/*     ( . . . ) (v1) = (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:*                     */
/* offset[0]                            = first point in box0       */
/* offset[0] + grid->box[0]->nnodes-1   =  last point in box0       */
/* offset[1]                            = first point in box1       */
/* offset[1] + grid->box[1]->nnodes-1   =  last point in box1       */
/* offset[2]                            = first point in box2       */
/* offset[2] + grid->box[2]->nnodes-1   =  last point in box2       */
/* ...                                                              */
/* where: offset[0] = 0,                                            */
/*        offset[b] = Sum_{c=1}^{b} (grid->box[c-1]->nnodes)*nvars  */
/* I.e. the total number of lines for nvars=vlx->n variables is:    */ 
/*   Sum_{b=0}^{nboxes-1} (grid->box[b]->nnodes)*nvars              */
/* **************************************************************** */


/* set matrix lines Aline[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int line,nlines=0;
    tSparseVector **Aline;
    forallboxes(grid,b)  nlines+=(grid->box[b]->nnodes)*nvars;
    Aline=calloc(nlines, sizeof(*Aline));
    for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();
*/
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixLines_slowly(tSparseVector **Aline,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  int b, col, line;

  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  /* go over all boxes, points and vars */
  col = 0;
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i,j, bb;

    if(pr)
    {
      printf("\n"); prdivider(0);  
      printf("SetMatrixLines_slowly: working on box%d\ncol=",b);
    }
    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *x = box->v[vlx->index[j]];

        // col  = offset[b] + i*vlx->n + j;  
        // where: offset[b] = grid->box[b-1]->nnodes
        if(pr) { printf("%d ",col); fflush(stdout);}

        /* put a single 1 into x at point i, i.e. in line=col */
        x[i]=1;
        
        /* evaluate Fx if x has a single 1 in line=col, in order to
           compute matrix column col */
        Fx(vlFx, vlx, vlc1, vlc2);

        /* check where in column col there are entries */
        line = 0;
        forallboxes(grid,bb)
        {
          tBox *box = grid->box[bb];
          int i,j;

          forallpoints(box,i)
            for(j = 0; j < vlx->n; j++)
            {
              double *Fx = box->v[vlFx->index[j]]; 
              if(Fx[i]!=0)  AddToSparseVector(Aline[line], col, Fx[i]);
              line++;
            }
        } /* end: forallboxes(grid,bb) */

        /* remove the 1 in x at point i */
        x[i]=0;

        /* increase column counter */
        col++;
      }
  } /* end: forallboxes(grid,b) */
  if(pr)
  {
    printf("\nSetMatrixLines_slowly: "
           "the %d*%d matrix Aline=%p is now set!\n",
           col, col, Aline);
    fflush(stdout);
  }
}
/* NOTE: 
In order to speed up SetMatrixLines_slowly we should make an Fx
that computes vlFx only in the following planes:
+the three coord-planes containing the point where x=1 in the curreent box
+the faces of each box adjacent to the current box.
Since only in those places vlFx can be non-zero if x=1 only at one point !!! */


/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* for serveral vars in varlist vlx                                 */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (u0)   (r0)                                        */
/*     ( . . . ) (u1)   (r1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. ) = (. )                                        */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (v1)   (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:*                     */
/* offset[0]                            = first point in box0       */
/* offset[0] + grid->box[0]->nnodes-1   =  last point in box0       */
/* offset[1]                            = first point in box1       */
/* offset[1] + grid->box[1]->nnodes-1   =  last point in box1       */
/* offset[2]                            = first point in box2       */
/* offset[2] + grid->box[2]->nnodes-1   =  last point in box2       */
/* ...                                                              */
/* where: offset[0] = 0,                                            */
/*        offset[b] = Sum_{c=1}^{b} (grid->box[c-1]->nnodes)*nvars  */
/* I.e. the total number of lines for nvars=vlx->n variables is:    */ 
/*   Sum_{b=0}^{nboxes-1} (grid->box[b]->nnodes)*nvars              */
/* **************************************************************** */

/* set matrix lines Aline[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int line,nlines=0;
    tSparseVector **Aline;
    forallboxes(grid,b)  nlines+=(grid->box[b]->nnodes)*nvars;
    Aline=calloc(nlines, sizeof(*Aline));
    for(line=0; line<nlines; line++)  Aline[line]=AllocateSparseVector();
*/
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixLines_forSortedVars_slowly(tSparseVector **Aline,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  int b, col, line;

  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  /* go over all boxes, points and vars */
  col = 0;
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i,j, bb;

    if(pr)
    {
      printf("\n"); prdivider(0);  
      printf("SetMatrixLines_forSortedVars_slowly: working on box%d\ncol=",b);
    }
    /* loop over vars */
    for(j = 0; j < vlx->n; j++)
    {
      double *x = box->v[vlx->index[j]];

      forallpoints(box,i)
      {
        // col  = offset[b] + j*grid->box[b]->nnodes + i; 
        // where: offset[b] = grid->box[b-1]->nnodes
        if(pr) { printf("%d ",col); fflush(stdout);}

        /* put a single 1 into x at point i, i.e. in line=col */
        x[i]=1;
        
        /* evaluate Fx if x has a single 1 in line=col, in order to
           compute matrix column col */
        Fx(vlFx, vlx, vlc1, vlc2);

        /* check where in column col there are entries */
        line = 0;
        forallboxes(grid,bb)
        {
          tBox *box = grid->box[bb];
          int i,j;

          for(j = 0; j < vlx->n; j++)
          {
            double *Fx = box->v[vlFx->index[j]]; 
            forallpoints(box,i)
            {
              if(Fx[i]!=0)  AddToSparseVector(Aline[line], col, Fx[i]);
              line++;
            }
          }
        } /* end: forallboxes(grid,bb) */

        /* remove the 1 in x at point i */
        x[i]=0;

        /* increase column counter */
        col++;
      }
    } /* end: loop over vars */
  } /* end: forallboxes(grid,b) */
  if(pr)
  {
    printf("\nSetMatrixLines_forSortedVars_slowly: "
           "the %d*%d matrix Aline=%p is now set!\n",
           col, col, Aline);
    fflush(stdout);
  }
}


/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* for serveral vars in varlist vlx                                 */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (u0)   (r0)                                        */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (u1)   (r1)                                        */
/*     ( . . . ) (v1) = (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:*                     */
/* offset[0]                            = first point in box0       */
/* offset[0] + grid->box[0]->nnodes-1   =  last point in box0       */
/* offset[1]                            = first point in box1       */
/* offset[1] + grid->box[1]->nnodes-1   =  last point in box1       */
/* offset[2]                            = first point in box2       */
/* offset[2] + grid->box[2]->nnodes-1   =  last point in box2       */
/* ...                                                              */
/* where: offset[0] = 0,                                            */
/*        offset[b] = Sum_{c=1}^{b} (grid->box[c-1]->nnodes)*nvars  */
/* I.e. the total number of lines for nvars=vlx->n variables is:    */ 
/*   Sum_{b=0}^{nboxes-1} (grid->box[b]->nnodes)*nvars              */
/* **************************************************************** */

/* set matrix columns Acol[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int col,ncols=0;
    tSparseVector **Acol;
    forallboxes(grid,b)  ncols+=(grid->box[b]->nnodes)*nvars;
    Acol=calloc(ncols, sizeof(*Acol));
    for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();
*/
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixColumns_slowly(tSparseVector **Acol,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  int b, col, line;

  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  /* go over all boxes, points and vars */
  col = 0;
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i,j, bb;

    if(pr)
    {
      printf("\n"); prdivider(0);  
      printf("SetMatrixColumns_slowly: working on box%d\ncol=",b);
    }
    forallpoints(box,i)
      for(j = 0; j < vlx->n; j++)
      {
        double *x = box->v[vlx->index[j]];

        // col  = offset[b] + i*vlx->n + j;  
        // where: offset[b] = grid->box[b-1]->nnodes
        if(pr) { printf("%d ",col); fflush(stdout);}

        /* put a single 1 into x at point i, i.e. in line=col */
        x[i]=1;
        
        /* evaluate Fx if x has a single 1 in line=col, in order to
           compute matrix column col */
        Fx(vlFx, vlx, vlc1, vlc2);

        /* check where in column col there are entries */
        line = 0;
        forallboxes(grid,bb)
        {
          tBox *box = grid->box[bb];
          int i,j;

          forallpoints(box,i)
            for(j = 0; j < vlx->n; j++)
            {
              double *Fx = box->v[vlFx->index[j]]; 
              if(Fx[i]!=0)  AddToSparseVector(Acol[col], line, Fx[i]);
              line++;
            }
        } /* end: forallboxes(grid,bb) */

        /* remove the 1 in x at point i */
        x[i]=0;

        /* increase column counter */
        col++;
      }
  } /* end: forallboxes(grid,b) */
  if(pr)
  {
    printf("\nSetMatrixColumns_slowly: "
           "the %d*%d matrix Acol=%p is now set!\n",
           col, col, Acol);
    fflush(stdout);
  }
}


/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* for serveral vars in varlist vlx                                 */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (u0)   (r0)                                        */
/*     ( . . . ) (u1)   (r1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. ) = (. )                                        */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (v1)   (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:*                     */
/* offset[0]                            = first point in box0       */
/* offset[0] + grid->box[0]->nnodes-1   =  last point in box0       */
/* offset[1]                            = first point in box1       */
/* offset[1] + grid->box[1]->nnodes-1   =  last point in box1       */
/* offset[2]                            = first point in box2       */
/* offset[2] + grid->box[2]->nnodes-1   =  last point in box2       */
/* ...                                                              */
/* where: offset[0] = 0,                                            */
/*        offset[b] = Sum_{c=1}^{b} (grid->box[c-1]->nnodes)*nvars  */
/* I.e. the total number of lines for nvars=vlx->n variables is:    */ 
/*   Sum_{b=0}^{nboxes-1} (grid->box[b]->nnodes)*nvars              */
/* **************************************************************** */

/* set matrix columns Acol[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int col,ncols=0;
    tSparseVector **Acol;
    forallboxes(grid,b)  ncols+=(grid->box[b]->nnodes)*nvars;
    Acol=calloc(ncols, sizeof(*Acol));
    for(col=0; col<ncols; col++)  Acol[col]=AllocateSparseVector();
*/
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixColumns_forSortedVars_slowly(tSparseVector **Acol,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  int b, col, line;

  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  /* go over all boxes, points and vars */
  col = 0;
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i,j, bb;

    if(pr)
    {
      printf("\n"); prdivider(0);  
      printf("SetMatrixColumns_forSortedVars_slowly: working on box%d\ncol=",b);
    }
    for(j = 0; j < vlx->n; j++)
    {
      double *x = box->v[vlx->index[j]];

      forallpoints(box,i)
      {
        // col  = offset[b] + i*vlx->n + j;  
        // where: offset[b] = grid->box[b-1]->nnodes
        if(pr) { printf("%d ",col); fflush(stdout);}

        /* put a single 1 into x at point i, i.e. in line=col */
        x[i]=1;
        
        /* evaluate Fx if x has a single 1 in line=col, in order to
           compute matrix column col */
        Fx(vlFx, vlx, vlc1, vlc2);

        /* check where in column col there are entries */
        line = 0;
        forallboxes(grid,bb)
        {
          tBox *box = grid->box[bb];
          int i,j;

          for(j = 0; j < vlx->n; j++)
          {
            double *Fx = box->v[vlFx->index[j]]; 

            forallpoints(box,i)
            {
              if(Fx[i]!=0)  AddToSparseVector(Acol[col], line, Fx[i]);
              line++;
            }
          }
        } /* end: forallboxes(grid,bb) */

        /* remove the 1 in x at point i */
        x[i]=0;

        /* increase column counter */
        col++;
      }
    }
  } /* end: forallboxes(grid,b) */
  if(pr)
  {
    printf("\nSetMatrixColumns_forSortedVars_slowly: "
           "the %d*%d matrix Acol=%p is now set!\n",
           col, col, Acol);
    fflush(stdout);
  }
}
