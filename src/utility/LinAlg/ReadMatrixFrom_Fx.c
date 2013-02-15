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

/* FIXME: make parallel
   see SetMatrixColumns_slowly on how to do it. */
printf("\nFIXME: SetMatrixLines_slowly is not OpenMP parallel!!!");

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

/* FIXME: make parallel
   see SetMatrixColumns_slowly on how to do it. */
printf("\nFIXME: SetMatrixLines_forSortedVars_slowly is not OpenMP parallel!!!");

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
  int b;
  int *offset;

  /* set offsets */
  offset = (int *) calloc( grid->nboxes, sizeof(int) );
  offset[0] = 0;
  for(b=1; b < grid->nboxes; b++)
    offset[b] = offset[b-1] + (grid->box[b-1]->nnodes) * (vlx->n); 

  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  /* go over all boxes, points and vars */
  forallboxes(grid,b)
  {
    if(pr)
    {
      printf("\n"); prdivider(0);  
      printf("SetMatrixColumns_slowly: working on box%d\ncol=",b);
    }

    SGRID_LEVEL6_Pragma(omp parallel)
    {
      int i;
      tGrid *grid_p = make_empty_grid(grid->nvariables, 0);
      tBox *box_p = grid_p->box[b];
      tVarList *vlFx_p = vlalloc(grid_p);
      tVarList *vlx_p  = vlalloc(grid_p);
      tVarList *vlc1_p = vlalloc(grid_p);
      tVarList *vlc2_p = vlalloc(grid_p);

      /* make local copy of grid and var lists */      
      copy_grid(grid, grid_p, 0);
      vlpushvl(vlFx_p, vlFx);
      vlpushvl(vlx_p, vlx);
      vlpushvl(vlc1_p, vlc1);
      vlpushvl(vlc2_p, vlc2);

      SGRID_LEVEL6_Pragma(omp for)
      forallpoints(box_p,i)
      {
        int j;

        for(j = 0; j < vlx_p->n; j++)
        {
          double *x = box_p->v[vlx_p->index[j]];
          int col   = offset[b] + i*(vlx_p->n) + j;
          int line;
          int bb;

          if(pr) { printf("%d ",col); fflush(stdout);}

          /* put a single 1 into x at point i, i.e. in line=col */
          x[i]=1;
          
          /* evaluate Fx if x has a single 1 in line=col, in order to
             compute matrix column col */
          Fx(vlFx_p, vlx_p, vlc1_p, vlc2_p);

          /* check where in column col there are entries */
          line = 0;
          forallboxes(grid_p,bb)
          {
            tBox *box = grid_p->box[bb];
            int i,j;

            forallpoints(box,i)
              for(j = 0; j < vlx_p->n; j++)
              {
                double *Fx = box->v[vlFx_p->index[j]]; 
                if(Fx[i]!=0)  AddToSparseVector(Acol[col], line, Fx[i]);
                line++;
              }
          } /* end: forallboxes(grid_p,bb) */

          /* remove the 1 in x at point i */
          x[i]=0;
        } /* end: for(j = 0; j < vlx_p->n; j++) */
      }
      /* free local copies */
      vlfree(vlFx_p);
      vlfree(vlx_p);
      vlfree(vlc1_p);
      vlfree(vlc2_p);
      free_grid(grid_p);
    }
  } /* end: forallboxes(grid0,b) */
  if(pr)
  {
    int nboxes = grid->nboxes;
    int lb_nnodes = grid->box[nboxes-1]->nnodes;
    int ncol = offset[nboxes-1] + (lb_nnodes-1)*vlx->n + vlx->n;
    printf("\nSetMatrixColumns_slowly: "
           "the %d*%d matrix Acol=%p is now set!\n",
           ncol, ncol, Acol);
    fflush(stdout);
  }
  free(offset);
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

/* FIXME: make parallel
   see SetMatrixColumns_slowly on how to do it. */
printf("\nFIXME: SetMatrixColumns_forSortedVars_slowly is not OpenMP parallel!!!");

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


/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* in one box for one var in varlist vlx                            */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (v1) = (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:                      */
/* 0              = first point in box                              */
/* box->nnodes-1  =  last point in box                              */
/* **************************************************************** */

/* set matrix columns Acol[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int ncols;
    tSparseVector **Acol;
    ncols=(grid->box->nnodes);
    Acol=AllocateSparseVectorArray(ncols);
  */
/* The index vlind of the var in vlx and the boxindex b are passed in. */
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixColumns_ForOneVarInOneBox_slowly(tSparseVector **Acol,
    int vlind, int b,
    void  (*Fx)(tVarList *Fdx,  tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  
  /* set x to zero */
  vladd(vlx, 0.0,NULL, 0.0,NULL);

  if(pr)
  {
    printf("\n"); prdivider(0);  
    printf("SetMatrixColumns_ForOneVarInOneBox_slowly: var%d in box%d\ncol=",vlind,b);
  }

  SGRID_LEVEL6_Pragma(omp parallel)
  {
    int i;
    tGrid *grid_p = make_empty_grid(grid->nvariables, 0);
    tBox *box_p = grid_p->box[b];
    tVarList *vlFx_p = vlalloc(grid_p);
    tVarList *vlx_p  = vlalloc(grid_p);
    tVarList *vlc1_p = vlalloc(grid_p);
    tVarList *vlc2_p = vlalloc(grid_p);

    /* make local copy of grid and var lists */      
    copy_grid(grid, grid_p, 0);
    vlpushvl(vlFx_p, vlFx);
    vlpushvl(vlx_p, vlx);
    vlpushvl(vlc1_p, vlc1);
    vlpushvl(vlc2_p, vlc2);

    SGRID_LEVEL6_Pragma(omp for)
    forallpoints(box_p,i)
    {
      double *x = box_p->v[vlx_p->index[vlind]];
      int col   = i;
      int line;

      if(pr) { printf("%d ",col); fflush(stdout);}

      /* put a single 1 into x at point i */
      x[i]=1;
      
      /* evaluate Fx if x has a single 1 in line=i, in order to
         compute matrix column col */
      Fx(vlFx_p, vlx_p, vlc1_p, vlc2_p);

      /* check where in column col there are entries */
      forallpoints(box_p, line)
      {
        double *Fx = box_p->v[vlFx_p->index[vlind]]; 
        if(Fx[line]!=0)  AddToSparseVector(Acol[col], line, Fx[line]);
      }

      /* remove the 1 in x at point i */
      x[i]=0;
    }
    /* free local copies */
    vlfree(vlFx_p);
    vlfree(vlx_p);
    vlfree(vlc1_p);
    vlfree(vlc2_p);
    free_grid(grid_p);
  }
  if(pr)
  {
    int ncol = grid->box[b]->nnodes;
    printf("\nSetMatrixColumns_ForOneVarInOneBox_slowly: "
           "the %d*%d matrix Acol=%p is now set!\n",
           ncol, ncol, Acol);
    fflush(stdout);
  }
}

/* **************************************************************** */
/* solve linear Matrix-Equations A x = b  			    */
/* in one box for one var in varlist vlx                            */
/* The system of equantions has the form:                           */
/*         A      x      b                                          */
/*     ( . . . ) (v0)   (s0)                                        */
/*     ( . . . ) (v1) = (s1)                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/*     ( . . . ) (. )   (. )                                        */
/* where vlx is comprised of the vars u and v and b is comprised of */
/* of r and s, and the index 0, 1, ... denotes the grid point.      */
/* Note the grid points are labeled like this:                      */
/* 0              = first point in box                              */
/* box->nnodes-1  =  last point in box                              */
/* **************************************************************** */

/* set matrix columns Acol[l] of type (tSparseVector *)  */
/* before we call this routine we need to do this:
    int ncols;
    tSparseVector **Acol;
    ncols=(grid->box->nnodes);
    Acol=AllocateSparseVectorArray(ncols);
  */
/* The index vlind of the var in vlx and the boxindex b are passed in. */
/* It's very slow!!! since we loop over all points and Fx has to 
   be called for each of these points!!! */
void SetMatrixColumns_ForOneVarInOneSubBox_slowly(tSparseVector **Acol,
    int vlind, int b, 
    int sbi, int sbj, int sbk,  int nsb1, int nsb2, int nsb3,
    void  (*Fx)(tVarList *Fdx, tVarList *dx,  tVarList *c1, tVarList *c2),
    tVarList *vlFx, tVarList *vlx, tVarList *vlc1, tVarList *vlc2, int pr)
{
  tGrid *grid = vlx->grid;
  tBox *box = grid->box[b];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i1,i2, j1,j2, k1,k2;

  /* set i1,i2, j1,j2, k1,k2 */
  IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
  //if(pr) printf("\ni: %d...%d  j: %d...%d  k:%d...%d\n", i1,i2-1, j1,j2-1, k1,k2-1);

  if(pr)
  {
    prdivider(0);  
    printf("SetMatrixColumns_ForOneVarInOneSubBox_slowly: var%d, box%d, "
           "subbox %d/%d,%d/%d,%d/%d\ncol=", 
           vlind,b, sbi,nsb1-1, sbj,nsb2-1, sbk,nsb3-1);
  }

  SGRID_LEVEL6_Pragma(omp parallel)
  {
    tGrid *grid_p = make_empty_grid(grid->nvariables, 0);
    tBox *box_p = grid_p->box[b];
    tVarList *vlFx_p = vlalloc(grid_p);
    tVarList *vlx_p  = vlalloc(grid_p);
    tVarList *vlc1_p = vlalloc(grid_p);
    tVarList *vlc2_p = vlalloc(grid_p);
    tVarBoxSubboxIndices blockinfo[1];
    int i,j,k;

    /* make local copy of grid and var lists */      
    copy_grid(grid, grid_p, 0);
    vlpushvl(vlFx_p, vlFx);
    vlpushvl(vlx_p, vlx);
    vlpushvl(vlc1_p, vlc1);
    vlpushvl(vlc2_p, vlc2);

    /* put block info into vlx_p: can be used by func Fx called below to
       loop only over relevant block */
    blockinfo->vli = vlind;
    blockinfo->bi  = b;
    blockinfo->sbi = sbi;
    blockinfo->sbj = sbj;
    blockinfo->sbk = sbk;
    blockinfo->nsb1 = nsb1;
    blockinfo->nsb2 = nsb2;
    blockinfo->nsb3 = nsb3;
    blockinfo->vari = vlx_p->index[vlind];
    vlx_p->vlPars = (void *) blockinfo;

    /* set x to zero, and call Fx to initialize all to zero in Fx */
    vlsetconstant(vlx_p, 0.0);
    Fx(vlFx_p, vlx_p, vlc1_p, vlc2_p);

    /* loop over subbox */
    SGRID_LEVEL6_Pragma(omp for)
    for(i=i1; i<i2; i++) /* do i loop first since usually i2-i1 > k2-k1 */
    for(k=k1; k<k2; k++)
    for(j=j1; j<j2; j++)
    {
      double *x = box_p->v[vlx_p->index[vlind]];
      int ijk = Index(i,j,k); 
      int col = Ind_n1n2( (i-i1),(j-j1),(k-k1), (i2-i1),(j2-j1) );
      int ii,jj,kk;

      //if(pr) { printf("ijk=%d ", ijk); }
      if(pr) { printf("%d ", col); fflush(stdout);}

      /* put a single 1 into x at point ijk */
      x[ijk]=1;
      
      /* evaluate Fx if x has a single 1 in line=i, in order to
         compute matrix column col */
      Fx(vlFx_p, vlx_p, vlc1_p, vlc2_p);

      /* check where in column col there are entries */
      for(kk=k1; kk<k2; kk++)
      for(jj=j1; jj<j2; jj++)
      for(ii=i1; ii<i2; ii++)
      {
        double *Fx = box_p->v[vlFx_p->index[vlind]];
        int iijjkk = Index(ii,jj,kk);
        int line = Ind_n1n2( (ii-i1),(jj-j1),(kk-k1), (i2-i1),(j2-j1) );
        if(Fx[iijjkk]!=0)  AddToSparseVector(Acol[col], line, Fx[iijjkk]);
      }

      /* remove the 1 in x at point ijk */
      x[ijk]=0;
    }
    /* free local copies */
    vlfree(vlFx_p);
    vlfree(vlx_p);
    vlfree(vlc1_p);
    vlfree(vlc2_p);
    free_grid(grid_p);
  }

  if(pr)
  {
    int ncols = (i2-i1)*(j2-j1)*(k2-k1);
    printf("\nSetMatrixColumns_ForOneVarInOneSubBox_slowly: "
           "the %d*%d matrix Acol=%p is now set!\n", ncols,ncols, Acol);
    prdivider(0);
    fflush(stdout);
  }
}
