/* populate_bfaces_AR.c */
/* Alireza Rashti 2/2018 */

#include "sgrid.h"
#include "Coordinates.h"
#include "populate_bfaces_AR.h"

#define EPS 1E-5
#define EPS2 1E-3


/*populating bface structure*/
int populate_bfaces(tGrid *grid)
{
  struct FACE_POINT_S ***FacePoint; // Format is FacePoint[box][face]->
  char *maxits = strdup(Gets("Coordinates_newtMAXITS")); /* save par */
  char *guess = strdup(Gets("Coordinates_XYZ_of_xyz_Guess")); /* save par */
  int b;

  /* free all bfaces that someone else may have made before */
  forallboxes(grid, b)
    free_all_bfaces(grid->box[b]);

  /* Visualize boxes */
  if (1)
    visualize_boxes(grid);

  /* Operation */
  printf("\n***Populating Bfaces***\n");

  /* Reduce iteration number in newton_linesrch_itsP so that it fails
     faster if we look for X,Y,Z in at box that does not contain it.
     Also make it pick its own initial guess (which is very slow) */
  Seti("Coordinates_newtMAXITS", 75);
  Sets("Coordinates_XYZ_of_xyz_Guess","yes");

  /*Allocating memory for face point struct*/
  FacePoint = allc_FacePoint(grid);

  /*initializing the face point struct*/
  init_FacePoint(FacePoint,grid);

  /*Finding the adjacent structure for inner points*/
  find_adjacent_inner(FacePoint,grid);

  /*Finding thr adjacent structure for edge points*/
  find_adjacent_edge(FacePoint,grid);

  /*Grouping the similar points on each face to make bfaces*/
  group_similar_points(FacePoint,grid);

  /*freeing memory*/
  free_FacePoint(FacePoint);

  /* Testing bfaces */
  if (1)
    test_bfaces(grid);

  /* Visualize bfaces*/
  if (1)
    visualize_bfaces(grid);

  /* restore Coordinates_... pars to saved values */
  Sets("Coordinates_newtMAXITS", maxits);
  Sets("Coordinates_XYZ_of_xyz_Guess", guess);
  free(guess);
  free(maxits);

  return 0;
}

/*freeing memory*/
static void free_FacePoint(struct FACE_POINT_S ***FacePoint)
{
  int b;

  for (b = 0; FacePoint[b] != NULL; b++)
  {
    int f;

    for (f = 0; FacePoint[b][f] != NULL; f++)
    {
      if (FacePoint[b][f]->s > 0)
        free(FacePoint[b][f]->inner);

      if (FacePoint[b][f]->l > 0)
        free(FacePoint[b][f]->edge);

      if (FacePoint[b][f]->sh > 0)
      {
        int i;
        for (i = 0; i < FacePoint[b][f]->sh; i++)
        {
          assert(FacePoint[b][f]->shared[i].np);

          free(FacePoint[b][f]->shared[i].indx);
          free(FacePoint[b][f]->shared[i].kind);
        }
        free(FacePoint[b][f]->shared);
      }

      free(FacePoint[b][f]);
    }

    free(FacePoint[b]);

  }

  free(FacePoint);

}

/*Grouping the similar points on each face to make bfaces*/
static void group_similar_points(struct FACE_POINT_S ***const FacePoint,tGrid *grid)
{
  int b,b_max;

  b_max = tot_num(FacePoint);

  for (b = 0; b < b_max; b++)
  {
    int n1 = grid->box[b]->n1;
    int n2 = grid->box[b]->n2;
    int n3 = grid->box[b]->n3;
    int f, f_max;
    f_max = tot_num(FacePoint[b]);

    for (f = 0; f < f_max; f++)
    {
      FLAG_T flg;

      if (FacePoint[b][f]->internal_face == 1)
        continue;

      /*find outerbound bfaces*/
      flg = find_outerbound_bfaces(FacePoint,b,f,grid);

      if (flg == CONTINUE_F) continue;

      /*finding the remaining bfeces*/
      find_remaining_bfaces(FacePoint,b,f,grid);

    }//for (f = 0; f < f_max; f++)
  }//for (b = 0; b < b_max; b++)

  /*there are some flags which are not set yet,
  furtheremore, the indices of fpts needs to be orderly*/
  setting_remaining_flags(grid);

}


/* Setting the flags */
static void setting_remaining_flags(tGrid *grid)
{

  /* set the correct ofi index if it is needed */
  set_ofi_flag(grid);

  /* Finding all of the bfaces which touching each other and need copy,
     they are called pair, and then order the index of each ftps point */
  order_ftps_pair(grid);

  /* Setting the SameX/Y/Z flags */
  set_same_X_Y_Z_flg(grid);

  /* Setting oXi, oYi, oZi flags */
  set_oXi_oYi_oZi_flg(grid);

  /* make sure flags are consistent */
  set_consistent_flags_in_all_bfaces(grid);
}

/* Setting oXi, oYi, oZi flags */
static void set_oXi_oYi_oZi_flg(tGrid *grid)
{
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      tBface *bface = box->bface[bf];

      if( (bface->same_fpts == 0) && (bface->ob >= 0) )
      {
        bface->oXi = Ind("oX");
        bface->oYi = Ind("oY");
        bface->oZi = Ind("oZ");
      }
    }
  }

}

/*Setting the SameX/Y/Z flags*/
static void set_same_X_Y_Z_flg(tGrid *grid)
{
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      tBface *bface = box->bface[bf];
      
      if (bface->touch == 1 && bface->same_fpts == 0)
      {
        tBox *box2 = grid->box[bface->ob];
        tBface *bface2 = box2->bface[bface->ofi];

        /* Given the plane and bface,
            check if there is any X or Y or Z constant */
        bface->sameX = check_sameXYZ(SAME_X,bface2);
        bface->sameY = check_sameXYZ(SAME_Y,bface2);
        bface->sameZ = check_sameXYZ(SAME_Z,bface2);
      }
    }
  }//forallboxes(grid,b)

}


/* Given the plane and bface, check if there is any X or Y or Z constant */
static unsigned int check_sameXYZ(FLAG_T kind,tBface *bface)
{
  tBox *box;
  double X0[3];
  int b, m, i,slice;

  assert(bface->fpts != 0);

  b = bface->fpts->blist[0];
  m = bface->fpts->npoints[b];
  box = bface->fpts->grid->box[b];

  get_X_coord(X0,box,bface->fpts->point[b][0]);

  if (kind == SAME_X)
    slice = 0;

  else if (kind == SAME_Y)
    slice = 1;

  else if (kind == SAME_Z)
    slice = 2;

  for (i = 1; i < m ; i++)
  {
    double X[3];
    get_X_coord(X,box,bface->fpts->point[b][i]);
    
    if (!dequal(X[slice],X0[slice]))
      return 0;
  }

  return 1;
}

/*Finding all of the bfaces which
    touching each other and need copy, they are called pair, and
    then order the index of each ftps point*/
static void order_ftps_pair(tGrid *grid)
{
  /*those bafaces which are paired*/
  struct PAIR_S *pair = 0;
  tBox *box;
  FLAG_T flg;
  int b,np;

  np = 0;
  forallboxes(grid,b)
  {
    box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      if (box->bface[bf]->touch == 1 && box->bface[bf]->same_fpts == 1)
      {
        /*check if this bface is not considered previously*/
         flg = check_bface(pair,np,box->bface[bf]);

         if (flg == CONTINUE_F) continue;

         assert(box->bface[bf]->ofi >= 0);

         add_to_pair(&pair,box->bface[bf],grid->box[box->bface[bf]->ob]->bface[box->bface[bf]->ofi],&np);
      }
    }
  }

  /*Analyzing each pair*/
  int i;
  for (i = 0; i < np; i++)
  {
    /*arranging the point index in fpts*/
    assert(pair[i].bface1->fpts->nblist == pair[i].bface2->fpts->nblist);

    /*For each blist temproray take it 0*/
    int bln1 = pair[i].bface1->fpts->blist[0];
    int bln2 = pair[i].bface2->fpts->blist[0];
    int m1 = pair[i].bface1->fpts->npoints[bln1];
    int m2 = pair[i].bface2->fpts->npoints[bln2];
    int j1;


    if (m1 != m2)
    {
      errorexit("The number of points for two copying bfaces "
            "are not matched!\n");
    }

    for (j1 = 0; j1 < m1; j1++)
    {
      int ijk1 = pair[i].bface1->fpts->point[bln1][j1];
      double x1[3];
      int j2;

      get_x_coord(x1,grid->box[bln1],ijk1);

      for (j2 = 0; j2 < m2; j2++)
      {
        double x2[3];
        int ijk2 = pair[i].bface2->fpts->point[bln2][j2];

        get_x_coord(x2,grid->box[bln2],ijk2);

        if (dequal(norm(x1,x2),0))
        {
          if (j1 != j2)
          {
            int tmp;
            tmp = pair[i].bface2->fpts->point[bln2][j1];
            pair[i].bface2->fpts->point[bln2][j1] = ijk2;
            pair[i].bface2->fpts->point[bln2][j2] = tmp;

          }
          else
            break;
        }
      }

    }//for (j1 = 0; j1 < m1; j1++)

  }//for (i = 0; i < np; i++)

  if (np > 0)
    free(pair);
}


/*set the correct ofi index if it is needed*/
static void set_ofi_flag(tGrid *grid)
{
  /*those bafaces which are paired*/
  struct PAIR_S *pair = 0, *pair2 = 0;
  FLAG_T flg;
  int b, np = 0, np2 = 0;
  int i;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      /* Check if this bface has already been counted */
      flg = check_bface(pair,np,box->bface[bf]);

      if (flg == CONTINUE_F) continue;

      flg = NONE_F;

      /*If they both touch and collocated*/
      if (box->bface[bf]->touch == 1 && box->bface[bf]->same_fpts == 1)
      {
        
        tBox *box2 = grid->box[box->bface[bf]->ob];
        int bf2;
        
        for (bf2 = 0; bf2 < box2->nbfaces; bf2++)
        {

          flg = check_bface(pair,np,box2->bface[bf2]);

          if (flg == CONTINUE_F) continue;

          if (box2->bface[bf2]->touch == 1 && box2->bface[bf2]->same_fpts == 1)
          {

            if ( box2->bface[bf2]->ofi == box->bface[bf]->f   &&
                 box2->bface[bf2]->ob  == box->bface[bf]->b   &&
                 box2->bface[bf2]->f   == box->bface[bf]->ofi &&
                 box2->bface[bf2]->b   == box->bface[bf]->ob     )
            {

              box2->bface[bf2]->ofi = bf;
              box->bface[bf]->ofi   = bf2;

              add_to_pair(&pair,box->bface[bf],box2->bface[bf2],&np);

              flg = STOP_F;
              break;

            }
          }

        }//for (bf2 = 0; bf2 < grid->box[b2]->nbfaces; bf2++)

        if (flg != STOP_F)
        {
          /* Since couldn't be found, save this bface for further analysis */
          add_to_pair(&pair2,box->bface[bf],0,&np2);
          
          /* book keeping */
          add_to_pair(&pair,box->bface[bf],0,&np);
          
        }

      }
      /* If they both only touch */
      else if (box->bface[bf]->touch == 1 && box->bface[bf]->same_fpts != 1)
      {
        
        tBox *box2 = grid->box[box->bface[bf]->ob];
        int bf2;

        for (bf2 = 0; bf2 < box2->nbfaces; bf2++)
        {
          flg = check_bface(pair,np,box2->bface[bf2]);

          if (flg == CONTINUE_F) continue;

          if (box2->bface[bf2]->touch == 1 && box2->bface[bf2]->same_fpts != 1)
          {

            if ( box2->bface[bf2]->ofi == box->bface[bf]->f   &&
                 box2->bface[bf2]->ob  == box->bface[bf]->b   &&
                 box2->bface[bf2]->f   == box->bface[bf]->ofi &&
                 box2->bface[bf2]->b   == box->bface[bf]->ob    )

            {
              box2->bface[bf2]->ofi = bf;
              box->bface[bf]->ofi   = bf2;

              add_to_pair(&pair,box->bface[bf],box2->bface[bf2],&np);

              flg = STOP_F;
              break;

            }
          }

        }//for (bf2 = 0; bf2 < grid->box[b2]->nbfaces; bf2++)

        /* No interpolation bface could be found, so treat it as it goes inside
            the other box ob, like untouched one*/
        if (flg != STOP_F)
        {
          box->bface[bf]->touch = 0;
          box->bface[bf]->ofi = -1;
          add_to_pair(&pair,box->bface[bf],0,&np);
          
        }
      }// else if (box->bface[bf]->touch == 1 && box->bface[bf]->same_fpts != 1)

      /* If this bface is untouch */
      else if (box->bface[bf]->touch == 0 && box->bface[bf]->ofi == -1)
      {
        add_to_pair(&pair,box->bface[bf],box->bface[bf],&np);
      }// else if (box->bface[bf]->touch == 0 && box->bface[bf]->ofi == -1)

      else
      {
        errorexit("This case is not considered!\n");
      }

    }//for (bf = 0; bf < box->nbfaces; bf++)
  }

  /* Find the bfaces which couldn't be found */
  i = 0;
  while (i < np2)
  {
    
    set_ofi_flag_exclusively(&pair2[i],grid);
    
    i++;
  }
  
  if (np > 0)   free(pair);
  if (np2 > 0)  free(pair2);
}

/* set ofi flag exclusively by making extra bfaces, which is basically harmless*/
static void  set_ofi_flag_exclusively(struct PAIR_S *pair,tGrid *grid)
{
  tBface *const bface = pair->bface1;
  
  /*If they both touch and collocated*/
  if (bface->touch == 1 && bface->same_fpts == 1)
  {
    int b_new = bface->ob;
    int f_new = bface->ofi;
    tBox *box_new = grid->box[b_new];
    int fi_new = add_empty_bface(box_new,f_new);
    tBface *bface_new = box_new->bface[fi_new];
    int npoints;
    int ijk_new;
    int i;
    
    bface_new->b = b_new;
    bface_new->f = f_new;
    bface_new->ob = bface->b;
    bface_new->ofi = bface->fi;
    bface_new->touch = 1;
    bface_new->same_fpts = 1;
    
    bface->ofi = fi_new;
    
    if(bface_new->fpts == NULL)
      bface_new->fpts = AllocatePointList(grid);
    
    npoints = bface->fpts->npoints[bface->b];
    i = 0;
    while(i < npoints)
    {
      double x[3];
      
      get_x_coord(x,grid->box[bface->b],bface->fpts->point[bface->b][i]);
      ijk_new = find_ijk(box_new,x);
      assert(ijk_new >= 0);
      
      AddToPointList(bface_new->fpts, b_new, ijk_new);
      
      i++;
    }

    
  }

}

/*Given collocated point x and box, it returns the ijk of point*/
static int find_ijk(tBox *box, double x[3])
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double x2[3];
  int i,j,k;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < n2; j++)
      for (k = 0; k < n3; k++)
      {
        get_x_coord(x2,box,Index(i,j,k));

        if (dequal(norm(x,x2),0))
        {
          return Index(i,j,k);
        }
      }
  }

  return -1;
 
}

/*check if this bface is not considered previously*/
static FLAG_T check_bface(struct PAIR_S *pair,int np, tBface *bface)
{
  int i;

  for (i = 0; i < np; i++)
  {
    if (pair[i].bface1 == bface || pair[i].bface2 == bface)
      return CONTINUE_F;
  }

  return NONE_F;
}

/*finding the remaining bfeces*/
static void find_remaining_bfaces(struct FACE_POINT_S ***const FacePoint,int b, int f,tGrid *grid)
{
  int i;

  for (i = 0; i < FacePoint[b][f]->sh; i++)
  {
    struct SHARED_S *shared;
    struct SIMILAR_S P_touch, P_untouch;
    struct SIMILAR_S P_same_fpts, P_unsame_fpts;
    tBface *bface;

    P_touch.np = 0;
    P_touch.P = NULL;
    P_untouch.np = 0;
    P_untouch.P = NULL;

    P_same_fpts.np = 0;
    P_same_fpts.P = NULL;
    P_unsame_fpts.np = 0;
    P_unsame_fpts.P = NULL;

    shared = &FacePoint[b][f]->shared[i];

    group_point(TOUCH_F,&P_touch,&P_untouch,shared);

    if (P_touch.np > 0)
    {
      group_point(SAME_FTPS_F,&P_same_fpts,&P_unsame_fpts,&P_touch);
      free(P_touch.P);

      if (P_same_fpts.np > 0)
      {
        bface = make_bface(&P_same_fpts, grid);
        bface->touch = 1;
        bface->same_fpts = 1;
        bface->ob = P_same_fpts.P[0]->adjacent.box;
        bface->ofi = P_same_fpts.P[0]->adjacent.face;
        /* NOTE: here we set ofi as the other face and later we will fix it to other face index */
        free(P_same_fpts.P);
      }

      if (P_unsame_fpts.np > 0)
      {

        bface = make_bface(&P_unsame_fpts, grid);
        bface->touch = 1;
        bface->same_fpts = 0;
        bface->ob = P_unsame_fpts.P[0]->adjacent.box;
        bface->ofi = P_unsame_fpts.P[0]->adjacent.face;
        /* NOTE: here we set ofi as the other face and later we will fix it to other face index */
        free(P_unsame_fpts.P);
      }
    }

    if (P_untouch.np > 0)
    {
      bface = make_bface(&P_untouch,grid);
      bface->ob = P_untouch.P[0]->adjacent.box;
      bface->ofi = -1;
      free(P_untouch.P);
    }

  }

}

/*making the bfaces and return pointer to int*/
static tBface *make_bface(struct SIMILAR_S *P, tGrid *grid)
{
  int i, fi;
  int b = P->P[0]->geometry.b;
  int f = P->P[0]->geometry.f;
  tBox *box = grid->box[b];
  fi = add_empty_bface(box, f);
  tBface *bface = box->bface[fi];
  bface->f = f;
  bface->b = b;
  assert(f < TOT_NUM_FACE);

  if(bface->fpts==NULL)
    bface->fpts = AllocatePointList(box->grid);

  for (i = 0; i < P->np; i++)
  {
    AddToPointList(bface->fpts, b, P->P[i]->geometry.ijk);
  }

  return bface;
}

/*Given a flag, its group to gether all point with the flag feature
in to sturcture P1 and the rest in stucture P2*/
static void group_point(FLAG_T kind, \
  struct SIMILAR_S *P1, struct SIMILAR_S *P2, void *input)
{
  int i;

  if (kind == TOUCH_F)
  {
    struct SHARED_S *P = input;

    for (i = 0; i < P->np; i++)
    {
      struct POINT_S *pnt = P->pnt(P->FacePoint,P->kind[i],P->indx[i]);

      if (pnt->adjacent.touch == 1)
      {
        add_point(P1,pnt);
      }
      else
      {
        add_point(P2,pnt);
      }
    }
  }

  else if (kind == SAME_FTPS_F)
  {
    struct SIMILAR_S *P = input;
    for (i = 0; i < P->np; i++)
    {
      if (P->P[i]->adjacent.interpolation == 1)
      {
        add_point(P2,P->P[i]);
      }
      else
      {
        add_point(P1,P->P[i]);
      }
    }
  }

}

/*Returning the address of a point in FacePoint stuct*/
struct POINT_S *point_in_FacePoint(struct FACE_POINT_S *FacePoint,int kind, int indx)
{
  if (kind == INNER_F)
  {
    return &FacePoint->inner[indx];
  }
  else if (kind == EDGE_F)
  {
    return &FacePoint->edge[indx];
  }
}

/*Adding a needle to hay*/
static void add_point(struct SIMILAR_S *hay,struct POINT_S *needle)
{
  hay->np++;
  hay->P = realloc(hay->P,hay->np*sizeof(*hay->P));
  hay->P[hay->np-1] = needle;
}

/*find outerbound bfaces, if it found a outerbound in returns continue flg*/
static int find_outerbound_bfaces(struct FACE_POINT_S ***const FacePoint,
int b, int f,tGrid *grid)
{
  if (FacePoint[b][f]->outerbound == 1)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int bfacei,i;
    int l = 0, s = 0;

    bfacei = add_empty_bface(box,f);

    box->bface[bfacei]->outerbound = 1;

    /*Filling point list for inner*/
    for (i = 0; i < FacePoint[b][f]->s; i++)
    {
      if (FacePoint[b][f]->inner[i].adjacent.outerbound == 1)
      {
        add_point_to_bface_inbox(box,bfacei,\
            FacePoint[b][f]->inner[i].geometry.ijk,f);
        s++;
      }
    }

    /*Filling point list for edge*/
    for (i = 0; i < FacePoint[b][f]->l; i++)
    {

      if (FacePoint[b][f]->edge[i].adjacent.outerbound == 1)
      {
        add_point_to_bface_inbox(box,bfacei,\
          FacePoint[b][f]->edge[i].geometry.ijk,f);
        l++;
      }
    }

    if (l == FacePoint[b][f]->l && s == FacePoint[b][f]->s)
      return CONTINUE_F;
    else
      return NONE_F;
  }

  else
    return NONE_F;
}

/*Finding the adjacent structure for edge points*/
static void find_adjacent_edge(struct FACE_POINT_S ***const FacePoint,tGrid *grid)
{
  int *inbox = 0,nb = 0;
  int *blist = 0;
  double X[3] = {0};
  int b,b_max;
  
  b_max = tot_num(FacePoint);

  for (b = 0; b < b_max; b++)
  {
    int n1 = grid->box[b]->n1;
    int n2 = grid->box[b]->n2;
    int n3 = grid->box[b]->n3;
    int f, f_max;

    f_max = tot_num(FacePoint[b]);

    for (f = 0; f < f_max; f++)
    {
      const int l = FacePoint[b][f]->l;
      const int sh = FacePoint[b][f]->sh;
      int i;
      
      blist = realloc(blist,FacePoint[b][f]->sh*sizeof(*blist));
      
      /* if this face is internal face */
      if (FacePoint[b][f]->internal_face == 1)
      {
        continue;
      }

      
      for (i = 0; i < sh; i++)
      {
        blist[i] = FacePoint[b][f]->shared[i].box;
      }

      for (i = 0; i < l; i++)
      {
        struct POINT_S *edge = FacePoint[b][f]->edge;
       
        double *x = FacePoint[b][f]->edge[i].geometry.x;
        int k;

        /* if this point reach outerbound */
        double q[3] = { edge[i].geometry.x[0]+edge[i].geometry.N[0]*EPS,\
                          edge[i].geometry.x[1]+edge[i].geometry.N[1]*EPS,\
                          edge[i].geometry.x[2]+edge[i].geometry.N[2]*EPS};
        
        int out = b_XYZ_of_xyz(grid,&X[0],&X[1],&X[2],q[0],q[1],q[2]);
        
        if (out == b)
        {
          int *obl = 0; /* list */
          int ii,jj;
          
          /* make obl that contains all boxes */
          jj = 0;
          forallboxes(grid, ii) 
            if(ii != b) 
            {
              obl = realloc(obl,(jj+1)*sizeof(*obl));
              obl[jj++] = ii;
            }
          
          out = b_XYZ_of_xyz_inboxlist(grid,obl,ii-1,&X[0],&X[1],&X[2],q[0],q[1],q[2]);
          
          free(obl);
        }
        
        if (out < 0 && FacePoint[b][f]->outerbound == 1)
        {
          FacePoint[b][f]->edge[i].adjacent.outerbound = 1;
          continue;
        }

        inbox = 0;
        nb = 0;

        for (k = 0; k < sh; k++)
        {
          
          int b_ = b_XYZ_of_xyz_inboxlist(grid,&blist[k],1,&X[0],&X[1],&X[2],x[0],x[1],x[2]);
          
          if ( b_ >= 0 )
          {
            inbox = realloc(inbox,(nb+1)*sizeof(*inbox));
            inbox[nb] = b_;
            nb++;
          }
        }

        if (nb == 0) //Find the adjacent box;
        {
          blist = realloc(blist,(sh+1)*sizeof(*blist));
          assert(blist != 0);
          
          blist[sh] = b;
          
          inbox = b_xyz_in_exblist(grid,x,blist,sh+1,&nb);
          
          if (nb == 0 && out < 0)
          {
            FacePoint[b][f]->outerbound = 1;
            FacePoint[b][f]->edge[i].adjacent.outerbound = 1;
            continue;
          }
          else if (nb == 0)//if still couldn't be found
          {
            fprintf(stderr,"Point = (%f,%f,%f)\n",x[0],x[1],x[2]);
            errorexit("There is no box for the above point!\n");
          }
          
        }

        populate_adjacent(FacePoint,EDGE_F,b,f,i,inbox,nb,grid);

        if (nb > 0)
          free(inbox);
      }

      

    }//for (f = 0; f < f_max; f++)
  }//for (b = 0; b < b_max; b++)
  
  free (blist);
}

/*Finding the adjacent structure for inner points*/
static void find_adjacent_inner(struct FACE_POINT_S ***const FacePoint,tGrid *grid)
{
  double X[3] = {0};
  int blist;
  int b,b_max;
  
  b_max = tot_num(FacePoint);

  for (b = 0; b < b_max; b++)
  {
    int n1 = grid->box[b]->n1;
    int n2 = grid->box[b]->n2;
    int n3 = grid->box[b]->n3;
    int f, f_max;

    f_max = tot_num(FacePoint[b]);
    blist = b;

    for (f = 0; f < f_max; f++)
    {
      int *adj_box,n_adj_box;//list and number of adjacent boxes
      const int s = FacePoint[b][f]->s;
      int i;

      /* if this face is internal face */
      if (FacePoint[b][f]->internal_face == 1)
          continue;

      for (i = 0; i < s; i++)
      {
        struct POINT_S *inner = FacePoint[b][f]->inner;
        double q[3] = { inner[i].geometry.x[0]+inner[i].geometry.N[0]*EPS,
                        inner[i].geometry.x[1]+inner[i].geometry.N[1]*EPS,
                        inner[i].geometry.x[2]+inner[i].geometry.N[2]*EPS };
                        
        /* This set the priority to outerbound or innerbound */
        int out = b_XYZ_of_xyz(grid,X,X+1,X+2,q[0],q[1],q[2]);
        
        if (out == b)
        {
          int *obl = 0; /* list */
          int ii,jj;
          
          /* make obl that contains all boxes */
          jj = 0;
          forallboxes(grid, ii) 
            if(ii != b) 
            {
              obl = realloc(obl,(jj+1)*sizeof(*obl));
              obl[jj++] = ii;
            }
          
          out = b_XYZ_of_xyz_inboxlist(grid,obl,ii-1,&X[0],&X[1],&X[2],q[0],q[1],q[2]);
          
          free(obl);
        }
        
        if (out < 0)
        {
          FacePoint[b][f]->inner[i].adjacent.outerbound = 1;
          FacePoint[b][f]->outerbound = 1;// Flag this face for future uses
          continue;
        }

        adj_box = 0;
        n_adj_box = 0;

        adj_box = b_xyz_in_exblist(grid,inner[i].geometry.x,&blist,1,&n_adj_box);
        
        /*if the point couldn't be found, let's try the normal approach*/
        if (n_adj_box == 0)
        {
        
          adj_box = malloc(sizeof(*adj_box));
          assert(adj_box != 0);
          n_adj_box++;
          adj_box[0] = b_XYZ_of_xyz(grid,&X[0],&X[1],&X[2],q[0],q[1],q[2]);

          if (adj_box[0] < 0)//the point is considered as outer or inner bound
            n_adj_box = 0;
        }

        if (n_adj_box > 0)
          populate_adjacent(FacePoint,INNER_F,b,f,i,adj_box,n_adj_box,grid);
        else
        {
          FacePoint[b][f]->inner[i].adjacent.outerbound = 1;
          FacePoint[b][f]->outerbound = 1;// Flag this face for future uses
        }

        free(adj_box);
      }//for (i = 0; i < s; i++)
    }//for (f = 0; f < f_max; f++)

  }//for (b = 0; b < b_max; b++)
}

/*given x finding the pertinent box(es) excluding ex_blist
  and put the answers in inbox and nb; inbox[0..nb-1]*/
static int *b_xyz_in_exblist(tGrid *grid,double *x,int *ex_blist,int ex_bn,int *nb)
{
  int *inbox, b, *b_list, j;
  double X[3] = {0};

  forallboxes(grid,b);

  const int k = b-ex_bn;

  if (k == 0)
  {
    errorexit("There is no excluding b list!");
  }

  b_list = malloc(k*sizeof(*b_list));
  assert(b_list != NULL);

  j = 0;
  forallboxes(grid,b)
  {
    int i, flg;

    flg = 0;
    for (i = 0; i < ex_bn; i++)
      if (ex_blist[i] == b)
        flg = 1;

    if (flg == 0)
    {
      b_list[j] = b;
      assert(j < k);
      j++;
    }

  }

  *nb = 0;
  inbox = 0;
  j = 0;
  while (j < k)
  {

    b = b_XYZ_of_xyz_inboxlist(grid,&b_list[j],1,&X[0],&X[1],&X[2],x[0],x[1],x[2]);

    if (b >= 0)
    {
      inbox = realloc(inbox,(*nb+1)*sizeof(*inbox));
      assert(inbox != 0);

      inbox[*nb] = b;
      ++*nb;
    }

    j++;
  }

  free(b_list);
  
  if (*nb == 0)
    return 0;
    
  return inbox;

}

/*adding info to share structure*/
static void add_to_share(struct FACE_POINT_S *const FacePoint,int f, int b, FLAG_T kind, int indx)
{
  const int sh = FacePoint->sh;

  /*Make sure the identical box and face won't be double counted*/
  if (sh != 0)
  {
    int i = 0;

    while (i < sh)
    {
      if ( FacePoint->shared[i].box == b &&
           FacePoint->shared[i].face == f)
      {
        FacePoint->shared[i].indx = realloc(FacePoint->shared[i].indx,\
                (FacePoint->shared[i].np+1)*sizeof(*FacePoint->shared[i].indx));
        FacePoint->shared[i].indx[FacePoint->shared[i].np] = indx;

        FacePoint->shared[i].kind = realloc(FacePoint->shared[i].kind,\
                (FacePoint->shared[i].np+1)*sizeof(*FacePoint->shared[i].kind));
        FacePoint->shared[i].kind[FacePoint->shared[i].np] = kind;
        FacePoint->shared[i].np++;
        return;
      }

      i++;
    }
  }

  FacePoint->shared = realloc(FacePoint->shared,\
                  (FacePoint->sh+1)*sizeof(*FacePoint->shared));

  FacePoint->shared[sh].FacePoint = FacePoint;
  FacePoint->shared[sh].face = f;
  FacePoint->shared[sh].box = b;
  FacePoint->shared[sh].indx = malloc(sizeof(*FacePoint->shared[sh].indx));
  FacePoint->shared[sh].indx[0] = indx;
  FacePoint->shared[sh].kind = malloc(sizeof(*FacePoint->shared[sh].kind));
  FacePoint->shared[sh].kind[0] = kind;
  FacePoint->shared[sh].pnt = point_in_FacePoint;
  FacePoint->shared[sh].np = 1;

  FacePoint->sh++;

}

/*Absolute value*/
static double ABS(double x)
{
  return x > 0 ? x : -x;
}

/*Dot product*/
static double dot_product(double *x,double *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

/*If the found point is at a face and it
 also check if it is at more that one face and the faces number saved in f,
 furthermore if the f[i] < 0 it means that there is no more face*/
static int IsAtFace(int box_n, int *ijk,int *f,tGrid *grid,struct FACE_POINT_S ***const FacePoint)
{
  tBox *box = grid->box[box_n];
  int n[3] = {box->n1,box->n2,box->n3};
  struct FACE_POINT_S **fp = FacePoint[box_n];
  unsigned int flg = 0;
  int i;

  for (i = 0; i < TOT_NUM_FACE; i++)
    f[i] = -1;

  i = 0;
  if (ijk[0] == 0 && fp[X_FACE0]->internal_face == 0)
  {
    f[i++] = X_FACE0;
    flg = 1;
  }
  if (ijk[0] == n[0]-1 && fp[X_FACE1]->internal_face == 0)
  {
    f[i++] = X_FACE1;
    flg = 1;
  }
  if (ijk[1] == 0  && fp[Y_FACE0]->internal_face == 0)
  {
    f[i++] = Y_FACE0;
    flg = 1;
  }
  if (ijk[1] == n[1]-1 && fp[Y_FACE1]->internal_face == 0)
  {
    f[i++] = Y_FACE1;
    flg = 1;
  }
  if (ijk[2] == 0 &&  fp[Z_FACE0]->internal_face == 0)
  {
    f[i++] = Z_FACE0;
    flg = 1;
  }
  if (ijk[2] == n[2]-1 && fp[Z_FACE1]->internal_face == 0)
  {
    f[i++] = Z_FACE1;
    flg = 1;
  }

  return flg;

}

/*If the point x collocated with the same point in adjacent point*/
static int IsCollocated(double *x1, int box2, int *ijk_adj,double *X,tGrid *grid)
{
  tBox *box = grid->box[box2];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double x2[3];
  int i,j,k;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < n2; j++)
      for (k = 0; k < n3; k++)
      {
        get_x_coord(x2,box,Index(i,j,k));

        if (dequal(norm(x1,x2),0))
        {
          ijk_adj[0] = i;
          ijk_adj[1] = j;
          ijk_adj[2] = k;
          get_X_coord(X,box,Index(i,j,k));

          return 1;
        }
      }
  }

  return 0;
}

/*Calulcating the Euclidean norm*/
static double norm(double *x1,double *x2)
{
  return sqrt(SQ(x1[0]-x2[0])+SQ(x1[1]-x2[1])+SQ(x1[2]-x2[2]));
}

/*Allocating memory for FacePoint struct*/
static void *allc_FacePoint(tGrid *grid)
{
  struct FACE_POINT_S ***FacePoint;
  int b,f;

  forallboxes(grid,b);

  FacePoint = malloc((b+1)*sizeof(*FacePoint));
  assert(FacePoint != NULL);
  FacePoint[b] = NULL;

  forallboxes(grid,b)
  {
    FacePoint[b] = malloc((TOT_NUM_FACE+1)*sizeof(*FacePoint[b]));
    assert(FacePoint[b] != NULL);
    FacePoint[b][TOT_NUM_FACE] = NULL;
  }

  forallboxes(grid,b)
    for (f = 0; f < TOT_NUM_FACE; f++)
    {
      FacePoint[b][f] = calloc(1,sizeof(*FacePoint[b][f]));
      assert(FacePoint[b][f] != NULL);
    }

  return FacePoint;
}

/*Initializing inner and edge points*/
static void init_FacePoint(struct FACE_POINT_S ***FacePoint,tGrid *grid)
{
  int b, f;

  assert(TOT_NUM_FACE == 6);/*Since this algorithm only works for this*/

  /*Allocate memory for points*/
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n[3] = {box->n1,box->n2,box->n3};

    for (f = 0; f < TOT_NUM_FACE; f++)
    {

      int m[2];
      int i_l,j_l,k_l;/*lower range of i, j and k*/
      int i_u,j_u,k_u;/*upper range of i, j and k*/

      setup_range(f,n,&i_l,&j_l,&k_l,&i_u,&j_u,&k_u,m);

      /*total number of inner points*/
      int s = (m[0]-2)*(m[1]-2);
      /*total number of edge points*/
      int l = 2*(m[0]+m[1]-2);

      FacePoint[b][f]->s = s;
      FacePoint[b][f]->l = l;
      FacePoint[b][f]->sh = 0;

      FacePoint[b][f]->inner = \
        calloc(s,sizeof(*FacePoint[b][f]->inner));
      assert(FacePoint[b][f]->inner != NULL);

      FacePoint[b][f]->edge = \
        calloc(l,sizeof(*FacePoint[b][f]->edge));
      assert(FacePoint[b][f]->edge != NULL);

      FacePoint[b][f]->shared = NULL;
    }
  }

  /*Filling the geometry structure*/
  fill_geometry(FacePoint,grid);
}

/*Filling the geometry structure*/
static void fill_geometry(struct FACE_POINT_S ***FacePoint,tGrid *grid)
{
  int b, f;
  FLAG_T flg = NONE_F;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int extface[TOT_NUM_FACE]={0};
    int n[3] = {box->n1,box->n2,box->n3};
    int n1 = n[0] ,n2 = n[1];//They are needed for macro Index

    find_external_faces_of_box(box,extface,1);

    for (f = 0; f < TOT_NUM_FACE; f++)
    {
      int l = FacePoint[b][f]->l;
      int s = FacePoint[b][f]->s;
      struct POINT_S *inner = FacePoint[b][f]->inner;
      struct POINT_S *edge  = FacePoint[b][f]->edge;
      struct POINT_S *P;

      int i,j,k,e,in;/*dummy indices*/
      int i_l,j_l,k_l;/*lower range of i, j and k */
      int i_u,j_u,k_u;/*upper range of i, j and k */

      if ( extface[f] == 0 )
      {
        FacePoint[b][f]->internal_face = 1;
        continue;
      }

      /*Filling the edge and inner points*/
      setup_range(f,n,&i_l,&j_l,&k_l,&i_u,&j_u,&k_u,NULL);
      e = 0;
      in = 0;
      for (i = i_l; i <= i_u; i++)
      {
        for (j = j_l; j <= j_u; j++)
        {
         for (k = k_l; k <= k_u; k++)
         {

           if (IsEdge(f,n,i,j,k) == 1)/*If it is on edge*/
           {
             assert(e < l);
             P = &edge[e];
             e++;
           }
           else
           {
             assert(in < s);
             P = &inner[in];
             in++;
           }

           get_x_coord(P->geometry.x,box,Index(i,j,k));
           get_X_coord(P->geometry.X,box,Index(i,j,k));

           P->geometry.ijk = Index(i,j,k);
           P->geometry.b = b;
           P->geometry.f = f;

           get_normal(P->geometry.N,box,f,Index(i,j,k));


         }
        }
      }
    }
  }//End of forallboxes
}

/*Gettig the normalized normal in outward direction*/
static void get_normal(double *N,tBox *box,int face,long int ijk)
{
  double n[4];

  boxface_outwarddir_at_ijk(box,face,ijk,n);

  N[0] = n[1];
  N[1] = n[2];
  N[2] = n[3];

}

static void normalizing_N(double *N)
{
  double n = sqrt(SQ(N[0])+SQ(N[1])+SQ(N[2]));

  assert(!dequal(n,0));

  N[0] /= n;
  N[1] /= n;
  N[2] /= n;

}
/*Finding the total number of a row in a pointer to a structure*/
static int tot_num(void *p)
{
  int n;
  void **pp = p;

  for (n = 0; pp[n] != NULL; n++);

  return n;
}

/*Get the value of coordinate x,y,z*/
static void get_x_coord(double *x,tBox *box,long int ind)
{
  int ix, iy, iz;

  ix = Ind("x");
  iy = Ind("y");
  iz = Ind("z");

  x[0] = box->v[ix][ind];
  x[1] = box->v[iy][ind];
  x[2] = box->v[iz][ind];
}

/*Get the value of coordinate X,Y,Z*/
static void get_X_coord(double *X,tBox *box,long int ind)
{

  int iX, iY, iZ;

  iX = Ind("X");
  iY = Ind("Y");
  iZ = Ind("Z");

  X[0] = box->v[iX][ind];
  X[1] = box->v[iY][ind];
  X[2] = box->v[iZ][ind];
}

/*Determining the lower and upper ranges of i,j and k and also the number
of points in each direction of the face m[0] and m[1]*/
static void setup_range(int face,int *n,\
                                    int *i_l,int *j_l,int *k_l,\
                                    int *i_u,int *j_u,int *k_u, int *m)
{

 switch (face)
 {
    case X_FACE0:
      *i_l = 0;
      *i_u = 0;
      *j_l = 0;
      *j_u = n[1]-1;
      *k_l = 0;
      *k_u = n[2]-1;
      if (m != NULL)
      {
        m[0] = n[1];
        m[1] = n[2];
      }

      break;
    case X_FACE1:
      *i_l = n[0]-1;
      *i_u = n[0]-1;
      *j_l = 0;
      *j_u = n[1]-1;
      *k_l = 0;
      *k_u = n[2]-1;
      if (m != NULL)
      {
        m[0] = n[1];
        m[1] = n[2];
      }
      break;
    case Y_FACE0:
      *i_l = 0;
      *i_u = n[0]-1;
      *j_l = 0;
      *j_u = 0;
      *k_l = 0;
      *k_u = n[2]-1;
      if (m != NULL)
      {
        m[0] = n[0];
        m[1] = n[2];
      }
      break;
    case Y_FACE1:
      *i_l = 0;
      *i_u = n[0]-1;
      *j_l = n[1]-1;
      *j_u = n[1]-1;
      *k_l = 0;
      *k_u = n[2]-1;
      if (m != NULL)
      {
        m[0] = n[0];
        m[1] = n[2];
      }
      break;
    case Z_FACE0:
      *i_l = 0;
      *i_u = n[0]-1;
      *j_l = 0;
      *j_u = n[1]-1;
      *k_l = 0;
      *k_u = 0;
      if (m != NULL)
      {
        m[0] = n[0];
        m[1] = n[1];
      }
      break;
    case Z_FACE1:
      *i_l = 0;
      *i_u = n[0]-1;
      *j_l = 0;
      *j_u = n[1]-1;
      *k_l = n[2]-1;
      *k_u = n[2]-1;
      if (m != NULL)
      {
        m[0] = n[0];
        m[1] = n[1];
      }
      break;
  }

}

/*Check if it is on edge*/
static int IsEdge(int f, int *n, int i, int j, int k)
{
  if (f == X_FACE0 || f == X_FACE1)
  {
    if ( j == 0 || j == n[1]-1 || k == 0 || k == n[2]-1 )
      return 1;
    else
      return 0;
  }

  else if (f == Y_FACE0 || f == Y_FACE1)
  {
    if ( i == 0 || i == n[0]-1 || k == 0 || k == n[2]-1 )
      return 1;
    else
      return 0;
  }

  else if (f == Z_FACE0 || f == Z_FACE1)
  {
    if ( j == 0 || j == n[1]-1 || i == 0 || i == n[0]-1 )
      return 1;
    else
      return 0;
  }
}

/*Populating the adjacent structure for each kind*/
static void populate_adjacent(struct FACE_POINT_S ***const FacePoint,FLAG_T kind,int b, int f, long int i, int *adj_box, int n_adj_box,tGrid *grid)
{
  /*Adjacent geometry structure*/
  struct ADJ_GEO_S
  {
    int box;// adjacent box
    int face;// ajacent face
    int collocated; // 1 means yes 0 no
    int AtFace; // 1 means yes 0 no
    int interpolation;// 1 means yes 0 no
    int touch;// 1 means yes 0 no
    int ijk[3];
    int conditional; //if it is conditional 1, otherwise 0
    double nrm;// dot products of normal N1.N2
    
  } *record = 0;
  struct POINT_S *P;
  double X[3] = {0};
  int ijk_adj[3];//the index of the point in adjacent box
  int face_list[TOT_NUM_FACE];
  FLAG_T flg;
  int k;

  if (kind == INNER_F)
  {
    P = &FacePoint[b][f]->inner[i];
  }
  else if (kind == EDGE_F)
  {
    P = &FacePoint[b][f]->edge[i];
  }

  record = calloc(n_adj_box,sizeof(*record));
  k = 0;
  while (k < n_adj_box)
  {
    record[k].box = adj_box[k];
    
    /*If the point x collocated with the
       same point in adjacent point*/
    if (IsCollocated(P->geometry.x,record[k].box,record[k].ijk,P->adjacent.X,grid) == 1)
    {
      record[k].collocated = 1;

      /*If the found point is at face*/
      if (IsAtFace(record[k].box,record[k].ijk,face_list,grid,FacePoint) == 1)
      {
        double N[3];
        int i = 0;

        record[k].AtFace = 1;
        record[k].touch = 1;
        flg = NONE_F;

        i = 0;
        while (face_list[i] >= 0)
        {
          double nrm;

          get_normal(N,grid->box[record[k].box],\
            face_list[i],ijk_ind(grid->box[record[k].box],record[k].ijk));
          nrm = ABS(dot_product(P->geometry.N,N));

          if (dequal(nrm,1))
          {
            record[k].face = face_list[i];
            record[k].interpolation = 0;
            record[k].nrm = nrm;
            flg = FOUND_F;
            break;
          }
          i++;
        }
        /* if the point is a conditional point */
        if (flg == NONE_F)
        {
          int indx0;
          double nrm0 = 0;
          
          /*Find the fitest norm*/
          i = 0;
          while (face_list[i] >= 0)
          {
            double nrm;

            get_normal(N,grid->box[record[k].box],\
              face_list[i],ijk_ind(grid->box[record[k].box],record[k].ijk));
            nrm = ABS(dot_product(P->geometry.N,N));
            
            if (dgreatereq(nrm,nrm0))
            {
              nrm0 = nrm;
              indx0 = i;
            }
            
            i++;
          }
          
          record[k].face = face_list[indx0];
          record[k].interpolation = 0;
          record[k].nrm = nrm0;
          record[k].conditional = 1;
        }

      }
      /* It's not at face so it is inside the mesh and it
       * needs iterpolation and it is considered as untouch one */
      else
      {
        record[k].AtFace = -1;
        record[k].touch = 0;
      }

    }
    else/*Since the found point is not collocated it needs interpolation*/
    {
      /* check if point is on a face or not, i.e. touching or not */
      /* if it is on corner point pick out the correct face; since in this case we don't have
      collocated point, the appropriate normal will be chosen based on the closest point*/
      /* add this info into the share structure */

     
      int face[TOT_NUM_FACE];
      int nf = 0,i;

      XYZ_of_xyz(grid->box[record[k].box],X,X+1,X+2,
                P->geometry.x[0],P->geometry.x[1],P->geometry.x[2]);

      nf =
          XYZ_on_face(grid->box[record[k].box],face,X[0],X[1],X[2]);
      
      /* make sure you avoid the internal face*/    
      
      for (i = 0; i < TOT_NUM_FACE; i++)
      {
        if (face[i] == 1 && FacePoint[record[k].box][i]->internal_face == 1)
        {
          face[i] = 0;
          nf--;
        }
      }
      
      record[k].interpolation = 1;
      if (nf == 0)
      {
        record[k].touch = 0;
        record[k].face = -1;
      }
      else
      {
        struct FACE_NORMAL_S *fc_nr = calloc(nf,sizeof(*fc_nr));
        assert(fc_nr != 0);

        int j = 0;
        for (i = 0; i < TOT_NUM_FACE; i++)
        {
          if ( face[i] == 1 )
          {
            fc_nr[j].face = i;
            /* getting the approximate normal base on closest point */
            get_apprx_normal(&fc_nr[j],grid->box[record[k].box],i,P->geometry.x);

            j++;

          }
        }

        assert(j > 0);

        double nrm = -1;
        int keep_face = -1;
        for (i = 0; i < j; i++)
        {
          double nrm2 = ABS(dot_product(P->geometry.N,fc_nr[i].N));
          if (dgreatereq(nrm2,nrm))
          {
            keep_face = fc_nr[i].face;
            nrm = nrm2;
            record[k].nrm = nrm;
          }
        }

        assert(keep_face >= 0 );
        record[k].face = keep_face;
        record[k].touch = 1;

        free(fc_nr);
      }

    }

    k++;
  }

  /* Selecting the fitest adjacent */
  /* The fitest is:
  be collocated
  be AtFace
  match normal
  needs intepolation
  the matches normal
  finally the untouch one*/

  /* First look for collocation point */

  int n_identical = 0;/* For identical cases in which we have let say
                       * three bfaces with a common line with collocated points */
  int *desired_indx = 0;
  k = 0;
  flg = NONE_F;
  while(k < n_adj_box)
  {
    if (record[k].collocated == 1 &&
        record[k].AtFace     == 1 &&
        record[k].conditional == 0)
    {
      n_identical++;
      desired_indx = realloc(desired_indx,(n_identical)*sizeof(*desired_indx));
      desired_indx[n_identical-1] = k;
      flg = FOUND_F;
    }
    k++;
  }

  if (n_identical >= 1)
  {
    k = desired_indx[0];
    P->adjacent.box = record[k].box;
    P->adjacent.ijk[0] = record[k].ijk[0];
    P->adjacent.ijk[1] = record[k].ijk[1];
    P->adjacent.ijk[2] = record[k].ijk[2];
    P->adjacent.touch = record[k].touch;
    P->adjacent.face = record[k].face;
    P->adjacent.interpolation = record[k].interpolation;
    assert(P->adjacent.touch == 1);
    assert(P->adjacent.interpolation == 0);
    assert(dequal(record[k].nrm,1));
    add_to_share(FacePoint[b][f],P->adjacent.face,P->adjacent.box,kind,i);

    if (n_identical > 1) // we have identical cases
    {
      int l,end;
      struct POINT_S *P1 = duplicate_points(FacePoint,kind,b,f,i,n_identical,&end);

      for (l = 1; l < n_identical; l++)
      {
        k = desired_indx[l];
        P1[l-1].adjacent.box = record[k].box;
        P1[l-1].adjacent.ijk[0] = record[k].ijk[0];
        P1[l-1].adjacent.ijk[1] = record[k].ijk[1];
        P1[l-1].adjacent.ijk[2] = record[k].ijk[2];
        P1[l-1].adjacent.touch = record[k].touch;
        P1[l-1].adjacent.face = record[k].face;
        P1[l-1].adjacent.interpolation = record[k].interpolation;
        assert(P1[l-1].adjacent.touch == 1);
        assert(P1[l-1].adjacent.interpolation == 0);
        assert(dequal(record[k].nrm,1));
        add_to_share(FacePoint[b][f],P1[l-1].adjacent.face,P1[l-1].adjacent.box,kind,end+l-1);

      }

    }

    free(desired_indx);
  }

  /* Then, if collocation point cannot be found, look for interpolation points */
  if (flg != FOUND_F)
  {
    int index = -1;
    k = 0;
    while(k < n_adj_box)
    {
      if (record[k].interpolation == 1 &&
          record[k].touch         == 1 &&
          record[k].face          >= 0   )
      {
        if (index != -1 && dgreatereq(record[k].nrm,record[index].nrm))
          index = k;
        else if (index == -1)
          index = k;
      }
      k++;
    }

    if (index != -1)
    {
      /* If there is an ambiguity between outerbound and interpolation,
       * one has to study the neighbor of the point to realize the status */
      if (dequal(record[index].nrm,0))
      {
        struct FACE_NORMAL_S *fc_nr = 0;
        double N[3] = {P->geometry.N[0],P->geometry.N[1],P->geometry.N[2]};
        int face[TOT_NUM_FACE];
        int i, nf = 0;

        XYZ_of_xyz(grid->box[record[index].box],X,X+1,X+2,
                P->geometry.x[0],P->geometry.x[1],P->geometry.x[2]);

        nf = XYZ_on_face(grid->box[record[index].box],face,X[0],X[1],X[2]);
        assert(nf != 0);

        fc_nr = calloc(nf,sizeof(*fc_nr));
        assert(fc_nr != 0);

        int j = 0;
        for (i = 0; i < TOT_NUM_FACE; i++)
        {
          if ( face[i] == 1 )
          {
            fc_nr[j].face = i;
            /* getting the approximate normal base on closest point */
            get_apprx_normal(&fc_nr[j],grid->box[record[index].box],i,P->geometry.x);

            N[0] += EPS2*fc_nr[j].N[0];
            N[1] += EPS2*fc_nr[j].N[1];
            N[2] += EPS2*fc_nr[j].N[2];

            j++;

          }
        }

        free(fc_nr);

        double Nsqr = sqrt(SQ(N[0])+SQ(N[1])+SQ(N[2]));
        assert(!dequal(Nsqr,0));
        N[0] /= Nsqr;
        N[1] /= Nsqr;
        N[2] /= Nsqr;

        double q[3] = { P->geometry.x[0]+N[0]*EPS,
                        P->geometry.x[1]+N[1]*EPS,
                        P->geometry.x[2]+N[2]*EPS };
        int out;
        out = b_XYZ_of_xyz(grid,&X[0],&X[1],&X[2],q[0],q[1],q[2]);

        if (out < 0)
        {
          P->adjacent.outerbound = 1;
        }
        else
        {
          P->adjacent.box = record[index].box;
          P->adjacent.touch = 0;
          P->adjacent.face = -1;
          P->adjacent.interpolation = 1;
          add_to_share(FacePoint[b][f],P->adjacent.face,P->adjacent.box,kind,i);
        }

      }
      else
      {
        P->adjacent.box = record[index].box;
        P->adjacent.touch = record[index].touch;
        P->adjacent.face = record[index].face;
        P->adjacent.interpolation = record[index].interpolation;
        add_to_share(FacePoint[b][f],P->adjacent.face,P->adjacent.box,kind,i);
      }

      flg = FOUND_F;
    }

  }

  /* If it still is not found check for untouch one*/
  if (flg != FOUND_F)
  {
    k = 0;
    while(k < n_adj_box)
    {
      if (record[k].interpolation == 1 &&
          record[k].touch         == 0 &&
          record[k].face          <  0   )
      {
        P->adjacent.box = record[k].box;
        P->adjacent.touch = record[k].touch;
        P->adjacent.face = record[k].face;
        P->adjacent.interpolation = record[k].interpolation;
        add_to_share(FacePoint[b][f],P->adjacent.face,P->adjacent.box,kind,i);
        flg = FOUND_F;
        break;
      }
      k++;
    }

  }
  
  /* If it still is not found check for collocated point 
  and compromise the normal */
  if (flg != FOUND_F)
  {
    int indx;
    double nrm = -1;
    
    k = 0;
    while(k < n_adj_box)
    {
      /* find the greatest nrm */
      if ( record[k].collocated == 1 && 
           record[k].AtFace     == 1 &&
           dgreatereq(record[k].nrm,nrm) )
      {
        nrm = record[k].nrm;
        indx = k;
      }
      
      k++;
    }
    
    k = indx;  
    P->adjacent.box = record[k].box;
    P->adjacent.ijk[0] = record[k].ijk[0];
    P->adjacent.ijk[1] = record[k].ijk[1];
    P->adjacent.ijk[2] = record[k].ijk[2];
    P->adjacent.touch = record[k].touch;
    P->adjacent.face = record[k].face;
    P->adjacent.interpolation = record[k].interpolation;
    add_to_share(FacePoint[b][f],P->adjacent.face,P->adjacent.box,kind,i);
    
    assert(P->adjacent.touch == 1);
    assert(P->adjacent.interpolation == 0);
    
    if (dgreater(record[k].nrm,0))
      flg = FOUND_F;
      
  }
  
  /* check if the pertinent adjacent has been found */
  if (flg != FOUND_F)
  {
    double *x = P->geometry.x;
    fprintf(stderr,"ERROR:\nThe point (%f,%f,%f) in box: %d and at face: %d cannot be matched!\n"
        "%s,%d\n",x[0],x[1],x[2],b,f,__FILE__, __LINE__);
    abort();
  }

  free(record);

}

/* duplicating the points with the same properties regarding bface and return a pointer to the beginning of new structure*/
static void *duplicate_points(struct FACE_POINT_S ***const FacePoint,FLAG_T kind,int b, int f,int indx,int n_dup, int *end)
{
  struct POINT_S **pp;
  int *n,i;

  if (kind == EDGE_F)
  {
    pp = &FacePoint[b][f]->edge;
    n = &FacePoint[b][f]->l;
  }
  else if (kind == INNER_F)
  {
    pp = &FacePoint[b][f]->inner;
    n = &FacePoint[b][f]->s;
  }

  *end = *n;
  *n += n_dup-1;
  (*pp) = realloc(*pp,*n*sizeof(**pp));

  /* copying structures */
  for (i = *n-n_dup+1; i < *n; i++)
  {
    (*pp)[i].geometry = (*pp)[indx].geometry;
  }

  return &(*pp)[*n-n_dup+1];
}

/* getting the approximate normal base on closest point */
static void get_apprx_normal(struct FACE_NORMAL_S *fc_nr,tBox *box, int face,double *x1)
{
  /* find the closest point to ijk in the box and on the face */
  fc_nr->ijk = closest_point(box,face,x1);

  get_normal(fc_nr->N,box,face,fc_nr->ijk);

}

/* find the closest point to x1 in the box and on the face */
static int closest_point(tBox *box,int face,double *x1)
{
  int n[3] = {box->n1,box->n2,box->n3};
  double nrm = DBL_MAX;
  int ijk2 = -1;
  int i_l,j_l,k_l;/*lower range of i, j and k*/
  int i_u,j_u,k_u;/*upper range of i, j and k*/
  int i,j,k;

  setup_range(face,n,&i_l,&j_l,&k_l,&i_u,&j_u,&k_u,0);

  for (i = i_l; i <= i_u; i++)
    for (j = j_l; j <= j_u; j++)
      for (k = k_l; k <= k_u; k++)
      {
        int ind[3] = {i,j,k};
        double x2[3];
        get_x_coord(x2,box,ijk_ind(box,ind));
        double nrm2 = norm(x1,x2);

        if (dlesseq(nrm2,nrm))
        {
          ijk2 = ijk_ind(box,ind);
          nrm = nrm2;
        }
      }

  assert(ijk2 >= 0);

  return ijk2;
}

/*Returning the index at box*/
static int ijk_ind(tBox *box, int *i)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;

  return Index(i[0],i[1],i[2]);
}

/* visualizing bfaces */
static void visualize_bfaces(tGrid *grid)
{
  struct PAIR_S *pair = 0;
  int b,np = 0;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      tBface *bface = box->bface[bf];

      if (bface->outerbound == 1)
      {
        print_bface(bface,NULL,"outerbound",pair,np);

        add_to_pair(&pair,bface,0,&np);
      }
      else if (bface->touch == 1 && bface->same_fpts == 1)
      {
        tBface *bface2 = grid->box[bface->ob]->bface[bface->ofi];
        print_bface(bface,bface2,"touch_and_same_fpts",pair,np);

        add_to_pair(&pair,bface,bface2,&np);
      }
      else if (bface->touch == 1 && bface->same_fpts == 0)
      {
        tBface *bface2 = grid->box[bface->ob]->bface[bface->ofi];
        print_bface(bface,bface2,"touch_and_dffrnt_fpts",pair,np);

        add_to_pair(&pair,bface,bface2,&np);
      }
      else if (bface->touch == 0)
      {
        print_bface(bface,NULL,"Untouch",pair,np);

        add_to_pair(&pair,bface,0,&np);
      }

    }
  }
  
  free(pair);
}

/* printing bface */
static void print_bface(tBface *bface1,tBface *bface2,const char *str,struct PAIR_S *pair, int np)
{

  FLAG_T flg;

  flg = check_bface(pair,np,bface1);

  if (flg == CONTINUE_F) return;

  flg = check_bface(pair,np,bface2);

  if (flg == CONTINUE_F) return;

  FILE *fp1,*fp2;
  char fname1[100] = {0},fname2[100] = {0};// file name
  char dir1[100] = {0},dir2[100] = {0};
  char *folder = Gets("outdir");
  char *slsh = "/";
  int f1,f2,b1,b2; // face and box number

  strcpy(dir1,folder);
  strcpy(dir2,folder);
  strcat(dir1,slsh);
  strcat(dir2,slsh);

  f1 = bface1->f;
  b1 = bface1->b;

  if (bface2 != 0)
  {
    tGrid *grid = bface1->grid;
    char appn1[20],appn2[20];
    int i,m;

    appn1[0] = '\0';
    appn2[0] = '\0';

    if (bface1->touch == 1 && bface1->setnormalderiv == 1 
        && bface1->same_fpts == 1)
    {
      sprintf(appn1,"ND");
      sprintf(appn2,"NND");
    }
    else if (bface2->touch == 1 && bface2->setnormalderiv == 1
        && bface1->same_fpts == 1)
    {
      sprintf(appn2,"ND");
      sprintf(appn1,"NND");
    }
    else if (bface1->touch == 1 && bface1->same_fpts == 0)
    {
      if (bface1->sameX == 1)  strcat(appn1,"X_\0");
      if (bface1->sameY == 1)  strcat(appn1,"Y_\0");
      if (bface1->sameZ == 1)  strcat(appn1,"Z_\0");

      if (bface2->sameX == 1)  strcat(appn2,"X_\0");
      if (bface2->sameY == 1)  strcat(appn2,"Y_\0");
      if (bface2->sameZ == 1)  strcat(appn2,"Z_\0");
    }
    else
    {
      sprintf(appn2," ");
      sprintf(appn1," ");
    }


    f2 = bface2->f;
    b2 = bface2->b;
    sprintf(fname1,"b1:%d f1:%d fi1:%d_b2:%d f2:%d fi2:%d _%s_%s1",
                    b1,f1,bface1->fi,b2,f2,bface2->fi,appn1,str);
    sprintf(fname2,"b1:%d f1:%d fi1:%d_b2:%d f2:%d fi2:%d _%s_%s2",
                    b1,f1,bface1->fi,b2,f2,bface2->fi,appn2,str);
    strcat(dir1,fname1);
    strcat(dir2,fname2);
    fp1 = fopen(dir1,"w+");
    fp2 = fopen(dir2,"w+");
    assert(fp1 != 0);
    assert(fp2 != 0);


    m = bface1->fpts->npoints[bface1->b];
    assert(m > 0);
    for(i = 0; i < m; i++)
    {
      double x1[3];

      get_x_coord(x1,grid->box[bface1->b],
                  bface1->fpts->point[bface1->b][i]);
      fprintf(fp1,"%f %f %f\n",x1[0],x1[1],x1[2]);

    }

    fclose(fp1);

    m = bface2->fpts->npoints[bface2->b];
    assert(m > 0);
    for(i = 0; i < m; i++)
    {
      double x2[3];

      get_x_coord(x2,grid->box[bface2->b],
                  bface2->fpts->point[bface2->b][i]);
      fprintf(fp2,"%f %f %f\n",x2[0],x2[1],x2[2]);

    }

    fclose(fp2);
  }
  else
  {
    tGrid *grid = bface1->grid;
    int i,m;

    sprintf(fname1,"b1:%d f1:%d _%s",b1,f1,str);
    strcat(dir1,fname1);
    fp1 = fopen(dir1,"w+");

    m = bface1->fpts->npoints[bface1->b];
    assert(m > 0);
    for(i = 0; i < m; i++)
    {
      double x1[3];

      get_x_coord(x1,grid->box[bface1->b],
                  bface1->fpts->point[bface1->b][i]);

      fprintf(fp1,"%f %f %f\n",x1[0],x1[1],x1[2]);

    }

    fclose(fp1);
  }


}

/* test bfaces*/
static void test_bfaces(tGrid *grid)
{
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;

    for (bf = 0; bf < box->nbfaces; bf++)
    {
      tBface *bface = box->bface[bf];

      /* Check the pair bfaces */
      if (bface->touch == 1)
      {
        tBface *bface2 = grid->box[bface->ob]->bface[bface->ofi];

        if (grid->box[bface2->ob]->bface[bface2->ofi] != bface)
        {
          errorexit("The ob or ofi flags is not set correctly!\n");
        }
      }// if (bface->touch == 1)

      /* check the indices of pair bfaces */
      if (bface->touch == 1 && bface->same_fpts == 1)
      {
        tBface *bface2 = grid->box[bface->ob]->bface[bface->ofi];
        int nb1, nb2;

        if (bface->setnormalderiv != 1 && bface2->setnormalderiv != 1)
        {
          errorexit("Flag setnormalderiv has not been set correctly\n");
        }

        nb1 = bface->fpts->nblist;
        nb2 = bface2->fpts->nblist;

        assert(nb1 == 1);
        assert(nb2 == 1);

        int np1 = bface->fpts->npoints[bface->fpts->blist[nb1-1]];
        int np2 = bface2->fpts->npoints[bface2->fpts->blist[nb2-1]];

        if ( np1 != np2 )
        {
          errorexit("The number of points "
              "for paired bfaces are not matched!\n");
        }

        int i;
        for (i = 0; i < np1; i++)
        {
          tBox *box2 = grid->box[bface->ob];
          double x1[3];
          double x2[3];

          int *npoints = bface->fpts->npoints;
          int blist = bface->fpts->blist[0];

          get_x_coord(x1,box,bface->fpts->point[blist][i]);

          npoints = bface2->fpts->npoints;
          blist = bface2->fpts->blist[0];

          get_x_coord(x2,box2,bface2->fpts->point[blist][i]);

          if (!dequal(norm(x1,x2),0))
          {
            errorexit("The indices of points "
              "for paired bfaces are not match!\n");
          }

        }
      }// if (bface->touch == 1 && bface->same_fpts == 1)

      if (bface->touch == 0)
      {
        if (bface->ofi >= 0)
          errorexit("The untouch bface is not set correctly");
      }
    }// for (bf = 0; bf < box->nbfaces; bf++)

  }// forallboxes(grid,b)

}

/* visualize boxes */
static void visualize_boxes(tGrid *grid)
{
  FILE *fp;
  char fname[100] = {0}, dir[100] = {0};
  char *folder = Gets("outdir");
  char *slsh = "/";
  int b;



  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i;

    dir[0] = '\0';
    strcpy(dir,folder);
    strcat(dir,slsh);

    fname[0] = 0;
    sprintf(fname,"box:%d",b);
    strcat(dir,fname);

    fp = fopen(dir,"w+");
    assert(fp != 0);


    forallpoints(box,i)
    {
      double x[3];
      get_x_coord(x,box,i);

      fprintf(fp,"%f %f %f\n",x[0],x[1],x[2]);
    }


    fclose(fp);
  }

}

/*add bface to pair structure*/
static void add_to_pair(struct PAIR_S **pair,tBface *bface1,tBface *bface2,int *np)
{
  assert(bface1 != 0 || bface2 != 0);

  *pair = realloc(*pair,(*np+1)*sizeof(**pair));

  if (bface1 != 0 && bface2 != 0)
  {
    (*pair)[*np].bface1 = bface1;
    (*pair)[*np].bface2 = bface2;
  }
  else if (bface1 == 0)
  {
    (*pair)[*np].bface1 = bface2;
    (*pair)[*np].bface2 = bface2;
  }
  else
  {
    (*pair)[*np].bface1 = bface1;
    (*pair)[*np].bface2 = bface1;
  }

  ++*np;
}
