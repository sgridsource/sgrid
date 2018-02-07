#include "sgrid.h"
#include "Coordinates.h"
#include "populate_bfaces_AR.h"

/*
NOTE2:
It's only work for boxes with 6 faces
*/


/*populating bface structure*/
void populate_bface(tGrid *grid)
{
  struct FACE_POINT_S ***FacePoint;//Format is FacePoint[box][face]->
  int b;
  
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
      if (FacePoint[b][f]->inner != NULL)
        free(FacePoint[b][f]->inner);
        
      if (FacePoint[b][f]->edge != NULL)
        free(FacePoint[b][f]->edge);
        
      if (FacePoint[b][f]->shared != NULL)
        free(FacePoint[b][f]->shared);
        
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


/*there are some flags which are not set yet,
  furtheremore, the indices of fpts needs to be orderly*/
static void setting_remaining_flags(tGrid *grid)
{
  
  /*set the correct ofi index if it is needed*/
  //set_ofi_flag(grid);
  set_ofi_in_all_bfaces(grid);
  
  /*Finding all of the bfaces which touching each other and need copy, 
  they are called pair, and then order the index of each ftps point*/
  order_ftps_pair(grid);
  
  /*set bit fields in bfaces for each box on the grid, 
    once, ob and ofi are known*/
  set_bits_in_all_bfaces(grid);
  
}

/*Finding all of the bfaces which 
    touching each other and need copy, they are called pair, and
    then order the index of each ftps point*/
static void order_ftps_pair(tGrid *grid)
{
  /*those bafaces which are paired*/
  struct PAIR_S *pair;
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
         
         pair = realloc(pair,(np+1)*sizeof(*pair));
         pair[np].bface1 = box->bface[bf];
         pair[np].bface2 = grid->box[box->bface[bf]->ob]->bface[box->bface[bf]->ofi];
         np++;
      }
    }
  }
  
  /*Analyzing each pair*/
  int i;
  for (i = 0; i < np; i++)
  {
    pair[i].bface1->setnormalderiv = 1;
    
    /*arranging the point index in fpts*/
    assert(pair[i].bface1->fpts->nblist == pair[i].bface2->fpts->nblist);
    
    /*For each blist temproray take it 0*/
    int bln1 = pair[i].bface1->fpts->blist[0];
    int bln2 = pair[i].bface2->fpts->blist[0];
    int m1 = pair[i].bface1->fpts->npoints[bln1];
    int m2 = pair[i].bface2->fpts->npoints[bln2];
    int j1;
    
    assert(m1 == m2);
    
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
}  


/*set the correct ofi index if it is needed*/
static void set_ofi_flag(tGrid *grid)
{
  int b;
  
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int bf;
    
    for (bf = 0; bf < box->nbfaces; bf++)
    {
      /*If they both touch and collocated*/
      if (box->bface[bf]->touch == 1 && box->bface[bf]->same_fpts == 1)
      {
        int b2;
        int bf2;
        
        forallboxes(grid,b2)
        {
          tBox *box2 = grid->box[b2];
          for (bf2 = 0; bf2 < box2->nbfaces; bf2++)
          {
            if (box2->bface[bf2]->touch == 1 && box2->bface[bf2]->same_fpts == 1)
            {
              if (box2->bface[bf2]->ofi == box->bface[bf]->f && box2->bface[bf2]->ob == box->bface[bf]->b)
              {
                box2->bface[bf2]->ofi = bf;
                box->bface[bf]->ofi   = bf2; 
                break;
              }
            }
            
          }//for (bf2 = 0; bf2 < grid->box[b2]->nbfaces; bf2++)
        }
        
      }
      /*If they only toch*/
      else if (box->bface[bf]->touch == 1)
      {
        int b2;
        int bf2;
        
        forallboxes(grid,b2)
        {
          tBox *box2 = grid->box[b2];
          for (bf2 = 0; bf2 < box2->nbfaces; bf2++)
          {
            if (box2->bface[bf2]->touch == 1)
            {
              if (box2->bface[bf2]->ofi == box->bface[bf]->f && box2->bface[bf2]->ob == box->bface[bf]->b)
              {
                box2->bface[bf2]->ofi = bf;
                box->bface[bf]->ofi   = bf2; 
                break;
              }
            }
            
          }//for (bf2 = 0; bf2 < grid->box[b2]->nbfaces; bf2++)
        }
        
      }
    }
  
  }
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
    }
    
    if (P_untouch.np > 0)
    {
      bface = make_bface(&P_untouch,grid);
      free(P_untouch.P);
    }
    
    if (P_same_fpts.np > 0)
    {
      bface = make_bface(&P_same_fpts, grid);
      bface->touch = 1;
      bface->same_fpts= 1;
      free(P_same_fpts.P);
    }
      
    if (P_unsame_fpts.np > 0)
    {
      bface = make_bface(&P_unsame_fpts, grid);
      bface->touch = 1;
      free(P_unsame_fpts.P);
    }
    
    if (P_untouch.np > 0)
    {
      bface = make_bface(&P_untouch,grid);
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
      if (P->P[i]->adjacent.touch == 1)
      {
        add_point(P1,P->P[i]);
      }
      else
      {
        add_point(P2,P->P[i]);
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
        add_point(P1,P->P[i]);
      }
      else
      {
        add_point(P2,P->P[i]);
      }
    }
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
    
    bfacei = add_empty_bface(box,f);
    
    box->bface[bfacei]->outerbound = 1;
    
    /*Filling point list for inner*/
    for (i = 0; i < FacePoint[b][f]->s; i++)
    {
      add_point_to_bface_inbox(box,bfacei,\
          FacePoint[b][f]->inner[i].geometry.ijk,f);
    }
    
    /*Filling point list for edge*/
    for (i = 0; i < FacePoint[b][f]->l; i++)
    {
      add_point_to_bface_inbox(box,bfacei,\
          FacePoint[b][f]->edge[i].geometry.ijk,f);
    }
    
    return CONTINUE_F;
  }
  
  else
    return NONE_F;
}

/*Finding the adjacent structure for edge points*/
static void find_adjacent_edge(struct FACE_POINT_S ***const FacePoint,tGrid *grid)
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
      if (FacePoint[b][f]->outerbound == 1)
        continue;
      
      //assert(FacePoint[b][f]->sh != 0);
      if (FacePoint[b][f]->sh == 0)
        printf("%d %d\n",b,f);
      
      int l = FacePoint[b][f]->l;
      struct POINT_S *edge = FacePoint[b][f]->edge;
      int *blist = malloc(FacePoint[b][f]->sh*sizeof(*blist));
      int i;
      
      for (i = 0; i < FacePoint[b][f]->sh; i++)
      {
        blist[i] = FacePoint[b][f]->shared[i].box;
      }
      
      for (i = 0; i < l; i++)
      {
        double X[3];
        double *x = FacePoint[b][f]->edge[i].geometry.x;
        
        int adjacent_box = b_XYZ_of_xyz_inboxlist(grid,blist,FacePoint[b][f]->sh,X,X+1,X+2,x[0],x[1],x[2]);
        
        if (adjacent_box < 0) //Find the adjacent box;
        {
          blist = realloc(blist,(FacePoint[b][f]->sh+1)*sizeof(*blist));
          blist[FacePoint[b][f]->sh] = b;
          adjacent_box = b_xyz_in_exblist(grid,x,blist,FacePoint[b][f]->sh+1);
          
          if (adjacent_box < 0)//if still couldn't be found
          {
            fprintf(stderr,"ERROR:\n There is no box for this point!\n");
            abort();
          }
        }
        
        populate_adjacent(FacePoint,EDGE_F,b,f,i,adjacent_box,grid);
        
      }
      
      free (blist);
    }//for (f = 0; f < f_max; f++)
  }//for (b = 0; b < b_max; b++)
}

/*Finding the adjacent structure for inner points*/
static void find_adjacent_inner(struct FACE_POINT_S ***const FacePoint,tGrid *grid)
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
      int s = FacePoint[b][f]->s;
      struct POINT_S *inner = FacePoint[b][f]->inner;
      int i;
      
      for (i = 0; i < s; i++)
      {
        double q[3] = { inner[i].geometry.x[0]+inner[i].geometry.N[0]*EPS,\
                        inner[i].geometry.x[1]+inner[i].geometry.N[1]*EPS,\
                        inner[i].geometry.x[2]+inner[i].geometry.N[2]*EPS};
        double X[3];
        int adjacent_box = b_XYZ_of_xyz(grid,X,X+1,X+2,q[0],q[1],q[2]);
        
        if (adjacent_box < 0)//the face is considered as outer or inner bound so the loop is terminated
        {
          FacePoint[b][f]->outerbound = 1;
          break;
        }
        else
        {
          populate_adjacent(FacePoint,INNER_F,b,f,i,adjacent_box,grid);
        }
          
        
      }//for (i = 0; i < s; i++)
    }//for (f = 0; f < f_max; f++)
  }//for (b = 0; b < b_max; b++)
}

/*given x finding the pertinent box excluding ex_blist*/
static int b_xyz_in_exblist(tGrid *grid,double *x,int *ex_blist,int ex_bn)
{
  int b, found,*b_list,j;
  double X[3];
  
  forallboxes(grid,b);
  
  const int k = b-ex_bn;
  
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
  
  b = b_XYZ_of_xyz_inboxlist(grid,b_list,k,X,X+1,X+2,x[0],x[1],x[2]);
  
  return b;
}

/*adding info to share structure*/
static void add_info(struct FACE_POINT_S *const FacePoint,int f, int b, struct POINT_S *P)
{
  const int sh = FacePoint->sh;
  
  /*Make sure the identical box and face won't be double counted*/
  if (sh != 0)
  {
    int i = 0;
    
    while (i < sh)
    {
      if ( FacePoint->shared[i].box == b && FacePoint->shared[i].face == f )
      {
        FacePoint->shared[i].P = realloc(FacePoint->shared[i].P,\
                (FacePoint->shared[i].np+1)*sizeof(*FacePoint->shared[i].P));
        FacePoint->shared[i].P[FacePoint->shared[i].np] = P;
        FacePoint->shared[i].np++;
        return;
      }
        
      i++;
    }
  }
  
  FacePoint->shared = realloc(FacePoint->shared,\
                  (FacePoint->sh+1)*sizeof(*FacePoint->shared));
                  
  FacePoint->shared[sh].face = f;
  FacePoint->shared[sh].box = b;
  FacePoint->shared[sh].P = malloc(sizeof(*FacePoint->shared[sh].P));
  FacePoint->shared[sh].P[0] = P;
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
static int IsAtFace(int box_n, int *ijk,int *f,tGrid *grid)
{
  tBox *box = grid->box[box_n];
  int n[3] = {box->n1,box->n2,box->n3};
  unsigned int flg = 0;
  int i;
  
  for (i = 0; i < TOT_NUM_FACE; i++)
    f[i] = -1;
  
  i = 0;
  if (ijk[0] == 0)
  {
    f[i++] = X_FACE0;
    flg = 1;
  }
  if (ijk[0] == n[0]-1)
  {
    f[i++] = X_FACE1;
    flg = 1;
  }
  if (ijk[1] == 0)
  {
    f[i++] = Y_FACE0;
    flg = 1;
  }
  if (ijk[1] == n[1]-1)
  {
    f[i++] = Y_FACE1;
    flg = 1;
  }
  if (ijk[2] == 0)
  {
    f[i++] = Z_FACE0;
    flg = 1;
  }
  if (ijk[2] == n[2]-1)
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
  
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n[3] = {box->n1,box->n2,box->n3};
    int n1 = n[0] ,n2 = n[1];//They are needed for macro Index
    
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
      
      /*Filling the edge and inner points*/
      setup_range(f,n,&i_l,&j_l,&k_l,&i_u,&j_u,&k_u,NULL);
      e = 0;
      in = 0;
      for (i = i_l; i <= i_u; i++)
      {
        for (j = j_l; j <= j_u; j++)
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
           get_normal(P->geometry.N,box,f,Index(i,j,k));
           P->geometry.ijk = Index(i,j,k);
           P->geometry.b = b;
           P->geometry.f = f;
           
         }
      }
    }
  }//End of forallboxes
}

/*Gettig the normalized normal in outward direction*/
static void get_normal(double *N,tBox *box,int face,long int ijk)
{
  int idXd = Ind("dXdx");
  int idYd = Ind("dYdx");
  int idZd = Ind("dZdx");
  int j;

  if(box->v[idXd]==NULL)//If the box using cartesian coords
  {
    if (face == X_FACE0 || face == X_FACE1)
    {
      N[0] = 1;
      N[1] = 0;
      N[2] = 0;
    }
    else if (face == Y_FACE0 || face == Y_FACE1)
    {
      N[0] = 0;
      N[1] = 1;
      N[2] = 0;
    }
    else if (face == Z_FACE0 || face == Z_FACE1)
    {
      N[0] = 0;
      N[1] = 0;
      N[2] = 1;
    }
  }
  
  else
  {
    if (face == X_FACE0 || face == X_FACE1)
    {
      for(j = 0; j < 3; j++)
        N[j] = box->v[idXd+j][ijk];
    }
    else if (face == Y_FACE0 || face == Y_FACE1)
    {
      for(j = 0; j < 3; j++)
        N[j] = box->v[idYd+j][ijk];
    }
    else if (face == Z_FACE0 || face == Z_FACE1)
    {
      for(j = 0; j < 3; j++)
        N[j] = box->v[idZd+j][ijk];
    }
    
    normalizing_N(N);
    
  }
  
  /*Setting the sign in order to have outward direction*/
  int s = 2*(face%2)-1;
  N[0] *= s;
  N[1] *= s;
  N[2] *= s;
  
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
static void populate_adjacent(struct FACE_POINT_S ***const FacePoint,FLAG_T kind,int b, int f, long int i, int adjacent_box, tGrid *grid)
{
  struct POINT_S *P;
  
  int ijk_adj[3];//the index of the point in adjacent box
  int face_list[TOT_NUM_FACE];
  
  if (kind == INNER_F)
  {
    P = &FacePoint[b][f]->inner[i];
  }
  else if (kind == EDGE_F)
  {
    P = &FacePoint[b][f]->edge[i];
  }
  
  P->adjacent.box = adjacent_box;
  
  /*If the point x collocated with the 
     same point in adjacent point*/
  if (IsCollocated(P->geometry.x,adjacent_box,ijk_adj,P->adjacent.X,grid) == 1)
  {
    P->adjacent.ijk[0] = ijk_adj[0];
    P->adjacent.ijk[1] = ijk_adj[1];
    P->adjacent.ijk[2] = ijk_adj[2];
    
    /*If the found point is at face*/
    FLAG_T flg = ERROR_F;
    
    if (IsAtFace(adjacent_box,ijk_adj,face_list,grid) == 1)
    {
      double N[3];
      int i = 0;
      
      P->adjacent.touch = 1;
      
      i = 0;
      while (face_list[i] >= 0)
      {
        int n1 = grid->box[adjacent_box]->n1;
        int n2 = grid->box[adjacent_box]->n2;
        int n3 = grid->box[adjacent_box]->n3;
        
        get_normal(N,grid->box[adjacent_box],\
          face_list[i],ijk_ind(grid->box[adjacent_box],ijk_adj));
        
        if (dequal(ABS(dot_product(P->geometry.N,N)),1))
        {
          P->adjacent.face = face_list[i];
          P->adjacent.interpolation = 0;
          add_info(FacePoint[b][f],face_list[i],adjacent_box,P);
          
          flg = NONE_F;
          break;
        }
        i++;
      }
      
      if (flg == ERROR_F)//ERROR
      {
        fprintf(stderr,"ERROR:\n There is a collocation point in box: %d and at face: %d,\n"
          "but its normal isn't match!\n",b,f);
        abort();
      }
    }
    
  }
  else/*Since the found point is not collocated it needs interpolation*/
  {
    int face[TOT_NUM_FACE];
    int nf = 0,i;
    
    nf = 
        XYZ_on_face(grid->box[adjacent_box],face,P->geometry.X[0],P->geometry.X[1],P->geometry.X[2]);
    
    if (nf == 0)
      P->adjacent.touch = 0;//doesn't matter it is touching or not
    else
    {
      for (i = 0; i < TOT_NUM_FACE; i++)
      {
        if ( face[i] == 1 )
        {
          double N[3];
          
          get_normal(N,grid->box[adjacent_box],i,ijk_ind(grid->box[adjacent_box],ijk_adj));
          
          if (dequal(ABS(dot_product(P->geometry.N,N)),1))
          P->adjacent.face = i;
        }
      }
    }
      
    P->adjacent.interpolation = 1;
    P->adjacent.touch = 1;
    
  }
}

/*Returning the index at box*/
static int ijk_ind(tBox *box, int *i)
{
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  
  return Index(i[0],i[1],i[2]);
}