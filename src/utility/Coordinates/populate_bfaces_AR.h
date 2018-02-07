#include <assert.h>
#define SQ(x) (x)*(x)
#define EPS 1E-7

/*Info about the geometry of a point*/
struct GEOMETRY_S
{
  double x[3];//x_coord
  double X[3];//X_coord
  double N[3];//Normal vector
  int b;//box the point in it
  int f;//face the point on it
  int ijk;
};

/*Info about the adjacent of a point*/
struct ADJACENT_S
{
  double X[3];//X_coord
  int box;
  int face;
  int ijk[3];
  int oXi[3];
  unsigned int touch: 1;//1 if it touches, 0 if it overlaps
  unsigned int interpolation: 1;//1 if needs interpolation, 0 otherwise
};

/*All info about a point*/
struct POINT_S
{
  struct GEOMETRY_S geometry;
  struct ADJACENT_S adjacent;
};

/*Shared info*/
struct SHARED_S
{
  int face;
  int box;
  struct POINT_S **P;/*The colloction of point which are 
        in a same box and on same face and 
          having the same adjacent box and face*/
  int np;//total num of points above
};

/*Face points*/
struct FACE_POINT_S
{
  /*all points at a face but edge points*/
  struct POINT_S *inner;
  int s;//total num of inner points
  
  /*all edge points in a face*/
  struct POINT_S *edge;
  int l;//total num of edge points
  
  struct SHARED_S *shared;
  int sh;//total num of shared
  
  unsigned int outerbound;//if the face is outerbound it is 1, otherwise 0
};

/*Similiar points*/
struct SIMILAR_S
{
  struct POINT_S **P;
  int np;
};

/*those bafaces which are paired*/
struct PAIR_S
{
  tBface *bface1;
  tBface *bface2;
};

/*Convention for index of each face*/
enum FACE
{
  X_FACE0 = 0,/*when i = 0*/
  X_FACE1,/*when i = n1-1*/
  Y_FACE0,/*when j = 0*/
  Y_FACE1,/*when j = n2-1*/
  Z_FACE0,/*when k = 0*/
  Z_FACE1,/*when k = n3-1*/
  TOT_NUM_FACE/*Total number of faces*/
};

typedef enum FLAG
{
    INNER_F,
    EDGE_F,
    CONTINUE_F,
    NONE_F,
    STOP_F,
    TOUCH_F,
    SAME_FTPS_F,
    ERROR_F
    
}FLAG_T;


/*Functions*/
static void *allc_FacePoint(tGrid *grid);
static void init_FacePoint(struct FACE_POINT_S ***FacePoint,tGrid *grid);
static void fill_geometry(struct FACE_POINT_S ***FacePoint,tGrid *grid);
static int tot_num(void *p);
static void get_x_coord(double *x,tBox *box,long int ind);
static void get_X_coord(double *X,tBox *box,long int ind);
static void setup_range(int face,int *n,int *i_l,int *j_l,int *k_l,int *i_u,int *j_u,int *k_u, int *m);
static int IsEdge(int face, int *n, int i, int j, int k);
static void get_normal(double *N,tBox *box,int face,long int ijk);
static void normalizing_N(double *N);
static void find_adjacent_inner(struct FACE_POINT_S ***const FacePoint,tGrid *grid);
static void find_adjacent_edge(struct FACE_POINT_S ***const FacePoint,tGrid *grid);
static double norm(double *x1,double *x2);
static int IsCollocated(double *x1, int box2, int *ijk_adj,double *X,tGrid *grid);
static int IsAtFace(int box, int *ijk,int *f,tGrid *grid);
static double ABS(double x);
static double dot_product(double *x,double *y);
static void add_info(struct FACE_POINT_S *const FacePoint,int f, int b,struct POINT_S *P);
static int b_xyz_in_exblist(tGrid *grid,double *x,int *ex_blist,int ex_bn);
static void populate_adjacent(struct FACE_POINT_S ***const FacePoint,FLAG_T kind,int b, int f, long int i, int adjacent_box, tGrid *grid);
static void group_similar_points(struct FACE_POINT_S ***const FacePoint,tGrid *grid);
static int find_outerbound_bfaces(struct FACE_POINT_S ***const FacePoint,int b, int f,tGrid *grid);
static void find_remaining_bfaces(struct FACE_POINT_S ***const FacePoint,int b, int f,tGrid *grid);
static tBface *make_bface(struct SIMILAR_S *P, tGrid *grid);
static void add_point(struct SIMILAR_S *hay,struct POINT_S *needle);
static FLAG_T check_bface(struct PAIR_S *pair,int np, tBface *bface);
static void setting_remaining_flags(tGrid *grid);
static void group_point(FLAG_T kind,struct SIMILAR_S *P1, struct SIMILAR_S *P2, void *input);
static void free_FacePoint(struct FACE_POINT_S ***FacePoint);
static void set_ofi_flag(tGrid *grid);
static void order_ftps_pair(tGrid *grid);
static int ijk_ind(tBox *box, int *i);
int b_XYZ_of_xyz(tGrid *grid, double *X, double *Y, double *Z,double x, double y, double z);
int b_XYZ_of_xyz_inboxlist(tGrid *grid, int *blist, int nb,double *X, double *Y, double *Z,double x, double y, double z);
int add_empty_bface(tBox *box, int f);
int set_bits_in_all_bfaces(tGrid *grid);
int set_ofi_in_all_bfaces(tGrid *grid);
int XYZ_on_face(tBox *box, int *face, double X, double Y, double Z);
void AddToPointList(tPointList *PL, int boxindex, int newpoint);