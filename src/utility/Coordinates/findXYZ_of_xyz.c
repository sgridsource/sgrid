/* findXYZ_of_xyz.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Coordinates.h"


/* global vars needed to pass info into xyz_VectorFunc */
double desired_x;
double desired_y;
double desired_z;
tBox *box_for_xyz_VectorFunc;


/* funtion to be passed into newton_linesrch_its by XYZ_of_xyz */
void xyz_VectorFunc(int n, double *XYZvec, double *fvec)
{
  double xg,yg,zg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  xg = box_for_xyz_VectorFunc->x_of_X[1]((void *) box_for_xyz_VectorFunc,
                                         ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  yg = box_for_xyz_VectorFunc->x_of_X[2]((void *) box_for_xyz_VectorFunc,
                                         ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  zg = box_for_xyz_VectorFunc->x_of_X[3]((void *) box_for_xyz_VectorFunc,
                                         ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  fvec[1] = xg-desired_x;
  fvec[2] = yg-desired_y;
  fvec[3] = zg-desired_z;
}

/* find X,Y,Z from x,y,z (Note: X,Y,Z also contains initial guess) */
int XYZ_of_xyz(tBox *box, double *X, double *Y, double *Z,
                double x, double y, double z)
{
  double XYZvec[4];
  int check, stat;

  if(box->x_of_X[1]==NULL)
    { *X = x;  *Y = y;  *Z = z;  return 0; }

  box_for_xyz_VectorFunc = box;
  desired_x = x;
  desired_y = y;
  desired_z = z;

  /* bad initial guess */
  /* XYZvec[1] = (box->bbox[0] + box->bbox[1]) * 0.5;
     XYZvec[2] = (box->bbox[2] + box->bbox[3]) * 0.5;
     XYZvec[3] = (box->bbox[4] + box->bbox[5]) * 0.5; */

  /* initial guess supplied by caller */
  XYZvec[1] = *X;
  XYZvec[2] = *Y;
  XYZvec[3] = *Z;
      
  /* do newton_linesrch_its iterations: */
  stat = newton_linesrch_its(XYZvec, 3, &check, xyz_VectorFunc, 
                             Geti("Coordinates_newtMAXITS"),
                             Getd("Coordinates_newtTOLF") );
  *X = XYZvec[1];
  *Y = XYZvec[2];
  *Z = XYZvec[3];

  if(check || stat<0) printf("XYZ_of_xyz: check=%d stat=%d\n", check, stat);
  return stat-check;
}

/* find X,Y,Z from x,y,z (Note: X,Y,Z also contains initial guess)
   return index of box in which X,Y,Z are found */
int b_XYZ_of_xyz(tGrid *grid, double *X, double *Y, double *Z,
                 double x, double y, double z)
{
  double X1,Y1,Z1;
  int stat, bi;

  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];

    X1=*X; Y1=*Y; Z1=*Z;
    stat = XYZ_of_xyz(box, &X1,&Y1,&Z1, x,y,z);
    if(stat<0) continue;
    if(dless(X1,box->bbox[0]) || dless(box->bbox[1],X1)) continue;
    if(dless(Y1,box->bbox[2]) || dless(box->bbox[3],Y1)) continue;
    if(dless(Z1,box->bbox[4]) || dless(box->bbox[5],Z1)) continue;
    /* round X1,Y1,Z1 inside box */
    if(X1 < box->bbox[0]) X1 = box->bbox[0];
    if(X1 > box->bbox[1]) X1 = box->bbox[1];
    if(Y1 < box->bbox[2]) Y1 = box->bbox[2];
    if(Y1 > box->bbox[3]) Y1 = box->bbox[3];
    if(Z1 < box->bbox[4]) Z1 = box->bbox[4];
    if(Z1 > box->bbox[5]) Z1 = box->bbox[5];
    *X=X1; *Y=Y1; *Z=Z1;
    return bi; /* return box index if success */
  }
  return -bi;
}

/* find X,Y,Z from x,y,z (Note: X,Y,Z also contains initial guess)
   search in all boxes in the list: bi = blist[0], blist[1], ..., blist[nb-1]
   return index of box in which X,Y,Z are found or erro code if failure */
int b_XYZ_of_xyz_inboxlist(tGrid *grid, int *blist, int nb,
                           double *X, double *Y, double *Z,
                           double x, double y, double z)
{
  double X1,Y1,Z1;
  int stat, bi, i;
  tBox *box;

  for(i=0; i<nb; i++)
  {
    bi = blist[i];
    box = grid->box[bi];

    X1=*X; Y1=*Y; Z1=*Z;
    stat = XYZ_of_xyz(box, &X1,&Y1,&Z1, x,y,z);
    if(stat<0) continue;
    if(dless(X1,box->bbox[0]) || dless(box->bbox[1],X1)) continue;
    if(dless(Y1,box->bbox[2]) || dless(box->bbox[3],Y1)) continue;
    if(dless(Z1,box->bbox[4]) || dless(box->bbox[5],Z1)) continue;
    /* round X1,Y1,Z1 inside box */
    if(X1 < box->bbox[0]) X1 = box->bbox[0];
    if(X1 > box->bbox[1]) X1 = box->bbox[1];
    if(Y1 < box->bbox[2]) Y1 = box->bbox[2];
    if(Y1 > box->bbox[3]) Y1 = box->bbox[3];
    if(Z1 < box->bbox[4]) Z1 = box->bbox[4];
    if(Z1 > box->bbox[5]) Z1 = box->bbox[5];
    *X=X1; *Y=Y1; *Z=Z1;
    return bi; /* return box index if success */
  }
  return -bi;
}

/* find nearest X,Y,Z in box from x,y,z */
double nearestXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                         double x, double y, double z)
{
  int i;
  double *pX = box->v[Ind("X")];
  double *pY = box->v[Ind("Y")];
  double *pZ = box->v[Ind("Z")];
  double *px = box->v[Ind("x")];
  double *py = box->v[Ind("y")];
  double *pz = box->v[Ind("z")];
  double dxi, dyi, dzi, r, rmin=-1.0;

  forallpoints(box, i)
  {
    dxi = px[i] - x;
    dyi = py[i] - y;
    dzi = pz[i] - z;
    r = dxi*dxi + dyi*dyi + dzi*dzi;
    if(r<=rmin || rmin<0.0 )
    {
      rmin = r;
      *ind = i;
      *X = pX[i];
      *Y = pY[i];
      *Z = pZ[i];
    } 
  }
  return rmin;
}

/* find nearest X,Y,Z on grid from x,y,z */
double nearest_b_XYZ_of_xyz(tGrid *grid,  int *b, int *ind,
                            double *X, double *Y, double *Z,
                            double x, double y, double z)
{
  double r, rmin=-1.0;
  double X1,Y1,Z1;
  int ind1, bi;

  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];

    r = nearestXYZ_of_xyz(box, &ind1, &X1,&Y1,&Z1, x,y,z);
    if(r<=rmin || rmin<0.0 )
    {
      rmin=r;  *b=bi;  *ind=ind1;   *X=X1;  *Y=Y1;  *Z=Z1;
    }
  }
  return rmin;
}

/* find nearest X,Y,Z in grid-plane from x,y,z */
double nearestXYZ_of_xyz_inplane(tBox *box, int *ind, 
                                 double *X, double *Y, double *Z,
                                 double x, double y, double z,
                                 int plane, int pind)
{
  int i,j,k, in;
  double *pX = box->v[Ind("X")];
  double *pY = box->v[Ind("Y")];
  double *pZ = box->v[Ind("Z")];
  double *px = box->v[Ind("x")];
  double *py = box->v[Ind("y")];
  double *pz = box->v[Ind("z")];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double dx, dy, dz, r, rmin=-1.0;

  if(box->x_of_X[1]==NULL)
    errorexit("nearestXYZ_of_xyz_inplane: box->x_of_X[1] = NULL");

  if(plane==1)
    forplane1(i,j,k, n1,n2,n3, pind)
    {
      in = Index(i,j,k);
      dx = px[in] - x;
      dy = py[in] - y;
      dz = pz[in] - z;
      r = dx*dx + dy*dy + dz*dz;
      if(r<=rmin || rmin<0.0 )
      {
        rmin = r;
        *ind = in;
        *X = pX[in];
        *Y = pY[in];
        *Z = pZ[in];
      } 
    }
  else if(plane==2)
    forplane2(i,j,k, n1,n2,n3, pind)
    {
      in = Index(i,j,k);
      dx = px[in] - x;
      dy = py[in] - y;
      dz = pz[in] - z;
      r = dx*dx + dy*dy + dz*dz;
      if(r<=rmin || rmin<0.0 )
      {
        rmin = r;
        *ind = in;
        *X = pX[in];
        *Y = pY[in];
        *Z = pZ[in];
      } 
    }
  else if(plane==3)
    forplane3(i,j,k, n1,n2,n3, pind)
    {
      in = Index(i,j,k);
      dx = px[in] - x;
      dy = py[in] - y;
      dz = pz[in] - z;
      r = dx*dx + dy*dy + dz*dz;
      if(r<=rmin || rmin<0.0 )
      {
        rmin = r;
        *ind = in;
        *X = pX[in];
        *Y = pY[in];
        *Z = pZ[in];
      } 
    }
  else
    errorexit("nearestXYZ_of_xyz_inplane: plane has to be 1, 2 or 3");

  return rmin;
}

/* function to be passed into newton_linesrch_its by X_of_x_forgiven_YZ */
void x_VectorFunc_YZ(int n, double *XYZvec, double *fvec)
{
  double xg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  xg = box_for_xyz_VectorFunc->x_of_X[1]((void *) box_for_xyz_VectorFunc,
                                          ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  fvec[1] = xg-desired_x;
}

/* find X from x for a given Y,Z (Note: X also contains initial guess) */
int X_of_x_forgiven_YZ(tBox *box, double *X, double x, double Y, double Z)
{
  double XYZvec[4];
  int check, stat;

  if(box->x_of_X[1]==NULL)  { *X = x;  return 0; }

  box_for_xyz_VectorFunc = box;
  desired_x = x;

  /* initial guess supplied by caller */
  XYZvec[1] = *X;
  
  /* Y,Z are fixed */
  XYZvec[2] = Y;
  XYZvec[3] = Z;
      
  /* do newton_linesrch_its iterations: */
  stat = newton_linesrch_its(XYZvec, 1, &check, x_VectorFunc_YZ, 
                             Geti("Coordinates_newtMAXITS"),
                             Getd("Coordinates_newtTOLF") );
  *X = XYZvec[1];

  if(check || stat<0) printf("XYZ_of_xyz: check=%d stat=%d\n", check, stat);
  return stat-check;
}

/* find X from x for a given Y,Z (Note: X also contains initial guess)
   return index of box in which X,Y,Z are found */
int b_X_of_x_forgiven_YZ(tGrid *grid, double *X, double x, double Y, double Z)
{
  double X1;
  int stat, bi;

  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];

    X1=*X;
    stat = X_of_x_forgiven_YZ(box, &X1, x, Y,Z);
//printf("bi=%d: stat=%d X1=%g Y=%g Z=%g\n", bi, stat, X1, Y,Z);
//printf("%g %g   %g %g   %g %g\n",
//box->bbox[0],box->bbox[1], box->bbox[2],box->bbox[3], box->bbox[4],box->bbox[5]);
    if(stat<0) continue;
    if(dless(X1,box->bbox[0]) || dless(box->bbox[1],X1)) continue;
    if(dless(Y,box->bbox[2]) || dless(box->bbox[3],Y)) continue;
    if(dless(Z,box->bbox[4]) || dless(box->bbox[5],Z)) continue;
    /* round X1 inside box */
    if(X1 < box->bbox[0]) X1 = box->bbox[0];
    if(X1 > box->bbox[1]) X1 = box->bbox[1];
    *X=X1;
    return bi; /* return box index if success */
  }
  return -bi;
}
