/* findXYZ_of_xyz.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Coordinates.h"


/* global vars needed to pass info into xyz_VectorFunc */
double desired_x;
double desired_y;
double desired_z;
tBox *box_for_xyz_VectorFunc;


/* funtion to be passed into newton_lnsrch by XYZ_of_xyz */
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
void XYZ_of_xyz(tBox *box, double *X, double *Y, double *Z,
                double x, double y, double z)
{
  double XYZvec[4];
  int check;

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
      
  /* do newton_lnsrch iterations: */
  newton_lnsrch(XYZvec, 3, &check, xyz_VectorFunc, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  *X = XYZvec[1];
  *Y = XYZvec[2];
  *Z = XYZvec[3];

  if(check) printf("XYZ_of_xyz: check=%d\n", check);  
}

/* find nearest X,Y,Z on grid from x,y,z */
double nearestXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                         double x, double y, double z)
{
  int i;
  double *pX = box->v[Ind("X")];
  double *pY = box->v[Ind("Y")];
  double *pZ = box->v[Ind("Z")];
  double dxi, dyi, dzi, r, rmin=-1.0;

  if(box->x_of_X[1]==NULL)
    errorexit("nearestXYZ_of_xyz: box->x_of_X[1] = NULL");

  forallpoints(box, i)
  {
    dxi = box->x_of_X[1]((void *) box, i, pX[i], pY[i], pZ[i]) - x;
    dyi = box->x_of_X[2]((void *) box, i, pX[i], pY[i], pZ[i]) - y;
    dzi = box->x_of_X[3]((void *) box, i, pX[i], pY[i], pZ[i]) - z;
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
      dx = box->x_of_X[1]((void *) box, in, pX[in],pY[in],pZ[in]) - x;
      dy = box->x_of_X[2]((void *) box, in, pX[in],pY[in],pZ[in]) - y;
      dz = box->x_of_X[3]((void *) box, in, pX[in],pY[in],pZ[in]) - z;
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
      dx = box->x_of_X[1]((void *) box, in, pX[in],pY[in],pZ[in]) - x;
      dy = box->x_of_X[2]((void *) box, in, pX[in],pY[in],pZ[in]) - y;
      dz = box->x_of_X[3]((void *) box, in, pX[in],pY[in],pZ[in]) - z;
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
      dx = box->x_of_X[1]((void *) box, in, pX[in],pY[in],pZ[in]) - x;
      dy = box->x_of_X[2]((void *) box, in, pX[in],pY[in],pZ[in]) - y;
      dz = box->x_of_X[3]((void *) box, in, pX[in],pY[in],pZ[in]) - z;
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

/* function to be passed into newton_lnsrch by X_of_x_forgiven_YZ */
void x_VectorFunc_YZ(int n, double *XYZvec, double *fvec)
{
  double xg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  xg = box_for_xyz_VectorFunc->x_of_X[1]((void *) box_for_xyz_VectorFunc,
                                          ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  fvec[1] = xg-desired_x;
}

/* find X from x for a given Y,Z (Note: X also contains initial guess) */
void X_of_x_forgiven_YZ(tBox *box, double *X, double x, double Y, double Z)
{
  double XYZvec[4];
  int check;

  box_for_xyz_VectorFunc = box;
  desired_x = x;

  /* initial guess supplied by caller */
  XYZvec[1] = *X;
  
  /* Y,Z are fixed */
  XYZvec[2] = Y;
  XYZvec[3] = Z;
      
  /* do newton_lnsrch iterations: */
  newton_lnsrch(XYZvec, 1, &check, x_VectorFunc_YZ, 
 		Geti("Coordinates_newtMAXITS"),
    		Getd("Coordinates_newtTOLF") );
  *X = XYZvec[1];

  if(check) printf("X_of_x_forgiven_YZ: check=%d\n", check);  
}
