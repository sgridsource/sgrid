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
  int ind=0; /* works only if the x_of_X[i] don't use ind */

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
