/* findXYZ_of_xyz.c */
/* Wolfgang Tichy 8/2007 */

#include "sgrid.h"
#include "Coordinates.h"


/* pars needed to pass info into xyz_VectorFuncP */
typedef struct T_grid_box_desired_xyz_struct {
  tGrid *grid; /* grid */
  tBox *box;   /* box */
  double desired_x;
  double desired_y;
  double desired_z;
} t_grid_box_desired_xyz_struct;   


/* funtion to be passed into newton_linesrch_itsP by XYZ_of_xyz */
void xyz_VectorFuncP(int n, double *XYZvec, double *fvec, void *p)
{
  t_grid_box_desired_xyz_struct *pars;
  tBox *box;
  double desired_x; 
  double desired_y;
  double desired_z;
  double xg,yg,zg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  /* get pars */
  pars = (t_grid_box_desired_xyz_struct *) p;
  box = pars->box;
  desired_x = pars->desired_x;
  desired_y = pars->desired_y;
  desired_z = pars->desired_z;

  xg = box->x_of_X[1]((void *) box, ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  yg = box->x_of_X[2]((void *) box, ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  zg = box->x_of_X[3]((void *) box, ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  fvec[1] = xg-desired_x;
  fvec[2] = yg-desired_y;
  fvec[3] = zg-desired_z;
  //printf("X,Y,Z=%g,%g,%g\n", XYZvec[1],XYZvec[2],XYZvec[3]);
  //printf("%g,%g,%g  xg,yg,zg=%g,%g,%g\n", fvec[1],fvec[2],fvec[3], xg,yg,zg);
}

/* find X,Y,Z from x,y,z (Note: X,Y,Z also contains initial guess) */
int XYZ_of_xyz(tBox *box, double *X, double *Y, double *Z,
                double x, double y, double z)
{
  double tol = Getd("Coordinates_newtTOLF");
  double XYZvec[4];
  //double fvec[4];
  //double err, r;
  t_grid_box_desired_xyz_struct pars[1];
  tSingInfo si[1];
  int check, stat;

  if(box->x_of_X[1]==NULL)
    { *X = x;  *Y = y;  *Z = z;  return 0; }

  pars->box = box;
  pars->desired_x = x;
  pars->desired_y = y;
  pars->desired_z = z;

  /* test if X,Y,Z are at singularity */
  stat=0;
  if(box->isSing!=NULL)  stat = box->isSing((void *) box, *X,*Y,*Z, 1, si);
  if(stat) /* if yes, try something: */
  {
    double fvec[4];
    double err, r;

    stat = recover_if_start_on_singularity(box, X,Y,Z, x,y,z, si);
    /* check if error is ok */
    XYZvec[1] = *X;
    XYZvec[2] = *Y;
    XYZvec[3] = *Z;
    xyz_VectorFuncP(3, XYZvec, fvec, (void *) pars);
    err = sqrt(fvec[1]*fvec[1] + fvec[2]*fvec[2] + fvec[3]*fvec[3]);
    r = sqrt(x*x + y*y + z*z);
    if(r>0.0)  err = err/r;
    if(err>tol*10.0 && stat>=0)
    {
      printf("XYZ_of_xyz: recover_if_start_on_singularity failed: err=%g\n",
             err);
      printf("X=%g Y=%g Z=%g\n", *X,*Y,*Z);
      stat = -stat;
    }
    return stat;
  }

  /* bad initial guess */
  /* XYZvec[1] = (box->bbox[0] + box->bbox[1]) * 0.5;
     XYZvec[2] = (box->bbox[2] + box->bbox[3]) * 0.5;
     XYZvec[3] = (box->bbox[4] + box->bbox[5]) * 0.5; */

  /* initial guess supplied by caller */
  XYZvec[1] = *X;
  XYZvec[2] = *Y;
  XYZvec[3] = *Z;

  /* do newton_linesrch_itsP iterations: */
  stat = newton_linesrch_itsP(XYZvec, 3, &check, xyz_VectorFuncP, (void *) pars,
                              Geti("Coordinates_newtMAXITS"), tol);
  *X = XYZvec[1];
  *Y = XYZvec[2];
  *Z = XYZvec[3];

  //if(check || stat<0) printf("XYZ_of_xyz: check=%d stat=%d\n", check, stat);

//  /* get error in x,y,z */
//  xyz_VectorFuncP(3, XYZvec, fvec, (void *) pars);
//  err = sqrt(fvec[1]*fvec[1] + fvec[2]*fvec[2] + fvec[3]*fvec[3]);
//  r = sqrt(x*x + y*y + z*z);
//  if(r>0.0)  err = err/r;
//  if(err>1e-5 && stat>=0)
//  {
//    //printf("XYZ_of_xyz: newton_linesrch_itsP failed! err=%g\n", err);
//    stat = -stat;
//  }

  if(check || stat<0)
  {
    //printf("XYZ_of_xyz: check=%d stat=%d err=%g\n", check, stat, err);
    printf("XYZ_of_xyz: check=%d stat=%d\n", check, stat);
    printf("            in box%d at x=%g y=%g z=%g\n", box->b, x,y,z);
  }

  return stat-check;
}

/* do something if we start on a coordinate singularity */
int recover_if_start_on_singularity(tBox *box,
        double *X, double *Y, double *Z,
        double x, double y, double z, tSingInfo *si)
{
  int stat;
  int face[6]; /* face where sing. is located, e.g. face[2]=1 => face2 */
  int dir[4];  /* direction info */
  int zc[4];   /* cols with zeros, e.g. zc[3]=1 => col3 has all zeros */
  int nf, i, j;

  /* find faces */
  nf = XYZ_on_face(box, face, *X,*Y,*Z);
  /* remove face in periodic dirs */
  for(i=0; i<6; i++)
    if(face[i] && box->periodic[1+i/2]) { face[i] = 0; nf--; }

  /* find dirs with sing. */
  for(i=1; i<=3; i++) dir[i]=0;
  if(face[0] || face[1])      dir[1]=1;
  else if(face[2] || face[3]) dir[2]=1;
  else if(face[4] || face[5]) dir[3]=1;
  else errorexit("singularity should be on one of the faces");

  /* check which cols in dx_dX have only zeros */
  zc[1] = zc[2] = zc[3] = 1;
  for(j=1; j<=3; j++)
    for(i=1; i<=3; i++)  zc[j] = zc[j] && si->dx_dX[i][j] == '0';

  /* catch some cases (so far only the ones for AnsorgNS) */
  if(1)
  {
    if(dir[1])
    {
      prSingInfo(si);
      printf("X=%g Y=%g Z=%g\n", *X,*Y,*Z);
      for(i=0; i<6; i++) printf("face[%d]=%d ", i, face[i]);
      printf(" nf=%d\n", nf);
      printf(" dir[1]=%d dir[2]=%d dir[3]=%d\n", dir[1], dir[2], dir[3]);
      printf("  zc[1]=%d  zc[2]=%d  zc[3]=%d\n", zc[1], zc[2], zc[3]);

      if(zc[3])
      {
        if(si->dx_dX[2][1] == '.')
          errorexit("implement Y_of_x_forgiven_XZ(box, Y, x, *X,*Z);");
        else if(si->dx_dX[2][2] == '.')
          Y_of_y_forgiven_XZ(box, Y, y, *X, *Z);
        else if(si->dx_dX[2][3] == '.')
          errorexit("implement Y_of_z_forgiven_XZ(box, Y, z, *X,*Z);");
      }
    }
    else if(dir[2])
    {
      if(zc[3])
      {
        if(si->dx_dX[1][1] == '.')
          stat = X_of_x_forgiven_YZ(box, X, x, *Y,*Z);
        else if(si->dx_dX[1][2] == '.')
          errorexit("implement X_of_y_forgiven_YZ(box, X, y, *Y,*Z);");
        else if(si->dx_dX[1][3] == '.')
          errorexit("implement X_of_z_forgiven_YZ(box, X, z, *Y,*Z);");
      }
    }
    else if(dir[3])
    {
      if(zc[2])
      {
        if(si->dx_dX[1][1] == '.')
          stat = X_of_x_forgiven_YZ(box, X, x, *Y,*Z);
        else if(si->dx_dX[1][2] == '.')
          errorexit("implement X_of_y_forgiven_YZ(box, X, y, *Y,*Z);");
        else if(si->dx_dX[1][3] == '.')
          errorexit("implement X_of_z_forgiven_YZ(box, X, z, *Y,*Z);");
      }
    }
  }
  return stat;
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
  return -bi-10000;  /* make sure it's never 0 */
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

//if(!finite(X1) || !finite(Y1) || !finite(Z1))
//errorexit("not finite");
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
  return -bi-10000;  /* make sure it's never 0 */
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

/* find nearest X,Y,Z on grid from x,y,z */
double nearest_b_XYZ_of_xyz_inboxlist(tGrid *grid, int *blist, int nb, 
                            int *b, int *ind,
                            double *X, double *Y, double *Z,
                            double x, double y, double z)
{
  double r, rmin=-1.0;
  double X1,Y1,Z1;
  int ind1, bi, i;
  tBox *box;

  for(i=0; i<nb; i++)
  {
    bi = blist[i];
    box = grid->box[bi];
    r = nearestXYZ_of_xyz(box, &ind1, &X1,&Y1,&Z1, x,y,z);
    if(r<=rmin || rmin<0.0 )
    {
      rmin=r;  *b=bi;  *ind=ind1;   *X=X1;  *Y=Y1;  *Z=Z1;
    }
  }
  return rmin;
}

/* find nearest X,Y,Z inside (i.e. away from boundary) box from x,y,z */
double nearestinnerXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                              double x, double y, double z)
{
  int i,j,k, ijk;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  double *pX = box->v[Ind("X")];
  double *pY = box->v[Ind("Y")];
  double *pZ = box->v[Ind("Z")];
  double *px = box->v[Ind("x")];
  double *py = box->v[Ind("y")];
  double *pz = box->v[Ind("z")];
  double dxi, dyi, dzi, r, rmin=-1.0;

  for(k=1; k<n3-1; k++)
  for(j=1; j<n2-1; j++)
  for(i=1; i<n1-1; i++)
  {
    ijk = Index(i,j,k);
    dxi = px[ijk] - x;
    dyi = py[ijk] - y;
    dzi = pz[ijk] - z;
    r = dxi*dxi + dyi*dyi + dzi*dzi;
    if(r<=rmin || rmin<0.0 )
    {
      rmin = r;
      *ind = ijk;
      *X = pX[ijk];
      *Y = pY[ijk];
      *Z = pZ[ijk];
    }
  }
  return rmin;
}

/* find the faces a point X,Y,Z is on, face[2]=1 if X,Y,Z is on face2  */
int XYZ_on_face(tBox *box, int *face, double X, double Y, double Z)
{
  double XYZ[4];
  double *bb=box->bbox;
  int f;
  int nf;

  /* init */
  XYZ[1]=X;
  XYZ[2]=Y;
  XYZ[3]=Z;
  for(f=0; f<6; f++) face[f]=0;

  /* find all faces we are on */
  for(nf=0, f=0; f<6; f++)
  {
    int d=f/2;
    int dir=1+d;
    if(dequal(XYZ[dir],bb[f])) { face[f]=1; nf++; }
  }
  return nf; /* number of faces point is on */
}

/* find guess is nearest X,Y,Z in box from x,y,z */
double guessXYZ_of_xyz(tBox *box, int *ind, double *X, double *Y, double *Z,
                       double x, double y, double z)
{
  double rmin = nearestXYZ_of_xyz(box, ind, X,Y,Z, x,y,z);
  int face[6];

  /* if we are on a face, move away a bit */
  if(XYZ_on_face(box, face, *X,*Y,*Z))
  {
    int f, d;
    double *bb = box->bbox;
    double LX[4];
    double dX[4];
    for(d=0; d<3; d++)
    {
      LX[1+d] = bb[2*d+1]-bb[2*d];
      dX[1+d] = 0.0;
    }

    for(f=0; f<6; f++)
      if(face[f])
      {
        dX[1+f/2] = (1e-6 * LX[1+f/2]) * (1-2*(f%2));
      }
    *X += dX[1];
    *Y += dX[2];
    *Z += dX[3];
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

/* function to be passed into newton_linesrch_itsP by X_of_x_forgiven_YZ */
void x_VectorFuncP_YZ(int n, double *XYZvec, double *fvec, void *p)
{
  t_grid_box_desired_xyz_struct *pars;
  tBox *box;
  double desired_x;
  double xg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  /* get pars */
  pars = (t_grid_box_desired_xyz_struct *) p;
  box = pars->box;
  desired_x = pars->desired_x;
      
  xg = box->x_of_X[1]((void *) box, ind, XYZvec[1],XYZvec[2],XYZvec[3]);
  fvec[1] = xg-desired_x;
}

/* find X from x for a given Y,Z (Note: X also contains initial guess) */
int X_of_x_forgiven_YZ(tBox *box, double *X, double x, double Y, double Z)
{
  double XYZvec[4];
  t_grid_box_desired_xyz_struct pars[1];
  int check, stat;

  if(box->x_of_X[1]==NULL)  { *X = x;  return 0; }

  pars->box = box;
  pars->desired_x = x;

  /* initial guess supplied by caller */
  XYZvec[1] = *X;
  
  /* Y,Z are fixed */
  XYZvec[2] = Y;
  XYZvec[3] = Z;
      
  /* do newton_linesrch_itsP iterations: */
  stat = newton_linesrch_itsP(XYZvec, 1, &check, x_VectorFuncP_YZ,
                              (void *) pars,
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
  return -bi-10000;  /* make sure it's never 0 */
}

/* find X from x for a given Y,Z (Note: X also contains initial guess)
   search in all boxes in the list: bi = blist[0], blist[1], ..., blist[nb-1]
   return index of box in which X,Y,Z are found or error code if failure */
int b_X_of_x_forgiven_YZ_inboxlist(tGrid *grid, int *blist, int nb, 
                                   double *X, double x, double Y, double Z)
{
  double X1;
  int stat, i, bi;

  for(i=0; i<nb; i++)
  {
    bi = blist[i];
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
  return -bi-10000;  /* make sure it's never 0 */
}

/* function to be passed into newton_linesrch_itsP by Y_of_y_forgiven_XZ */
void y_VectorFuncP_XZ(int n, double *XYZvec, double *fvec, void *p)
{
  t_grid_box_desired_xyz_struct *pars;
  tBox *box;
  double desired_y;
  double yg;
  int ind=-1; /* works only if the x_of_X[i] don't use ind */

  /* get pars */
  pars = (t_grid_box_desired_xyz_struct *) p;
  box = pars->box;
  desired_y = pars->desired_y;

  yg = box->x_of_X[2]((void *) box, ind, XYZvec[2],XYZvec[1],XYZvec[3]);
  fvec[1] = yg-desired_y;
}

/* find Y from y for a given X,Z (Note: Y also contains initial guess) */
int Y_of_y_forgiven_XZ(tBox *box, double *Y, double y, double X, double Z)
{
  double XYZvec[4];
  t_grid_box_desired_xyz_struct pars[1];
  int check, stat;

  if(box->x_of_X[1]==NULL)  { *Y = y;  return 0; }

  pars->box = box;
  pars->desired_y = y;

  /* initial guess supplied by caller */
  XYZvec[1] = *Y;

  /* X,Z are fixed */
  XYZvec[2] = X;
  XYZvec[3] = Z;

  /* do newton_linesrch_itsP iterations: */
  stat = newton_linesrch_itsP(XYZvec, 1, &check, y_VectorFuncP_XZ,
                              (void *) pars,
                              Geti("Coordinates_newtMAXITS"),
                              Getd("Coordinates_newtTOLF") );
  *Y = XYZvec[1];

  if(check || stat<0) printf("Y_of_y_forgiven_XZ: check=%d stat=%d\n", check, stat);
  return stat-check;
}
