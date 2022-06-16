/* BCs_from_bfaces.c */
/* Wolfgang Tichy 2/2018 */

#include "sgrid.h"
#include "boundary.h"



/* find normal vector (n[1],n[2],n[3]) of box face f at point ijk */
void boxface_normal_at_ijk(tBox *box, int f, int ijk, double n[4])
{
  int dir = 1 + f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  int iX   = Ind("X");
  int idXd = Ind("dXdx");
  int idYd = Ind("dYdx");
  int idZd = Ind("dZdx");
  double dXdx[4][4];
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(box->x_of_X[1]!=NULL) /* not Cartesian */
  {
    /* check if the var dXdx is available */
    if(box->v[idXd]!=NULL)
    {
      /* get dXdx */
      for(j=1; j<=3; j++)
      {
        dXdx[1][j] = box->v[idXd + j-1][ijk];
        dXdx[2][j] = box->v[idYd + j-1][ijk];
        dXdx[3][j] = box->v[idZd + j-1][ijk];
      }
      /* get normal from derivs */
      n[1] = dXdx[dir][1];
      n[2] = dXdx[dir][2];
      n[3] = dXdx[dir][3];
    }
    else /* call dX_dx functions */
    {
      double X = box->v[iX][ijk];
      double Y = box->v[iX+1][ijk];
      double Z = box->v[iX+2][ijk];
      for(j=1; j<=3; j++)
        n[j] = box->dX_dx[dir][j](box, -1, X,Y,Z);
    }
  }
  else  /* use Cartesian normal */
  {
    n[1] = n[2] = n[3] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if(smag == 0.0) smag = sig;
  n[1] /= smag;
  n[2] /= smag;
  n[3] /= smag;
/*
printf("x=(%g,%g,%g)\n",
box->v[Ind("x")][ijk], box->v[Ind("y")][ijk], box->v[Ind("z")][ijk]);
printf("X=(%g,%g,%g)\n",
box->v[Ind("X")][ijk], box->v[Ind("Y")][ijk], box->v[Ind("Z")][ijk]);
printf("f%d xc=%g n=(%g,%g,%g)\n", bface->f, box->CI->xc[1], *nx,*ny,*nz);
int i;
for(i=1; i<=3; i++)
{
for(j=1; j<=3; j++) printf("%g   ", dXdx[i][j]);
printf("\n");
}
*/
}

/* find vector (n[1],n[2],n[3]) along coord. direction out of box face f
   at point ijk */
void boxface_outwarddir_at_ijk(tBox *box, int f, int ijk, double n[4])
{
  int dir = 1 + f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  int iX   = Ind("X");
  /*  int idxd = Ind("dxdX"); // we do not have dxdX yet!!!
      int idyd = Ind("dydX");
      int idzd = Ind("dzdX"); */
  int idxd=-1, idyd=-1, idzd=-1;
  double dxdX[4][4];
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(box->x_of_X[1]!=NULL) /* not Cartesian */
  {
    /* check if the var dxdX is available */
    if(idxd>0 && box->v[idxd]!=NULL)
    {
      /* get dxdX */
      for(j=1; j<=3; j++)
      {
        dxdX[1][j] = box->v[idxd + j-1][ijk];
        dxdX[2][j] = box->v[idyd + j-1][ijk];
        dxdX[3][j] = box->v[idzd + j-1][ijk];
      }
      /* get direction vectors from derivs */
      n[1] = dxdX[1][dir];
      n[2] = dxdX[2][dir];
      n[3] = dxdX[3][dir];
    }
    else if(box->dx_dX[1][dir]!=NULL) /* call dx_dX functions */
    {
      double X = box->v[iX][ijk];
      double Y = box->v[iX+1][ijk];
      double Z = box->v[iX+2][ijk];
      for(j=1; j<=3; j++)
        n[j] = box->dx_dX[j][dir](box, -1, X,Y,Z);
    }
    else /* try normal instead */
    {
      boxface_normal_at_ijk(box, f, ijk, n);
      return;
      errorexit("box->dx_dX[1][dir] is NULL");
    }
  }
  else  /* use Cartesian normal */
  {
    n[1] = n[2] = n[3] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if(smag == 0.0) smag = sig;
  n[1] /= smag;
  n[2] /= smag;
  n[3] /= smag;
}

/* find normal vector (n[1],n[2],n[3]) of box face f at X,Y,Z */
void boxface_normal_at_XYZ(tBox *box, int f, double X,double Y,double Z,
                           double n[4])
{
  int dir = 1 + f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(box->x_of_X[1]!=NULL) /* not Cartesian */
  {
    for(j=1; j<=3; j++)
      n[j] = box->dX_dx[dir][j](box, -1, X,Y,Z);
  }
  else  /* use Cartesian normal */
  {
    n[1] = n[2] = n[3] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if(smag == 0.0) smag = sig;
  n[1] /= smag;
  n[2] /= smag;
  n[3] /= smag;
}

/* find vector (n[1],n[2],n[3]) along coord. direction out of box face f
   at point X,Y,Z */
void boxface_outwarddir_at_XYZ(tBox *box, int f, double X,double Y,double Z,
                               double n[4])
{
  int dir = 1 + f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(box->x_of_X[1]!=NULL) /* not Cartesian */
  {
    /* check if the func dx_dX is available */
    if(box->dx_dX[1][dir]!=NULL) /* call dx_dX functions */
    {
      for(j=1; j<=3; j++)
        n[j] = box->dx_dX[j][dir](box, -1, X,Y,Z);
    }
    else /* try normal instead */
    {
      boxface_normal_at_XYZ(box, f, X,Y,Z, n);
      return;
      errorexit("box->dx_dX[1][dir] is NULL");
    }
  }
  else  /* use Cartesian normal */
  {
    n[1] = n[2] = n[3] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if(smag == 0.0) smag = sig;
  n[1] /= smag;
  n[2] /= smag;
  n[3] /= smag;
}


/* interpolate onto points in box on bface 
   this uses 1d interplolation done in other box, in direction idir
   if bface->setnormalderiv == 0 we set:
     FPsi = Psi(box) - Psi_interp(obox)
   else
     FPsi[ind] = n[i] [ dPsi[i](box) - dPsi[i]_interp(obox) ]  */
void FPsi_1Dinterp_for_bface(int iFPsi, tBface *bface, int idir,
                             int iPsi, int idPsi[4])
{
  tGrid *grid = bface->grid;
  tBox *box = grid->box[bface->b];
  int n1 = box->n1;
  int n2 = box->n2;
  /* int n3 = box->n3;
     int dir = 1 + bface->f/2; */
  int ob  = bface->ob;
  int ofi = bface->ofi;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  tBox *obox = NULL;
  int oCi;
  tBface *obface;
  int of, odir, op, plN;
  int pi, ind;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *dPsi[4];
  double *P, *Pcoeffs[4], *dP[4];
  dPsi[1] = box->v[idPsi[1]];
  dPsi[2] = box->v[idPsi[2]];
  dPsi[3] = box->v[idPsi[3]];

  /* check if there is another box and another bface */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(ofi<0) errorexit("can only interpolate if there is another bface");

  /* set other box props */
  obox = grid->box[ob];
  obface = obox->bface[ofi];
  of = obface->f;
  odir = 1+of/2;

  P = obox->v[iPsi];                  /* values in box ob */
  dP[1] = obox->v[idPsi[1]];          /* values in box ob */
  dP[2] = obox->v[idPsi[2]];          /* values in box ob */
  dP[3] = obox->v[idPsi[3]];          /* values in box ob */
  Pcoeffs[1] = dmalloc(obox->nnodes); /* get mem for coeffs of P or dP */
  Pcoeffs[2] = dmalloc(obox->nnodes); /* get mem for coeffs */
  Pcoeffs[3] = dmalloc(obox->nnodes); /* get mem for coeffs */

  /* -1 means op and plN are not set yet */
  op = plN = -1;

  /* if interp-dir is in other face plane */
  if(odir!=idir)
  {
    plN = odir;
    if(odir==1)      op = ( (obox->n1-1) )*(of%2);
    else if(odir==2) op = ( (obox->n2-1) )*(of%2);
    else             op = ( (obox->n3-1) )*(of%2);
  }
  else errorexit("implement odir==idir case");

  /* set oCi to oXi, oYi, or oZi */
  if(idir==1)      oCi = oXi;
  else if(idir==2) oCi = oYi;
  else if(idir==3) oCi = oZi;
  else errorexit("1<=idir<=3 is required");
  if(oCi<=0) errorexit("oCi <= 0");

  /* get coeffs in dir idir and plane plN,op */
  if(bface->setnormalderiv == 0)
    spec_analysis1_inplaneN(obox, idir, plN, op, P, Pcoeffs[1]);
  else
  {
    spec_analysis1_inplaneN(obox, idir, plN, op, dP[1], Pcoeffs[1]);
    spec_analysis1_inplaneN(obox, idir, plN, op, dP[2], Pcoeffs[2]);
    spec_analysis1_inplaneN(obox, idir, plN, op, dP[3], Pcoeffs[3]);
  }

  /* loop over points in bface */
  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    int i1, i2;
    int oi, oj, ok;
    double oX = box->v[oXi][ind];
    double oY = box->v[oYi][ind];
    double oZ = box->v[oZi][ind];
    double oC = box->v[oCi][ind]; /* this can be oX, oY, or oZ */
    double Pinterp[4], n[4];

    errorexit("are the oi, oj, ok below actually correct???");

    if(odir==1)
    {
      oi = op;
      oj = find_j_Of_Y(obox, oY);
      ok = find_k_Of_Z(obox, oZ);
    }
    else if(odir==2)
    {
      oj = op;
      oi = find_i_Of_X(obox, oX);
      ok = find_k_Of_Z(obox, oZ);
    }
    else
    {
      ok = op;
      oi = find_i_Of_X(obox, oX);
      oj = find_j_Of_Y(obox, oY);
    }

    if(idir==1)      { i1 = oj; i2 = ok; }
    else if(idir==2) { i1 = oi; i2 = ok; }
    else             { i1 = oi; i2 = oj; }

    /* Do not set BC if i1 or i2 are negative. This can happen e.g.
       if this is not the last BC we set at this point, because oX,oY,oZ
       contains only the coords of the point of the other box we use
       to set the last BC at this point. */
    if(i1<0 || i2<0) continue;
    //errorexit("i1<0 || i2<0, one of these may be wrong: oi, oj, ok");

    if(bface->setnormalderiv == 0)
    {
      Pinterp[1] = spec_interpolate_in_dir_at_i1_i2(obox, idir, 
                                                    i1,i2, Pcoeffs[1], oC);
      FPsi[ind] = Psi[ind] - Pinterp[1];
    }
    else
    {
      Pinterp[1] = spec_interpolate_in_dir_at_i1_i2(obox, idir, 
                                                    i1,i2, Pcoeffs[1], oC);
      Pinterp[2] = spec_interpolate_in_dir_at_i1_i2(obox, idir, 
                                                    i1,i2, Pcoeffs[2], oC);
      Pinterp[3] = spec_interpolate_in_dir_at_i1_i2(obox, idir, 
                                                    i1,i2, Pcoeffs[3], oC);
      boxface_normal_at_ijk(box, bface->f, ind, n);
      FPsi[ind] = n[1] * (dPsi[1][ind] - Pinterp[1]) +
                  n[2] * (dPsi[2][ind] - Pinterp[2]) +
                  n[3] * (dPsi[3][ind] - Pinterp[3]);
    }
    if(!finit(Pinterp[1]))
    {
      printf("Pinterp[1]=%g  i1=%d oC=%.13g i2=%d  ind=%d\n",
             Pinterp[1], i1,oC,i2, ind);
      printbface(bface);
      printbface(obface);
      grid->time  = 42;
      write_grid(grid);
      errorexit("Pinterp[1] is not finite!");
    }
  }
  free(Pcoeffs[3]);
  free(Pcoeffs[2]);
  free(Pcoeffs[1]);
}

/* interpolate onto points in box on bface
   this uses 2d interplolation is done in other box, in plane plN
   if bface->setnormalderiv == 0 we set:
     FPsi = Psi(box) - Psi_interp(obox)
   else
     FPsi[ind] = n_i [ dPsi_i(box) - dPsi_i_interp(obox) ]  */
void FPsi_2Dinterp_for_bface(int iFPsi, tBface *bface, int plN,
                             int iPsi, int idPsi[4])
{
  tGrid *grid = bface->grid;
  tBox *box = grid->box[bface->b];
  /* int n1 = box->n1;
     int n2 = box->n2;
     int n3 = box->n3;
     int dir = 1 + bface->f/2; */
  int ob  = bface->ob;
  int ofi = bface->ofi;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  tBox *obox = NULL;
  int oCi1, oCi2;
  double rangeC1[2], rangeC2[2];
  tBface *obface;
  int of, op;
  int pi, ind;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *dPsi[4];
  double *P, *Pcoeffs[4], *dP[4];
  dPsi[1] = box->v[idPsi[1]];
  dPsi[2] = box->v[idPsi[2]];
  dPsi[3] = box->v[idPsi[3]];


  /* check if there is another box and another bface */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(ofi<0) errorexit("can only interpolate if there is another bface");

  /* set other box props */
  obox = grid->box[ob];
  obface = obox->bface[ofi];
  of = obface->f;
  //odir = 1+of/2;

  P = obox->v[iPsi];                  /* values in box ob */
  dP[1] = obox->v[idPsi[1]];          /* values in box ob */
  dP[2] = obox->v[idPsi[2]];          /* values in box ob */
  dP[3] = obox->v[idPsi[3]];          /* values in box ob */
  Pcoeffs[1] = dmalloc(obox->nnodes); /* get mem for coeffs of P or dP */
  Pcoeffs[2] = dmalloc(obox->nnodes); /* get mem for coeffs */
  Pcoeffs[3] = dmalloc(obox->nnodes); /* get mem for coeffs */

  /* find plane index and set oCi1/2 to oXi, oYi, or oZi */
  if(plN==1)
  {
    op = ( (obox->n1-1) )*(of%2);
    oCi1 = oYi;
    oCi2 = oZi;
    rangeC1[0] = obox->bbox[2];
    rangeC1[1] = obox->bbox[3];
    rangeC2[0] = obox->bbox[4];
    rangeC2[1] = obox->bbox[5];
  }
  else if(plN==2)
  {
    op = ( (obox->n2-1) )*(of%2);
    oCi1 = oXi;
    oCi2 = oZi;
    rangeC1[0] = obox->bbox[0];
    rangeC1[1] = obox->bbox[1];
    rangeC2[0] = obox->bbox[4];
    rangeC2[1] = obox->bbox[5];
  }
  else if(plN==3)
  {
    op = ( (obox->n3-1) )*(of%2);
    oCi1 = oXi;
    oCi2 = oYi;
    rangeC1[0] = obox->bbox[0];
    rangeC1[1] = obox->bbox[1];
    rangeC2[0] = obox->bbox[2];
    rangeC2[1] = obox->bbox[3];
  }
  else errorexit("1<=plN<=3 is required");
  if(oCi1<=0 || oCi2<=0) errorexit("oCi1/2 > 0 is required");

  /* get coeffs for 2d interpolation from box ob */
  if(bface->setnormalderiv == 0)
    spec_Coeffs_inplaneN(obox, plN,op, P, Pcoeffs[1]);
  else
  {
    spec_Coeffs_inplaneN(obox, plN,op, dP[1], Pcoeffs[1]);
    spec_Coeffs_inplaneN(obox, plN,op, dP[2], Pcoeffs[2]);
    spec_Coeffs_inplaneN(obox, plN,op, dP[3], Pcoeffs[3]);
  }

  /* loop over bface points */
  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    double C1 = box->v[oCi1][ind];
    double C2 = box->v[oCi2][ind];
    double Pinterp[4], n[4];

    /* skip point if C1 or C2 are out of their range */
    if(C1<rangeC1[0] || C1>rangeC1[1]) continue;
    if(C2<rangeC2[0] || C2>rangeC2[1]) continue;

    if(bface->setnormalderiv == 0)
    {
      Pinterp[1] = spec_interpolate_inplaneN(obox, plN, op,
                                             Pcoeffs[1], C1,C2);
      FPsi[ind] = Psi[ind] - Pinterp[1];
    }
    else
    {
      Pinterp[1] = spec_interpolate_inplaneN(obox, plN, op,
                                             Pcoeffs[1], C1,C2);
      Pinterp[2] = spec_interpolate_inplaneN(obox, plN, op,
                                             Pcoeffs[2], C1,C2);
      Pinterp[3] = spec_interpolate_inplaneN(obox, plN, op,
                                             Pcoeffs[3], C1,C2);
      boxface_normal_at_ijk(box, bface->f, ind, n);
      FPsi[ind] = n[1] * (dPsi[1][ind] - Pinterp[1]) +
                  n[2] * (dPsi[2][ind] - Pinterp[2]) +
                  n[3] * (dPsi[3][ind] - Pinterp[3]);
    }
    if(!isfinite(Pinterp[1]))
    {
      printf("Pinterp[1]=%g  C1=%.13g C2=%.13g  oCi1=%d oCi2=%d  ind=%d\n",
              Pinterp[1], C1,C2, oCi1,oCi2, ind);
      printf("obox->b=%d plN=%d op=%d Pcoeffs[1][0]=%g\n",
              obox->b, plN, op, Pcoeffs[1][0]);
      printf("obox->bbox = [%g,%g] [%g,%g] [%g,%g]\n", obox->bbox[0], obox->bbox[1],
              obox->bbox[2], obox->bbox[3], obox->bbox[4], obox->bbox[5]);
      printCI(obox);
      printf("box->b=%d\n", box->b);
      printf("box->bbox = [%g,%g] [%g,%g] [%g,%g]\n", box->bbox[0], box->bbox[1],
              box->bbox[2], box->bbox[3], box->bbox[4], box->bbox[5]);
      printCI(box);
      printf("outputting oX oY oZ with quick_Vars_output at time 42\n");
      quick_Vars_output(box, "oX oY oZ", 42.0, 42);
      printf("bface:\n");
      printbface(bface);
      printf("obface:\n");
      printbface(obface);
      grid->time  = 42.1;
      write_grid(grid);
      errorexit("Pinterp[1] is not finite!");
    }
  }
  free(Pcoeffs[3]);
  free(Pcoeffs[2]);
  free(Pcoeffs[1]);
}


/* interpolate onto points in box on bface fi 
   this uses 3d interplolation is done in other box
   if bface->setnormalderiv == 0 we set:
     FPsi = Psi(box) - Psi_interp(obox)
   else
     FPsi[ind] = n_i [ dPsi_i(box) - dPsi_i_interp(obox) ]  */
void FPsi_3Dinterp_for_bface(int iFPsi, tBface *bface, int iPsi, int idPsi[4])
{
  tGrid *grid = bface->grid;
  tBox *box = grid->box[bface->b];
  tBox *obox = NULL;
  int ob  = bface->ob;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  int pi, ind;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *dPsi[4];
  double *P, *Pcoeffs[4], *dP[4];
  dPsi[1] = box->v[idPsi[1]];
  dPsi[2] = box->v[idPsi[2]];
  dPsi[3] = box->v[idPsi[3]];

  /* check if there is another box */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(oXi<=0 || oYi <=0 || oZi<=0)
    errorexit("oXi, oYi, oZi are not all > 0");

  /* set other box props */
  obox = grid->box[ob];
  P = obox->v[iPsi];                  /* values in box ob */
  dP[1] = obox->v[idPsi[1]];          /* values in box ob */
  dP[2] = obox->v[idPsi[2]];          /* values in box ob */
  dP[3] = obox->v[idPsi[3]];          /* values in box ob */
  Pcoeffs[1] = dmalloc(obox->nnodes); /* get mem for coeffs of P or dP */
  Pcoeffs[2] = dmalloc(obox->nnodes); /* get mem for coeffs */
  Pcoeffs[3] = dmalloc(obox->nnodes); /* get mem for coeffs */

  /* get values by 3d interpolation from box ob */
  if(bface->setnormalderiv == 0)
    spec_Coeffs(obox, P, Pcoeffs[1]);
  else
  {
    spec_Coeffs(obox, dP[1], Pcoeffs[1]);
    spec_Coeffs(obox, dP[2], Pcoeffs[2]);
    spec_Coeffs(obox, dP[3], Pcoeffs[3]);
  }
  
  /* loop over bface points */
  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    double X = box->v[oXi][ind];
    double Y = box->v[oYi][ind];
    double Z = box->v[oZi][ind];
    double Pinterp[4], n[4];

    if(bface->setnormalderiv == 0)
    {
      Pinterp[1] = spec_interpolate(obox, Pcoeffs[1], X,Y,Z);
      FPsi[ind] = Psi[ind] - Pinterp[1];
    }
    else
    {
      Pinterp[1] = spec_interpolate(obox, Pcoeffs[1], X,Y,Z);
      Pinterp[2] = spec_interpolate(obox, Pcoeffs[2], X,Y,Z);
      Pinterp[3] = spec_interpolate(obox, Pcoeffs[3], X,Y,Z);
      boxface_normal_at_ijk(box, bface->f, ind, n);
      FPsi[ind] = n[1] * (dPsi[1][ind] - Pinterp[1]) +
                  n[2] * (dPsi[2][ind] - Pinterp[2]) +
                  n[3] * (dPsi[3][ind] - Pinterp[3]);
    }
    if(!finit(Pinterp[1]))
    {
      printf("Pinterp[1]=%g  X=%.13g Y=%.13g Z=%.13g  ind=%d\n",
             Pinterp[1], X,Y,Z, ind);
      printbox(obox);
      grid->time  = 42;
      write_grid(grid);
      errorexit("Pinterp[1] is not finite!");
    }
  }
  free(Pcoeffs[3]);
  free(Pcoeffs[2]);
  free(Pcoeffs[1]);
}


/* find index in other box in the very special case where we have:
   +same coords in both boxes
   +same number of points in all directions
   +all points in touching interface coincide   */
int ijk_in_other_box_if_touch_samePoints_sameCoords(tBface *bface, int pi)
{
  tGrid *grid = bface->grid;
  tBox *box = grid->box[bface->b];
  int ob  = bface->ob;
  int ofi = bface->ofi;
  tBox *obox = NULL;
  int n1 = box->n1;
  int n2 = box->n2;
  /* int n3 = box->n3; */
  int ijk, of, od, ok,oj,oi, oijk;

  /* get ijk from pi */
  if(pi>=bface->fpts->npoints[bface->b]) return -1;
  ijk = bface->fpts->point[bface->b][pi]; /* ijk in box for pi */

  if(ob<0) return -1;

  /* get pointer to other box */
  obox = box->grid->box[ob];

  of = (obox->bface[ofi])->f;
  od = 1+of/2;
  ok = kOfInd_n1n2(ijk, n1,n2);
  oj = jOfInd_n1n2_k(ijk, n1,n2, ok);
  oi = iOfInd_n1n2_jk(ijk, n1,n2, oj,ok);

  if(od==1)      oi = ( (obox->n1-1) )*(of%2);
  else if(od==2) oj = ( (obox->n2-1) )*(of%2);
  else if(od==3) ok = ( (obox->n3-1) )*(of%2);
  oijk = Ind_n1n2(oi,oj,ok, obox->n1,obox->n2);
  return oijk;
}

/* find index ijk in other box in case of same_fpts */
int ijk_in_other_box_if_same_fpts(tBface *bface, int pi)
{
  int ob  = bface->ob;
  int ofi = bface->ofi;
  tBox *obox = NULL;
  tPointList *ofpts;
  int oijk;

  if(ob<0 || ofi<0) return -1;

  /* get pointer to other box */
  obox = bface->grid->box[ob];

  /* other fpts */
  ofpts = (obox->bface[ofi])->fpts;
  if(ofpts==NULL) return -1;

  /* index */
  if(pi>=ofpts->npoints[ob]) return -1;
  oijk = ofpts->point[ob][pi];

  return oijk;
}

/* if bface->same_fpts we can copy values
   if bface->setnormalderiv == 0 we set:
     FPsi = Psi(box) - Psi(obox)
   else
     FPsi[ind] = n_i [ dPsi_i(box) - dPsi_i(obox) ]  */
void FPsi_copy_for_bface(int iFPsi, tBface *bface, int iPsi, int idPsi[4],
                         int (*get_oijk)(tBface *bface, int pi))
{
  tGrid *grid = bface->grid;
  tBox *box = grid->box[bface->b];
  tBox *obox = NULL;
  int ob  = bface->ob;
  int pi, ind;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *dPsi[4];
  double *P, *dP[4];
  dPsi[1] = box->v[idPsi[1]];
  dPsi[2] = box->v[idPsi[2]];
  dPsi[3] = box->v[idPsi[3]];

  /* check if there is another box */
  if(ob<0)  errorexit("can only copy if there is another box");

  /* set other box props */
  obox = grid->box[ob];
  P = obox->v[iPsi];          /* values in box ob */
  dP[1] = obox->v[idPsi[1]];  /* values in box ob */
  dP[2] = obox->v[idPsi[2]];  /* values in box ob */
  dP[3] = obox->v[idPsi[3]];  /* values in box ob */

  if(bface->setnormalderiv==0) /* fields are equal */
  {
    forPointList_inbox(bface->fpts, box, pi, ind)
    {
      int oind = get_oijk(bface, pi);
      FPsi[ind] = Psi[ind] - P[oind];
    }
  }
  else /* field derivs are equal */
  {
    forPointList_inbox(bface->fpts, box, pi, ind)
    {
      int oind = get_oijk(bface, pi);
      double n[4];
      boxface_normal_at_ijk(box, bface->f, ind, n);
      FPsi[ind] = n[1] * (dPsi[1][ind] - dP[1][oind]) +
                  n[2] * (dPsi[2][ind] - dP[2][oind]) +
                  n[3] * (dPsi[3][ind] - dP[3][oind]);
    }
  }
}

/* set interbox BCs for all bfaces in this box */
void set_interbox_BCs_for_bface(int iFPsi, tBface *bface, 
                                int iPsi, int idPsi[4])
{
  /* tGrid *grid = bface->grid;
     tBox *box = grid->box[bface->b]; */
  int ob = bface->ob;

  /* check if there is another box */
  if(ob>=0)
  {
    if(bface->touch) /* bface touches one other face */
    {
      if(bface->same_fpts) /* we can copy */
        FPsi_copy_for_bface(iFPsi, bface, iPsi, idPsi,
                            ijk_in_other_box_if_same_fpts);
      else if(bface->sameX && bface->sameY && bface->sameZ)
      {
        int trouble = 0;
        if(bface->fpts==NULL) trouble=1;
        //add more ifs here that check if we really can use
        // ijk_in_other_box_if_touch_samePoints_sameCoords, i.e. if we have:
        // +touching boxes where all points coincide and same num. of points
        if(trouble)
          errorexit("if sameX=sameY=sameZ=1 we should have same_fpts=1");
        FPsi_copy_for_bface(iFPsi, bface, iPsi, idPsi,
                            ijk_in_other_box_if_touch_samePoints_sameCoords);
      }
      else if(bface->sameX && bface->sameY)
      {
        int idir = 3;
        FPsi_1Dinterp_for_bface(iFPsi, bface, idir, iPsi, idPsi);
      }
      else if(bface->sameX && bface->sameZ)
      {
        int idir = 2;
        FPsi_1Dinterp_for_bface(iFPsi, bface, idir, iPsi, idPsi);
      }
      else if(bface->sameY && bface->sameZ)
      {
        int idir = 1;
        FPsi_1Dinterp_for_bface(iFPsi, bface, idir, iPsi, idPsi);
      }
      else if(bface->sameX || bface->sameY || bface->sameZ)
      {
        /* we need 2d interpolation because the touching surfaces,
           do not have any points in common */
        int plN;
        if(bface->sameX) plN=1;
        if(bface->sameY) plN=2;
        if(bface->sameZ) plN=3;
        FPsi_2Dinterp_for_bface(iFPsi, bface, plN, iPsi, idPsi);
      }
      else if(bface->fpts_off_face)
        errorexit("implement fpts_off_face=1 case!");
      else
      {
        int ofi = bface->ofi;
        printf("There may be some flags missing on this bface...\n");
        printf("bface:\n");
        printbface(bface);
        if(ob>=0 && ofi >=0)
        {
          printf("other bface:\n");
          printbface(bface->grid->box[ob]->bface[ofi]);
        }
        errorexit("we should not get here");
      }
    }
    else /* bface is not associated with just one other face */
    {
      /* get values by 3d interpolation from box ob */
      FPsi_3Dinterp_for_bface(iFPsi, bface, iPsi, idPsi);
      if(bface->setnormalderiv)
        errorexit("this non-touching bface, should have setnormalderiv=0");
    }
  }
  //else /* there is no other box */
  //{
  //  /* set far limit BC */
  //  if(bface->outerbound) do_something
  //}
}

/* set interbox BCs for all bfaces in this box */
void set_interbox_BCs_for_boxbfaces(tBox *box, int iFPsi,
                                    int iPsi, int idPsi[4])
{
  int fi;
  /* loop over bfaces */
  forallbfaces(box, fi)
  {
    tBface *bface = box->bface[fi];
    set_interbox_BCs_for_bface(iFPsi, bface, iPsi, idPsi);
  }
}

/* interbox BCs for one face of a box and select setnormalderiv flag */
/* this can be used to set normal derivs or field values the same,
   if setnormalderiv = -1 it uses bface->setnormalderiv */
void set_interbox_BC_onface(tBox *box, int iFPsi, int face, int iPsi,
                            int idPsi[4], int setnormalderiv)
{
  int fi;
  tBface bface_copy[1];

  forallbfaces(box, fi)
  {
    tBface *bface = box->bface[fi];

    /* do nothing for wrong face */
    if(bface->f != face) continue;

    /* shallow copy of bface */
    *bface_copy = *bface;

    /* adjust setnormalderiv flag  */
    if(setnormalderiv >= 0)
      bface_copy->setnormalderiv = setnormalderiv;

    /* call general BC for bface_copy */
    set_interbox_BCs_for_bface(iFPsi, bface_copy, iPsi, idPsi);
  }
}
