/* set BCs */
/* Wolfgang Tichy 2017 */


#include "sgrid.h"
#include "Poisson3.h"


/* functions in this file */
/* ... */





/* find normal vector of box faces for Cartesian coords. */
void CartBox_normal(tBface *bface, int ijk, double *nx, double *ny, double *nz)
{
  *nx = *ny = *nz = 0.0;

  switch(bface->f)
  {
    case 0:
      *nx = -1.0;  break;
    case 1:
      *nx =  1.0;  break;
    case 2:
      *ny = -1.0;  break;
    case 3:
      *ny =  1.0;  break;
    case 4:
      *nz = -1.0;  break;
    case 5:
      *nz =  1.0;  break;
    default: 
      errorexit("bface->f must be 0,1,2,3,4,5");
  }
}

/* find normal vector of box face given by bface at point ijk */
void bface_normal(tBface *bface, int ijk, double *nx, double *ny, double *nz)
{
  tBox *box = bface->grid->box[bface->b];
  int idXd = Ind("dXdx");
  int idYd = Ind("dYdx");
  int idZd = Ind("dZdx");
  double dXdx[4][4];
  double dxdX[4][4];
  double mag;
  int j;

  if(box->v[idXd]==NULL) 
  {
    CartBox_normal(bface, ijk, nx,ny,nz);
    return;
  }

  /* get dXdx */
  for(j=1; j<=3; j++)
  {
    dXdx[1][j] = box->v[idXd + j-1][ijk];
    dXdx[2][j] = box->v[idYd + j-1][ijk];
    dXdx[3][j] = box->v[idZd + j-1][ijk];
  }

  /* get dxdX */
  dxdX_from_dXdx(dxdX, dXdx);

  /* get normal from derivs */
  switch(bface->f)
  {
    case 0:
      *nx = -dxdX[1][1];
      *ny = -dxdX[2][1];
      *nz = -dxdX[3][1];
      break;
    case 1:
      *nx = dxdX[1][1];
      *ny = dxdX[2][1];
      *nz = dxdX[3][1];
      break;
    case 2:
      *nx = -dxdX[1][2];
      *ny = -dxdX[2][2];
      *nz = -dxdX[3][2];
      break;
    case 3:
      *nx = dxdX[1][2];
      *ny = dxdX[2][2];
      *nz = dxdX[3][2];
      break;
    case 4:
      *nx = -dxdX[1][3];
      *ny = -dxdX[2][3];
      *nz = -dxdX[3][3];
      break;
    case 5:
      *nx = dxdX[1][3];
      *ny = dxdX[2][3];
      *nz = dxdX[3][3];
      break;
    default: 
      errorexit("bface->f must be 0,1,2,3,4,5");
  }
  /* normalize */
  mag = sqrt(nx[0]*nx[0] + ny[0]*ny[0] + nz[0]*nz[0]);
  if(mag>0.0)
  { 
    *nx /= mag;
    *ny /= mag;
    *nz /= mag;
  }
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
for(i=1; i<=3; i++)
{
for(j=1; j<=3; j++) printf("   %g", dxdX[i][j]);
printf("\n");
}
*/
}


/* find index in other box in case same_fpts */
int ind_in_other_box_if_same_fpts(tBox *box, int fi, int pi)
{
  tBface *bface = box->bface[fi];
  int ob  = bface->ob;
  int ofi = bface->ofi;
  tBox *obox = NULL;
  tPointList *ofpts;
  int oind;

  if(ob<0) return -1;

  /* get pointer to other box */
  obox = box->grid->box[ob];

  /* other fpts and index */
  ofpts = (obox->bface[ofi])->fpts;
  oind = ofpts->point[ob][pi];

  return oind;
}

/* find index in other box in case sameX=sameY=sameZ=1 */
int ind_in_other_box_if_sameXYZ(tBox *box, int fi, int ind)
{
  tBface *bface = box->bface[fi];
  int ob  = bface->ob;
  int ofi = bface->ofi;
  tBox *obox = NULL;
  /* int oXi = bface->oXi;   int oYi = bface->oYi;  int oZi = bface->oZi; */
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int of, od, ok,oj,oi, oind;

  if(ob<0) return -1;

  /* get pointer to other box */
  obox = box->grid->box[ob];

  of = (obox->bface[ofi])->f;
  od = 1+of/2;
  ok = kOfInd_n1n2(ind, n1,n2);
  oj = jOfInd_n1n2_k(ind, n1,n2, ok);
  oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);

  if(od==1)      oi = ( (obox->n1-1) )*(of%2);
  else if(od==2) oj = ( (obox->n2-1) )*(of%2);
  else if(od==3) ok = ( (obox->n3-1) )*(of%2);
  oind = Ind_n1n2(oi,oj,ok, obox->n1,obox->n2);
  return oind;
}

/* set BC's between boxes and at outerbound */
void Poisson3_set_interbox_and_outerBCs(tBox *box, int iFPsi, int iPsi,
                                        int iPsix, int iPsiy, int iPsiz,
                                        double (*outerBC)(double X, double Y, double Z),
                                        int setOuterBCs)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *Psix = box->v[iPsix];
  double *Psiy = box->v[iPsiy];
  double *Psiz = box->v[iPsiz];
  int iX = Ind("X");
  double *X = box->v[iX];
  double *Y = box->v[iX+1];
  double *Z = box->v[iX+2];
  int ix = Ind("x");
  double *px = box->v[ix];
  double *pz = box->v[ix+1];
  double *py = box->v[ix+2];
  int fi;

  /* loop over bfaces */
  for(fi=0; fi<box->nbfaces; fi++)
  {
    tBface *bface = box->bface[fi];
    int ob  = bface->ob;
    int ofi = bface->ofi;
    tBox *obox = NULL;
    int oXi = bface->oXi;
    int oYi = bface->oYi;
    int oZi = bface->oZi;
    int pi, ind;

    /* make sure oXi,... are positive */
    if(oXi<0) oXi = Ind("X");
    if(oYi<0) oYi = Ind("Y");
    if(oZi<0) oZi = Ind("Z");

    /* check if there is another box */
    if(ob>=0)
    {
      double *P;
      double *dP[4];

      obox = grid->box[ob];

      if(bface->touch) /* bface touches one other face */
      {
        if(bface->same_fpts)
        {
          if(bface->setnormalderiv) /* field derivs are equal */
          {
            dP[1] = obox->v[iPsix]; /* derivs in other box */
            dP[2] = obox->v[iPsiy];
            dP[3] = obox->v[iPsiz];
            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              int oind = ind_in_other_box_if_same_fpts(box, fi, pi);
              double nx, ny, nz;
              bface_normal(bface, ind, &nx, &ny, &nz);
              FPsi[ind] = nx * (Psix[ind] - dP[1][oind]) +
                          ny * (Psiy[ind] - dP[2][oind]) +
                          nz * (Psiz[ind] - dP[3][oind]);
            }
          }
          else /* fields are equal */
          {
            P = obox->v[iPsi]; /* values in box ob */
            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              int oind = ind_in_other_box_if_same_fpts(box, fi, pi);
              FPsi[ind] = Psi[ind] - P[oind];
            }
          }
        }
        else if(bface->sameX && bface->sameY && bface->sameZ)
        {
          if(bface->setnormalderiv) /* field derivs are equal */
          {
            dP[1] = obox->v[iPsix]; /* derivs in other box */
            dP[2] = obox->v[iPsiy];
            dP[3] = obox->v[iPsiz];
            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              int oind = ind_in_other_box_if_sameXYZ(box, fi, ind);
              double nx, ny, nz;
              CartBox_normal(bface, ind, &nx, &ny, &nz);
              FPsi[ind] = nx * (Psix[ind] - dP[1][oind]) +
                          ny * (Psiy[ind] - dP[2][oind]) +
                          nz * (Psiz[ind] - dP[3][oind]);
            }
          }
          else /* fields are equal */
          {
            P = obox->v[iPsi]; /* values in box ob */
            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              int oind = ind_in_other_box_if_sameXYZ(box, fi, ind);
              FPsi[ind] = Psi[ind] - P[oind];
            }
          }
        }
        else if(bface->sameX && bface->sameZ)
        {
          int n1 = box->n1;
          int n2 = box->n2;
          int n3 = box->n3;
          int dir = 1 + bface->f/2;
          tBface *obface = obox->bface[ofi];
          int of = obface->f;
          int odir = 1+of/2;
          int op;

          if(odir==1)      op = ( (obox->n1-1) )*(of%2);
          else if(odir==2) op = ( (obox->n2-1) )*(of%2);
          else if(odir==3) op = ( (obox->n3-1) )*(of%2);

          if(dir==1) /* X-dir => interpolate in Y-dir */
          {
            if(bface->setnormalderiv) /* field derivs are equal */
            {
              /* get deriv values by 1d interpolation from box ob */
              double dPinterp[4];
              double *dcoeffs[4];
              dcoeffs[1] = obox->v[Ind("Poisson3_temp1")];
              dcoeffs[2] = obox->v[Ind("Poisson3_temp2")];
              dcoeffs[3] = obox->v[Ind("Poisson3_temp3")];
              dP[1] = obox->v[iPsix]; /* derivs in other box */
              dP[2] = obox->v[iPsiy];
              dP[3] = obox->v[iPsiz];
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[1], dcoeffs[1]);
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[2], dcoeffs[2]);
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[3], dcoeffs[3]);

              forPointList_inbox(bface->fpts, box, pi, ind)
              {
                int ok = kOfInd_n1n2(ind, n1,n2);
                int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
                int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
                //double X = box->v[oXi][ind];
                double Y = box->v[oYi][ind];
                //double Z = box->v[oZi][ind];
                double nx, ny, nz;

                if(odir==1)      oi = op;
                else if(odir==2) oj = op;
                else if(odir==3) ok = op;

                dPinterp[1] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[1], Y);
                dPinterp[2] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[2], Y);
                dPinterp[3] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[3], Y);

                CartBox_normal(bface, ind, &nx, &ny, &nz);
                FPsi[ind] = nx * (Psix[ind] - dPinterp[1]) +
                            ny * (Psiy[ind] - dPinterp[2]) +
                            nz * (Psiz[ind] - dPinterp[3]);
                if(!finite(dPinterp[1]) ||
                   !finite(dPinterp[2]) || !finite(dPinterp[3]))
                {
                  NumberChecker_CheckIfFinite(grid, "Poisson3_temp1");
                  NumberChecker_CheckIfFinite(grid, "Poisson3_temp2");
                  NumberChecker_CheckIfFinite(grid, "Poisson3_temp3");
                  printf("dPinterp[1]=%g  oi=%d Y=%.13g ok=%d  ind=%d\n",
                         dPinterp[1], oi,Y,ok, ind);
                  printf("dPinterp[2]=%g\n", dPinterp[2]);
                  printf("dPinterp[3]=%g\n", dPinterp[3]);
                  printbface(bface);
                  printbface(obface);
                  grid->time  = 42;
                  write_grid(grid);
                  errorexit("dPinterp[] is not finite!");
                }
              }
            }
            else /* fields are equal */
            {
              /* get values by 1d interpolation from box ob */
              double Pinterp;
              double *Pcoeffs = obox->v[Ind("Poisson3_temp1")];
              P = obox->v[iPsi]; /* values in box ob */
              spec_analysis1_inplaneN(obox, 2, odir, op, P, Pcoeffs);

              forPointList_inbox(bface->fpts, box, pi, ind)
              {
                int ok = kOfInd_n1n2(ind, n1,n2);
                int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
                int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
                //double X = box->v[oXi][ind];
                double Y = box->v[oYi][ind];
                //double Z = box->v[oZi][ind];

                if(odir==1)      oi = op;
                else if(odir==2) oj = op;
                else if(odir==3) ok = op;

                Pinterp =
                  spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok, Pcoeffs, Y);

                FPsi[ind] = Psi[ind] - Pinterp;
                if(!finite(Pinterp))
                {
                  NumberChecker_CheckIfFinite(grid, "Poisson3_temp1");
                  printf("Pinterp=%g  oi=%d Y=%.13g ok=%d  ind=%d\n",
                         Pinterp, oi,Y,ok, ind);
                  printbface(bface);
                  printbface(obface);
                  grid->time  = 42;
                  write_grid(grid);
                  errorexit("Pinterp is not finite!");
                }
              }
            }
          }
          else
            errorexit("implement more dir cases!");
        }
        else if(bface->sameX || bface->sameY || bface->sameZ)
        {
          /* we need 2d interpolation because the touching surfaces,
             do not have any points in common */
          //errorexit("fix bface->sameX || bface->sameY || bface->sameZ case!");
          int n1 = box->n1;
          int n2 = box->n2;
          int n3 = box->n3;
          int dir = 1 + bface->f/2;
          tBface *obface = obox->bface[ofi];
          int of = obface->f;
          int odir = 1+of/2;
          int op;

          if(odir==1)      op = ( (obox->n1-1) )*(of%2);
          else if(odir==2) op = ( (obox->n2-1) )*(of%2);
          else if(odir==3) op = ( (obox->n3-1) )*(of%2);

          if(bface->setnormalderiv) /* field derivs are equal */
          {
            /* get deriv values by 1d interpolation from box ob */
            double dPinterp[4];
            double *dcoeffs[4];
            dcoeffs[1] = obox->v[Ind("Poisson3_temp1")];
            dcoeffs[2] = obox->v[Ind("Poisson3_temp2")];
            dcoeffs[3] = obox->v[Ind("Poisson3_temp3")];
            dP[1] = obox->v[iPsix]; /* derivs in other box */
            dP[2] = obox->v[iPsiy];
            dP[3] = obox->v[iPsiz];
            spec_Coeffs_inplaneN(obox, odir,op, dP[1], dcoeffs[1]);
            spec_Coeffs_inplaneN(obox, odir,op, dP[2], dcoeffs[2]);
            spec_Coeffs_inplaneN(obox, odir,op, dP[3], dcoeffs[3]);

            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              double X = box->v[oXi][ind];
              double Y = box->v[oYi][ind];
              double Z = box->v[oZi][ind];
              double X1,X2, nx, ny, nz;

              if(odir==1)      { X1 = Y; X2 = Z; }
              else if(odir==2) { X1 = X; X2 = Z; }
              else if(odir==3) { X1 = X; X2 = Y; }

              dPinterp[1] = spec_interpolate_inplaneN(obox, odir, op,
                                                      dcoeffs[1], X1,X2); 
              dPinterp[2] = spec_interpolate_inplaneN(obox, odir, op,
                                                      dcoeffs[2], X1,X2); 
              dPinterp[3] = spec_interpolate_inplaneN(obox, odir, op,
                                                      dcoeffs[3], X1,X2); 

              CartBox_normal(bface, ind, &nx, &ny, &nz);
              FPsi[ind] = nx * (Psix[ind] - dPinterp[1]) +
                          ny * (Psiy[ind] - dPinterp[2]) +
                          nz * (Psiz[ind] - dPinterp[3]);
              if(!finite(dPinterp[1]) ||
                 !finite(dPinterp[2]) || !finite(dPinterp[3]))
              {
                NumberChecker_CheckIfFinite(grid, "Poisson3_temp1");
                NumberChecker_CheckIfFinite(grid, "Poisson3_temp2");
                NumberChecker_CheckIfFinite(grid, "Poisson3_temp3");
                printf("dPinterp[1]=%g  X1=%.13g X2=%.13g  ind=%d\n",
                       dPinterp[1], X1,X2, ind);
                printf("dPinterp[2]=%g\n", dPinterp[2]);
                printf("dPinterp[3]=%g\n", dPinterp[3]);
                printbface(bface);
                printbface(obface);
                grid->time  = 42;
                write_grid(grid);
                errorexit("dPinterp[] is not finite!");
              }
            }
          }
          else /* fields are equal */
          {
            /* get values by 1d interpolation from box ob */
            double Pinterp;
            double *Pcoeffs = obox->v[Ind("Poisson3_temp1")];
            P = obox->v[iPsi]; /* values in box ob */
            spec_Coeffs_inplaneN(obox, odir,op, P, Pcoeffs);

            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              double X = box->v[oXi][ind];
              double Y = box->v[oYi][ind];
              double Z = box->v[oZi][ind];
              double X1,X2;

              if(odir==1)      { X1 = Y; X2 = Z; }
              else if(odir==2) { X1 = X; X2 = Z; }
              else if(odir==3) { X1 = X; X2 = Y; }

              Pinterp = spec_interpolate_inplaneN(obox, odir, op,
                                                  Pcoeffs, X1,X2); 
              FPsi[ind] = Psi[ind] - Pinterp;
              if(!finite(Pinterp))
              {
                NumberChecker_CheckIfFinite(grid, "Poisson3_temp1");
                printf("Pinterp=%g  X1=%.13g X2=%.13g  ind=%d\n",
                       Pinterp, X1,X2, ind);
                printbface(bface);
                printbface(obface);
                grid->time  = 42;
                write_grid(grid);
                errorexit("Pinterp is not finite!");
              }
            }
          }
        }
        else
          errorexit("implement more bface->same? cases!");
      }
      else /* bface is not associated with just one other face */
      {
        /* get values by 3d interpolation from box ob */
        double *P = obox->v[iPsi]; /* values in other box ob */
        double *Pcoeffs = obox->v[Ind("Poisson3_temp1")];
        spec_Coeffs(obox, P, Pcoeffs);

        forPointList_inbox(bface->fpts, box, pi, ind)
        {
          double X = box->v[oXi][ind];
          double Y = box->v[oYi][ind];
          double Z = box->v[oZi][ind];
          double Pinterp = spec_interpolate(obox, Pcoeffs, X,Y,Z);
          if(!finite(Pinterp))
          {
            printf("Pinterp=%g  X=%.13g Y=%.13g Z=%.13g  ind=%d\n",
                   Pinterp, X,Y,Z, ind);
            NumberChecker_CheckIfFinite(grid, "Poisson3_temp1");
            printbox(obox);
            grid->time  = 42;
            write_grid(grid);
            errorexit("Pinterp is not finite!");
          }
          FPsi[ind] = Psi[ind] - Pinterp;
        }
      }
    }
    else /* there is no other box */
    {
      /* set far limit BC */
      if(bface->outerbound && setOuterBCs)
        forPointList_inbox(bface->fpts, box, pi, ind)
        {
          double x,y,z;
          /* get x,y,z of point ind */
          if(px!=NULL)
          {
            x = px[ind];  y = py[ind];  z = pz[ind];
          }
          else
          {
            x = X[ind];   y = Y[ind];   z = Z[ind];
          }
          FPsi[ind] = Psi[ind] - outerBC(x,y,z);
        }
    }
  }
}
/* set BC's between boxes and at outerbound */
void Poisson3_set_interbox_and_outerBCs__new(tBox *box, int iFPsi, int iPsi,
                                        int iPsix, int iPsiy, int iPsiz,
                                        double (*outerBC)(double X, double Y, double Z),
                                        int setOuterBCs)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  int idPsi[4];
  int fi;
  idPsi[1] = iPsix;
  idPsi[2] = iPsiy;
  idPsi[3] = iPsiz;

  /* set BCs for cases where there is another box */
  set_interbox_BCs_for_bfaces(box, iFPsi, iPsi, idPsi);

  /* loop over bfaces */
  forallbfaces(box, fi)
  {
    tBface *bface = box->bface[fi];
    int ob  = bface->ob;
    int pi, ind;
    /* check if there is no box */
    if(ob<0)
    {
      /* set far limit BC */
      if(bface->outerbound && setOuterBCs)
      {
        int iX = Ind("X");
        double *X = box->v[iX];
        double *Y = box->v[iX+1];
        double *Z = box->v[iX+2];
        int ix = Ind("x");
        double *px = box->v[ix];
        double *pz = box->v[ix+1];
        double *py = box->v[ix+2];

        forPointList_inbox(bface->fpts, box, pi, ind)
        {
          double x,y,z;
          /* get x,y,z of point ind */
          if(px!=NULL)
          {
            x = px[ind];  y = py[ind];  z = pz[ind];
          }
          else
          {
            x = X[ind];   y = Y[ind];   z = Z[ind];
          }
          FPsi[ind] = Psi[ind] - outerBC(x,y,z);
        }
      }
    }
  } /* end of of forallbfaces */
}

/* outer BC for Psi */
double Poisson3_Psi_outerBC(double x, double y, double z)
{
  double r = sqrt(x*x + y*y + z*z);
//return 0.0;
  return 1.0/(4*PI*r)*erf(r);
}
/* outer BC is equal to 1 */
double Poisson3_outerBC_1(double x, double y, double z)
{
  return 1.0;
}
/* outer BC is equal to 0 */
double Poisson3_outerBC_0(double x, double y, double z)
{
  return 0.0;
}

/* new main BC routine, replaces set_BNSdata_BCs__old */
void Poisson3_set_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluAll,
                      int nonlin)
{
  tGrid *grid = vlu->grid;
  int vind;
  int vindDerivs=1;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  for(vind=0; vind<vlu->n; vind++)
  {
    int b;
    int iFPsi = vlFu->index[vind];
    int iPsi  = vlu->index[vind];
    int iPsix = vluAll->index[vindDerivs];
    int iPsiy = vluAll->index[vindDerivs+1];
    int iPsiz = vluAll->index[vindDerivs+2];
    int ncomp = VarNComponents(iPsi);
    double (*outerBC)(double X, double Y, double Z);

    /* do nothing and goto end of loop if var with vind is not the one
       of the current block */
    if(blkinfo!=NULL) if(vlu->index[vind] != blkinfo->vari)
                        goto Increm_vindDerivs;
    /* set func pointer for outer BC */
    if(nonlin)
    {
      if(vind==0) outerBC = Poisson3_Psi_outerBC;
      else        outerBC = Poisson3_outerBC_1;
    }
    else
      outerBC = Poisson3_outerBC_0;

    /* box loop */
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[iFPsi];
      double *Psi  = box->v[iPsi];
      double *Psix = box->v[iPsix];
      double *Psiy = box->v[iPsiy];
      double *Psiz = box->v[iPsiz];

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* set some BCs for each box */
      Poisson3_set_interbox_and_outerBCs(box, iFPsi, iPsi,
                                         iPsix,iPsiy,iPsiz, outerBC,1);
      //Poisson3_set_interbox_and_outerBCs__new(box, iFPsi, iPsi,
      //                                   iPsix,iPsiy,iPsiz, outerBC,1);
    } /* end forallboxes */

    Increm_vindDerivs:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp + 1;
  } /* end loop over vars */
}
