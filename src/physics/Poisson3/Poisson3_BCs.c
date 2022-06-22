/* Poisson3_BCs.c */
/* set BCs, Wolfgang Tichy 2017 */


#include "sgrid.h"
#include "Poisson3.h"


/* functions in this file */
/* ... */

/************************************************************************/
/* Functions that we really need in the new version of Poisson3_set_BCs */
/************************************************************************/

/* set BC's between boxes and at outerbound */
void Poisson3_set_interbox_and_outerBCs(tBox *box, int iFPsi, int iPsi,
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

  /* loop over bfaces */
  forallbfaces(box, fi)
  {
    tBface *bface = box->bface[fi];
    int ob  = bface->ob;
    int pi, ind;

    if(ob>=0)
    {
      /* set BCs for cases where there is another box */
      set_interbox_BCs_for_bface(iFPsi, bface, iPsi, idPsi);
    }
    else  /* there is no box */
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
      //Poisson3_set_interbox_and_outerBCs__old(box, iFPsi, iPsi,
      //                                   iPsix,iPsiy,iPsiz, outerBC,1);
      Poisson3_set_interbox_and_outerBCs(box, iFPsi, iPsi,
                                         iPsix,iPsiy,iPsiz, outerBC,1);
    } /* end forallboxes */

    Increm_vindDerivs:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp + 1;
  } /* end loop over vars */
}
