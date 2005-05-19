/* BCs.c */
/* Wolfgang Tichy 5/2005 */

#include "sgrid.h"
#include "boundary.h"


/* set variable ui to value value on boundary */
void set_boundary_value(tPointList *PL, int ui, double value)
{
  tGrid *grid = PL->grid;
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *u = box->v[ui];

    forPointList_inbox(PL, box, pi , ijk)
      u[ijk] = value;
  }
}

/* set normal derivative of variable ui to value deriv on left boundary */
void set_boundary_normalderiv_leftBound(tPointList *PL, int direc, int ui, 
                                        double deriv)
{
  tGrid *grid = PL->grid;
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *u = box->v[ui];
    double *du= box->v[Ind("temp1")];
    double D00;

    spec_Deriv1(box, direc, u, du);
    if(direc==1)      D00=box->D1[0];
    else if(direc==2) D00=box->D2[0];
    else if(direc==3) D00=box->D3[0];
    else errorexit("set_boundary_normalderiv_leftBound: direc has to be 1,2,3");

    /* deriv = D00 u[0] + sum_{j=1...n-1} D0j u[j] */ 
    forPointList_inbox(PL, box, pi , ijk)
      u[ijk] = (deriv - du[ijk] + D00*u[ijk])/D00;
  }
}

/* set normal derivative of variable ui to value deriv on right boundary */
void set_boundary_normalderiv_rightBound(tPointList *PL, int direc, int ui, 
                                         double deriv)
{
  tGrid *grid = PL->grid;
  int bi, pi, ijk;
  
  forallboxes(grid,bi)
  {
    tBox *box=grid->box[bi];
    double *u = box->v[ui];
    double *du= box->v[Ind("temp1")];
    double Dmm;

    spec_Deriv1(box, direc, u, du);
    if(direc==1)      Dmm=box->D1[(box->n1)*(box->n1)-1];
    else if(direc==2) Dmm=box->D2[(box->n2)*(box->n2)-1];
    else if(direc==3) Dmm=box->D3[(box->n3)*(box->n3)-1];
    else errorexit("set_boundary_normalderiv_rightBound: direc has to be 1,2,3");

    /* deriv = Dmm u[m] + sum_{j=0...n-2} D0j u[j] */ 
    forPointList_inbox(PL, box, pi , ijk)
      u[ijk] = (deriv - du[ijk] + Dmm*u[ijk])/Dmm;
  }
}
