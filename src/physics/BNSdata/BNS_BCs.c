/* set BCs */
/* Wolfgang Tichy 2010 */


#include "sgrid.h"
#include "BNSdata.h"

#define Power pow


/* functions in this file */
/* ... */


/* standard BCs for all fields */
/* set BCs for a varlist */
void set_BNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  int vind;
  int vindDerivs=0;

  if( grid->box[0]->n1 != grid->box[1]->n1 ||
      grid->box[3]->n1 != grid->box[2]->n1 ||
      grid->box[1]->n1 != grid->box[2]->n1 ) 
    errorexit("all n1 in boxes0-3 must be the same because we currently use "
              "lines like:\n"
              "FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];\n"
              "where Psi = grid->box[0]->v[vi], P = grid->box[1]->v[vi]");
                      
  for(vind=0; vind<vlu->n; vind++)
  {
    int ncomp = VarNComponents(vlu->index[vind]);
    double PsiFarLimit = VarFarLimit(vlu->index[vind])*nonlin;
    int BCs_atInf = 1;
    int BCs_box1_2 = 1;
    int BCs_AxisAtOuterInterfaces = 1;
    //int BCs_box0_1 = 1;
    //int BCs_box3_2 = 1;
    char *varname = VarName(vlu->index[vind]);

    if(strstr(varname, "BNSdata_Sigma"))
    { 
      BCs_atInf = 0;
      BCs_box1_2 = 0;
      BCs_AxisAtOuterInterfaces = 0;
      /* printf("varname=%s\n", varname); */
    }

    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[vlFu->index[vind]];
      double *Psi  = box->v[vlu->index[vind]];
      double *Psix = box->v[vluDerivs->index[vindDerivs]];
      double *Psiy = box->v[vluDerivs->index[vindDerivs+1]];
      double *Psiz = box->v[vluDerivs->index[vindDerivs+2]];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;

      /* BCs */
      if(Getv("BNSdata_grid", "SphericalDF"))
      {
        forplane1(i,j,k, n1,n2,n3, 0)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*(vind+1)*nonlin;

        forplane1(i,j,k, n1,n2,n3, n1-1)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*(vind+1)*nonlin;
      }
      else if (Getv("BNSdata_grid", "AnsorgNS"))
      {
        double *P;
        double *dP[4];
        double *BM;

        /* special rho=0 case??? */
        if(b==0 || b==1 || b==2 || b==3)
        {
          int pl;
          char str[1000];
          snprintf(str, 999, "box%d_basis2", b);
          if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* loop over rho=0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                  FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                {
                  double Psi_rho = Psiy[Index(i,j,k)];
                  double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                  FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                }
              }
          }
          /* same as before, but also interpolate to rho=0 */
          else if(Getv("BNSdata_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
          /* same as before, but do it only in box0/3 at A,B=1,0 and A,B=1,1 */
          else if(Getv("BNSdata_regularization", 
                       "regularity_on_axis_at_center") && (b==0 || b==3) )
          {
            double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
            double *temp3 = box->v[Ind("BNSdata_temp3")];
            double *temp4 = box->v[Ind("BNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              i=n1-1;   /* do it only at A=1 */
              j=pl;
              for(k=0; k<n3; k++)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
        } /* end: special rho=0 case??? */

        if(b==0)  /* in box0 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box0 and box1 */
          P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box0 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==3)  /* in box3 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box3 and box2 */
          P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box3 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==1)
        {
          BM = (double *) calloc(max3(grid->box[0]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[0], 1, BM, 0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box0 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[0]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box1_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else if(b==2)
        {
          BM = (double *) calloc(max3(grid->box[3]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* values at A=1 are equal in box1 and box2 */
          P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[3], 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box3 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[3]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box2_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else errorexiti("b=%d should be impossible!", b);
      } /* end: else if (Getv("BNSdata_grid", "AnsorgNS")) */
      else if (Getv("BNSdata_grid", "4ABphi_2xyz"))
      {
        double *P;
        double *dP[4];
        double *X, *Y, *Z,  *xp, *yp, *zp;
        double *Pcoeffs;
        double Pinterp;
        double x,y,z;

        /* special rho=0 case??? */
        BNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                         Psi, Psix,Psiy,Psiz);
        /* cases for each box */
        if(b==0)  /* in box0 */
        {
          /* values at A=0 are equal in box0 and box1 */
          P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* values at A=Amin are interpolated from box5 */
          xp = box->v[Ind("x")];
          yp = box->v[Ind("y")];
          zp = box->v[Ind("z")];
          P = grid->box[5]->v[vlu->index[vind]]; /* values in box5 */
          Pcoeffs = grid->box[5]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[5], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
          {
            int ind=Index(i,j,k);
            x = xp[ind]; 
            y = yp[ind]; 
            z = zp[ind];
            Pinterp = spec_interpolate(grid->box[5], Pcoeffs, x,y,z);
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==3)  /* in box3 */
        {
          /* values at A=0 are equal in box3 and box2 */
          P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* values at A=Amin are interpolated from box4 */
          xp = box->v[Ind("x")];
          yp = box->v[Ind("y")];
          zp = box->v[Ind("z")];
          P = grid->box[4]->v[vlu->index[vind]]; /* values in box4 */
          Pcoeffs = grid->box[4]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[4], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amin */
          {
            int ind=Index(i,j,k);
            x = xp[ind]; 
            y = yp[ind]; 
            z = zp[ind]; 
            Pinterp = spec_interpolate(grid->box[4], Pcoeffs, x,y,z);
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==1)  /* in box1 */
        {
          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          if(BCs_box1_2)
          {
            dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
            forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
              FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          }
          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }
          
          /* Psi=0 at infinity */
          if(BCs_atInf)
          {
            if(Getv("box1_basis2", "ChebExtrema"))
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
                FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
            else // B=0 is not on grid for ChebZeros!!!
            {
              int l;
              double U0;
              double *BM = (double *) calloc(n2, sizeof(double));
              double *line = (double *) calloc(n2, sizeof(double));

              /* obtain BM vector for interpolation along B */
              spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              {
                /* find value of Psi at A=1, B=0 */
                get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
                FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
              }
              free(line);
              free(BM);
            }
          }
        }
        else if(b==2)  /* in box2 */
        {
          /* values at A=1 are equal in box1 and box2 */
          if(BCs_box1_2)
          {
            P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
            forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
              FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          }
          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }

          /* Psi=0 at infinity */
          if(BCs_atInf)
          {
            /* Psi=0 at infinity */
            if(Getv("box2_basis2", "ChebExtrema"))
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
                FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
            else // B=0 is not on grid for ChebZeros!!!
            {
              int l;
              double U0;
              double *BM = (double *) calloc(n2, sizeof(double));
              double *line = (double *) calloc(n2, sizeof(double));

              /* obtain BM vector for interpolation along B */
              spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
              {
                /* find value of Psi at A=1, B=0 */
                get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
                FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
              }
              free(line);
              free(BM);
            }
          }
        }
        else if(b==5)  /* in box5 */
        {
          /* values at border are interpolated from box0 */
          double A,B,phi;
          int pl; //, k_phi;
          double *pA = box->v[Ind("BNSdata_A")];
          double *pB = box->v[Ind("BNSdata_B")];
          double *pphi = box->v[Ind("BNSdata_phi")];
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[0]->v[vlu->index[vind]]; /* values in box0 */
          Pcoeffs = grid->box[0]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[0], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int ind=Index(pl,0,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int ind=Index(0,pl,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];
                             
              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int ind=Index(0,0,pl);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else if(b==4)  /* in box4 */
        {
          /* values at border are interpolated from box3 */
          double A,B,phi;
          int pl; //, k_phi;
          double *pA = box->v[Ind("BNSdata_A")];
          double *pB = box->v[Ind("BNSdata_B")];
          double *pphi = box->v[Ind("BNSdata_phi")];
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[3]->v[vlu->index[vind]]; /* values in box3 */
          Pcoeffs = grid->box[3]->v[Ind("BNSdata_temp1")];
          spec_Coeffs(grid->box[3], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int ind=Index(pl,0,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int ind=Index(0,pl,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int ind=Index(0,0,pl);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else errorexiti("b=%d should be impossible!", b);

        /* special rho=0 case again at A=0 and A=1 ???  */
        if( ( (b==0 || b==3) || 
              (b==1 || b==2) && BCs_AxisAtOuterInterfaces ) && 
            Getv("BNSdata_regularization",
                 "regularity_on_axis_at_interfaces") )
        {
          BNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                           Psi, Psix,Psiy,Psiz);
        } /* end: special rho=0 case again ??? */

      } /* end: else if (Getv("BNSdata_grid", "4ABphi_2xyz")) */

    } /* end forallboxes */
    /* increase index for derivs */
    vindDerivs += 3;
    if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}


/* treat rho=0 case: */
void BNSdata_RegularityConditions_for_Var_at_rho_eq_0(tBox *box, double *FPsi,
                        double *Psi, double *Psix, double *Psiy, double *Psiz)
{
  int b = box->b;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;

  /* special rho=0 case??? */
  if(b==0 || b==1 || b==2 || b==3)
  {
    int pl;
    char str[1000];
    snprintf(str, 999, "box%d_basis2", b);
    if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
    {
      double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
      double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
      double *temp3 = box->v[Ind("BNSdata_temp3")];
      double *temp4 = box->v[Ind("BNSdata_temp4")];

      /* get u_phi_phi */
      spec_Deriv2(box, 3, Psi, Psi_phi_phi);
      
      /* get u_rho_phi_phi at phi=0 */
      /* d/drho = dx^i/drho d/dx^i, 
         dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
         ==> d/drho u = d/dy u  at phi=0 */           
      /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
      cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

      /* loop over rho=0 boundary */
      for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
        forplane2(i,j,k, n1,n2,n3, pl)
        {
          if(k>0) /* phi>0: impose u_ijk = u_ij0 (not u_phi_phi=0) */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)]-Psi[Index(i,j,0)];
          else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
          {
            double Psi_rho = Psiy[Index(i,j,k)];
            double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
            FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
          }
        }
    }
    /* same as before, but also interpolate to rho=0 */
    else if(Getv("BNSdata_regularization", "regularity_on_axis"))
    {
      double *Psi_phi_phi = box->v[Ind("BNSdata_temp1")];
      double *Psi_y_phi_phi = box->v[Ind("BNSdata_temp2")];
      double *temp3 = box->v[Ind("BNSdata_temp3")];
      double *temp4 = box->v[Ind("BNSdata_temp4")];
      double *line = (double *) calloc(n2, sizeof(double));
      double *BM[2];
      BM[0] = (double *) calloc(n2, sizeof(double));
      BM[1] = (double *) calloc(n2, sizeof(double));

      /* get u_phi_phi */
      spec_Deriv2(box, 3, Psi, Psi_phi_phi);
      
      /* get u_rho_phi_phi at phi=0 */
      /* d/drho = dx^i/drho d/dx^i, 
         dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
         ==> d/drho u = d/dy u  at phi=0 */           
      /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
      cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

      /* obtain BM vectors for interpolation along B */
      spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
      spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

      /* loop over rho~0 boundary */
      for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
      {
        int l;
        double U0, V0;

        forplane2(i,j,k, n1,n2,n3, pl)
        {
          if(k>0) /* phi>0: impose u_phi_phi=0 */
          {
            /* find value Psi_phi_phi at B=0 or 1 */
            get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
            FPsi[Index(i,j,k)] = U0;
          }
          else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
          { /* Psi_rho = Psiy  
               Psi_rho_phi_phi = Psi_y_phi_phi */
            /* find value Psi_rho at B=0 or 1 */
            get_memline(Psiy, line, 2, i,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

            /* find value Psi_rho_phi_phi at B=0 or 1 */
            get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
            for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

            FPsi[Index(i,j,k)] = U0 + V0;
          }
        }
      }
      free(BM[0]);
      free(BM[1]);
      free(line);
    }
  } /* end: special rho=0 case */
}
