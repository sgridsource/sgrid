/* minimize_PN-NR_diff.c */
/* Wolfgang Tichy 5/2010 */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"

#define STRLEN 4000


/* global vars for func_to_minimize_for_numrec */
tBox *box__func_to_minimize_for_numrec;
double **ReNRPsi4mode__func_to_minimize_for_numrec;
double **ImNRPsi4mode__func_to_minimize_for_numrec;
int p_format__func_to_minimize_for_numrec;



/* read modes modes of NR Psi4 */
void read_sYlmModes_of_NR_Psi4(char *prefix,
                               double **ReNRPsi4mode, double **ImNRPsi4mode,
                               int lmax, int s, double t1, double t2, double dt)
{
  int i,n, l,m;
  int n1 = (lmax + 1)*(lmax + 1);
  FILE *in;
  char name[STRLEN];
  char str[STRLEN];
  int ndata = (t2-t1)/dt + 1;
  double t, time, Re, Im;

  printf("read_sYlmModes_of_NR_Psi4:\n");

  /* set all modes to zero at all times */
  for(i=0; i<n1; i++)
  for(n=0; n<ndata; n++)
    ReNRPsi4mode[i][n] = ImNRPsi4mode[i][n] = 0.0;

  i=4;  /* (1+1)^2 */
  for(l=2; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* open file */
    sprintf(name, "%sl%dm%d_s%d.t", prefix, l,m,s);
    in = fopen(name, "r");
    if(in) /* if file exists */
    {
      /* put 1e300 into ReNRPsi4mode to check later if we got data */
      double baddata=1e300;
      for(n=0; n<ndata; n++)
        ReNRPsi4mode[i][n] = baddata;

      /* read from file */
      printf("  reading %s\n", name);
      while(fscanline(in, str)!=EOF)
      {
        if(str[0] == '#' || str[0] == '\n') continue;
        sscanf(str, "%le %le %le", &time, &Re, &Im);
        //printf("%g\t %g\t %g\n", time, Re, Im);
        for(n=0; n<ndata; n++)
        {
          t = t1 + n*dt;
          if(dequal( (time-t)/t, 0.0 ))
          {
            ReNRPsi4mode[i][n] = Re;
            ImNRPsi4mode[i][n] = Im;
            break;
          }
        }
      } /* end while */
      fclose(in);

      /* check if we got data */
      for(n=0; n<ndata; n++)
      {
        t = t1 + n*dt;
        if(ReNRPsi4mode[i][n]==baddata && ImNRPsi4mode[i][n]==0.0)
        {
          printf("    could not read data at t=%.12g"
                 " --> setting NRPsi4mode=0.\n", t);
          ReNRPsi4mode[i][n] = ImNRPsi4mode[i][n] = 0.0;
        }
      }
    }
    else /* if file does not exist */
    {
      printf("  skipping %s\n", name);

      /* set mode to zero */
      for(n=0; n<ndata; n++)
        ReNRPsi4mode[i][n] = ImNRPsi4mode[i][n] = 0.0;
    }
    i++;
  }
}


/* compute NRPsi4 from its modes at time index n,
   set NRPsi4 on all radial indecies i=imin...imax */
void set_NRPsi4_inbox_from_modes(tBox *box, int Re_NRPsi4ind, int Im_NRPsi4ind,
                                 int Re_sYlmind, int Im_sYlmind,
                                 double **ReNRPsi4mode, double **ImNRPsi4mode,
                                 int n, int imin, int imax)
{
  int i,j,k, ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *Re_NRPsi4p = box->v[Re_NRPsi4ind];
  double *Im_NRPsi4p = box->v[Im_NRPsi4ind];
  double *Re_sYlmp = box->v[Re_sYlmind];
  double *Im_sYlmp = box->v[Im_sYlmind];
  double Rmode,Imode, RsYlm,IsYlm;
  double R,I;

  if(imax>=n1) errorexit("set_NRPsi4_inbox_from_modes: imax>=n1");

  /* go over angles */
  for(k=0; k<n3; k++)
  for(j=0; j<n2; j++)
  {
    R=I=0.0;
    for(i=0; i<n1; i++) /* this sums over modes */
    {
      ijk=Index(i,j,k);
    
      /* sYlm and modes of NRPsi4 */
      RsYlm = Re_sYlmp[ijk];
      IsYlm = Im_sYlmp[ijk];
      Rmode = ReNRPsi4mode[i][n];  /* n = time index */
      Imode = ImNRPsi4mode[i][n];

      /* set Re and Im part of Psi4 at angles given by j,k */
      R += (Rmode * RsYlm - Imode * IsYlm);
      I += (Rmode * IsYlm + Imode * RsYlm);
    }
    /* now set NRPsi4 on all radial indices i=imin...imax */
    for(i=imin; i<=imax; i++)
    {
      ijk=Index(i,j,k);
      Re_NRPsi4p[ijk] = R;
      Im_NRPsi4p[ijk] = I;
    }
  }
}


/* compute difference between PN and NR at radial index i=0 */
double PNPsi4_NRPsi4_diff_inbox(tBox *box, int Re_PNPsi4ind, int Im_PNPsi4ind,
                                int Re_NRPsi4ind, int Im_NRPsi4ind)
{
  int i,j,k, ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *Re_PNPsi4p = box->v[Re_PNPsi4ind];
  double *Im_PNPsi4p = box->v[Im_PNPsi4ind];
  double *Re_NRPsi4p = box->v[Re_NRPsi4ind];
  double *Im_NRPsi4p = box->v[Im_NRPsi4ind];
  double *PN_NR_diff = box->v[Ind("PN_CircularOrbit_GWs_PN_NR_diff")];
  double R,I;

  /* set integrand at i=0 */
  i=0;
  for(k=0; k<n3; k++)
  for(j=0; j<n2; j++)
  {
    ijk=Index(i,j,k);
    
    /* Re and Im of diff */
    R = Re_PNPsi4p[ijk] - Re_NRPsi4p[ijk];
    I = Im_PNPsi4p[ijk] - Im_NRPsi4p[ijk];

    /* compute diff^2 integrand */
    PN_NR_diff[ijk] = R*R + I*I;
  }
  
  /* do surface integral over diff^2 integrand */
  i=0;
  spec_sphericalDF2dIntegral_at_radial_index_i(box, PN_NR_diff, PN_NR_diff, i);
  
  return PN_NR_diff[0];
}


/* given some par compute total difference of PN-NR */
double PNPsi4_NRPsi4_totaldiff(tBox *box,
                               double **ReNRPsi4mode, double **ImNRPsi4mode,
                               double t1, double t2, double dt)
{
  int n;
  int ndata = (t2-t1)/dt + 1;
  double *diff;
  double time, tdiff;
  int Re_NRPsi4ind = Ind("PN_CircularOrbit_GWs_Re_NRPsi4"); /* NR Psi4 */
  int Im_NRPsi4ind = Ind("PN_CircularOrbit_GWs_Im_NRPsi4");
  int Re_PNPsi4ind = Ind("PN_CircularOrbit_GWs_Re_Psi4"); /* PN Psi4 */
  int Im_PNPsi4ind = Ind("PN_CircularOrbit_GWs_Im_Psi4");
  int hpind = Ind("PN_CircularOrbit_GWs_hplus");
  int hxind = Ind("PN_CircularOrbit_GWs_hcross");
  int Re_sYlmind = Ind("PN_CircularOrbit_GWs_Re_sYlm");
  int Im_sYlmind = Ind("PN_CircularOrbit_GWs_Im_sYlm");
  double yvec[12];  /* state vector used in odeint */
  double m1, m2, D; /* masses and distance of observer of h+,hx */
  double chi1x, chi1y, chi1z;  /* x,y,z comp of chi1 = S1/m1^2 */
  double chi2x, chi2y, chi2z;  /* x,y,z comp of chi2 = S2/m2^2 */
  double ti, tf;       /* initial and final time for integrator */
  int ImHmodesign;

  /* flip sign of Im H mode in case of plus_cross format (as for NINJA) */
  if(Getv("PN_CircularOrbit_GWs_HmodeOutputFormat", "plus_cross"))
    ImHmodesign=-1;
  else
    ImHmodesign=+1;

  //printf("PNPsi4_NRPsi4_totaldiff: Computing PN-NR diff\n");

  /* initial values for PN */
  m1 = Getd("PN_CircularOrbit_GWs_m1");
  m2 = Getd("PN_CircularOrbit_GWs_m2");
  chi1x = Getd("PN_CircularOrbit_GWs_chi1x");
  chi1y = Getd("PN_CircularOrbit_GWs_chi1y");
  chi1z = Getd("PN_CircularOrbit_GWs_chi1z");
  chi2x = Getd("PN_CircularOrbit_GWs_chi2x");
  chi2y = Getd("PN_CircularOrbit_GWs_chi2y");
  chi2z = Getd("PN_CircularOrbit_GWs_chi2z");
  /* chi1 = sqrt(chi1x*chi1x + chi1y*chi1y + chi1z*chi1z);
     chi2 = sqrt(chi2x*chi2x + chi2y*chi2y + chi2z*chi2z); */
  D = Getd("PN_CircularOrbit_GWs_D");
  yvec[0] = 0.0;  // you don't need to initialize this field it will be used as the storage of the separation
  yvec[1] = Getd("PN_CircularOrbit_GWs_omega")/(m1+m2); // initial orbital frequency
  yvec[2] = chi1x*m1*m1;  // S1x not Scap1! S1 = hi1*m1^2*Scap1
  yvec[3] = chi1y*m1*m1;  // S1y
  yvec[4] = chi1z*m1*m1;  // S1z
  yvec[5] = chi2x*m2*m2;  // S2x
  yvec[6] = chi2y*m2*m2;  // S2y
  yvec[7] = chi2z*m2*m2;  // S2z 
  yvec[8] = Getd("PN_CircularOrbit_GWs_Lnx");  // Lnx
  yvec[9] = Getd("PN_CircularOrbit_GWs_Lny");  // Lny
  yvec[10]= Getd("PN_CircularOrbit_GWs_Lnz");  // Lnz
  yvec[11]= Getd("PN_CircularOrbit_GWs_Phi");  // Phi orbital phase         
  xodeint(m1, m2, t1, t1, yvec); /* do this to initialize yvec[0] */

  /* allocate array for diff at different times */
  diff = dvector(0,ndata-1);
  
  /* loop over time */
  for(n=0; n<ndata; n++)
  {
    time = n*dt + t1;

    /* compute PNPsi4 at i=0 */
    compute_FDpsi4_and_hplus_hcross_on_sphere(box, Re_PNPsi4ind, Im_PNPsi4ind,
                                              hpind, hxind, yvec, D,m1,m2,
                                              0, dt*0.001, 0,0, 1);
//    /* get modes of Psi4 */
//    compute_sYlmModes_of_H(box, Re_Psi4ind, Im_Psi4ind,
//                           Re_sYlmind, Im_sYlmind, lmax,
//                           Re_Psi4mind, Im_Psi4mind, +1); 
//    /* output Psi4 modes */
//    output_sYlmModes_of_PN_H(Gets("PN_CircularOrbit_GWs_Psi4file_prefix"), 
//                             time, Re_Psi4modep, Im_Psi4modep,
//                             lmax,s, +1);

//    /* get modes of H = h+ - i hx  <-- sign of Im H is neg */
//    compute_sYlmModes_of_H(box, hpind, hxind, Re_sYlmind, Im_sYlmind, lmax,
//                           Re_Hmind, Im_Hmind, -1); 
//    /* output H modes */
//    output_sYlmModes_of_PN_H(Gets("PN_CircularOrbit_GWs_hfile_prefix"),
//                             time, Re_Hmodep, Im_Hmodep,
//                             lmax,s, ImHmodesign);

    /* advance orbit to time+dt */
    ti=time;
    tf=time+dt;
    xodeint(m1, m2, ti, tf, yvec);

    /* set NRPsi4 at radial index 0 */
    set_NRPsi4_inbox_from_modes(box, Re_NRPsi4ind, Im_NRPsi4ind,
                                Re_sYlmind, Im_sYlmind,
                                ReNRPsi4mode, ImNRPsi4mode, n, 0,0);
/* get modes of Psi4 */
//set_NRPsi4_inbox_from_modes(box, Re_NRPsi4ind, Im_NRPsi4ind,
//                            Re_sYlmind, Im_sYlmind,
//                            ReNRPsi4mode, ImNRPsi4mode, n, 0,box->n1-1);
//compute_sYlmModes_of_H(box, Re_NRPsi4ind, Im_NRPsi4ind,
//                       Re_sYlmind, Im_sYlmind, 2, hpind, hxind, +1); 
//double *Re_p = box->v[hpind];
//double *Im_p = box->v[hxind];
//output_sYlmModes_of_PN_H("NRpsi4_",time, Re_p, Im_p, 2,-2, +1);

    /* compute diff between PN and NR */
    diff[n] = PNPsi4_NRPsi4_diff_inbox(box, Re_PNPsi4ind, Im_PNPsi4ind,
                                       Re_NRPsi4ind, Im_NRPsi4ind);
  }

  /* integrate diff[n] over time to get tdiff */
  tdiff=integrate_simpson_1Dgrid(diff, dt, 0, ndata-1);

  /* free array for diff at different times */
  free_dvector(diff, 0,ndata-1);
  return tdiff;
}


/* thetax and phix are angles in speherical coords, 
   but with respect to x-axis (not z-axis) */
/* get nx,ny,nz from thetax and phix*/
void nx_ny_nz_from_thetax_phix(double *nx, double *ny, double *nz,
                                double thetax, double phix)
{
  *nx = cos(thetax);
  *ny = sin(thetax)*cos(phix);
  *nz = sin(thetax)*sin(phix);
}
/* get thetax and phix from nx,ny,nz */
void thetax_phix_from_nx_ny_nz(double *thetax, double *phix,
                                  double nx, double ny, double nz)
{
  *thetax = acos(nx);
  *phix = Arg(ny, nz);
}

/* set pars p[1...n] from sgrid pars */
int set_p_array_from_sgridpars(double *p, int p_format)
{
  double ax,ay,az, a;
  int npar;

  switch(p_format)
  {
    case 5:
      p[1] = Getd("PN_CircularOrbit_GWs_Phi");
      p[2] = Getd("PN_CircularOrbit_GWs_omega");
      p[3] = Getd("PN_CircularOrbit_GWs_D");
      p[4] = Getd("PN_CircularOrbit_GWs_m1");
      p[5] = Getd("PN_CircularOrbit_GWs_m2");
      npar = 5;
      break;
    case 7:
      p[1] = Getd("PN_CircularOrbit_GWs_Phi");
      p[2] = Getd("PN_CircularOrbit_GWs_omega");
      p[3] = Getd("PN_CircularOrbit_GWs_D");
      p[4] = Getd("PN_CircularOrbit_GWs_m1");
      p[5] = Getd("PN_CircularOrbit_GWs_m2");
      p[6] = Getd("PN_CircularOrbit_GWs_chi1z");
      p[7] = Getd("PN_CircularOrbit_GWs_chi2z");
      npar = 7;
      break;
    case 109:
      p[1] = Getd("PN_CircularOrbit_GWs_Phi");
      p[2] = Getd("PN_CircularOrbit_GWs_omega");
      p[3] = Getd("PN_CircularOrbit_GWs_D");

      ax = Getd("PN_CircularOrbit_GWs_chi1x");
      ay = Getd("PN_CircularOrbit_GWs_chi1y");
      az = Getd("PN_CircularOrbit_GWs_chi1z");
      a = sqrt(ax*ax + ay*ay + az*az);
      if(a==0.0) errorexit("p_format=109 does not work with chi1=0");
      thetax_phix_from_nx_ny_nz(p+4, p+5, ax/a, ay/a, az/a);

      ax = Getd("PN_CircularOrbit_GWs_chi2x");
      ay = Getd("PN_CircularOrbit_GWs_chi2y");
      az = Getd("PN_CircularOrbit_GWs_chi2z");
      a = sqrt(ax*ax + ay*ay + az*az);
      if(a==0.0) errorexit("p_format=109 does not work with chi2=0");
      thetax_phix_from_nx_ny_nz(p+6, p+7, ax/a, ay/a, az/a);

      thetax_phix_from_nx_ny_nz( p+8, p+9, 
                                 Getd("PN_CircularOrbit_GWs_Lnx"),
                                 Getd("PN_CircularOrbit_GWs_Lny"),
                                 Getd("PN_CircularOrbit_GWs_Lnz") );
      npar = 9;
      break;
    default :
      p[1] = Getd("PN_CircularOrbit_GWs_Phi");
      p[2] = Getd("PN_CircularOrbit_GWs_omega");
      p[3] = Getd("PN_CircularOrbit_GWs_D");

      p[4] = Getd("PN_CircularOrbit_GWs_chi1x");
      p[5] = Getd("PN_CircularOrbit_GWs_chi1y");
      p[6] = Getd("PN_CircularOrbit_GWs_chi1z");
      p[7] = Getd("PN_CircularOrbit_GWs_chi2x");
      p[8] = Getd("PN_CircularOrbit_GWs_chi2y");
      p[9] = Getd("PN_CircularOrbit_GWs_chi2z");
      /* p[10] = thetax; 
         p[11] = phix;  */
      thetax_phix_from_nx_ny_nz( p+10, p+11, 
                                 Getd("PN_CircularOrbit_GWs_Lnx"),
                                 Getd("PN_CircularOrbit_GWs_Lny"),
                                 Getd("PN_CircularOrbit_GWs_Lnz") );
      p[12] = Getd("PN_CircularOrbit_GWs_m1");
      p[13] = Getd("PN_CircularOrbit_GWs_m2");
      npar = 13;
  }
  return npar;
}

/* function we minimize with numrec's powell */
double func_to_minimize_for_numrec(double *p)
{
  double t1, t2, dt, tdiff;
  int n, npar;
  double nx,ny,nz, ax,ay,az, a;

  /* set sgrid pars from p array */
  switch(p_format__func_to_minimize_for_numrec)
  {
    case 1:
      npar = 1;
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      npar = 1;
      break;
    case 2:
      
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      npar = 2;
      break;
    case 3:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      npar = 3;
      break;
    case 5:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      Setd("PN_CircularOrbit_GWs_m1", p[4]);
      Setd("PN_CircularOrbit_GWs_m2", p[5]);
      npar = 5;
      break;
    case 7:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      Setd("PN_CircularOrbit_GWs_m1", p[4]);
      Setd("PN_CircularOrbit_GWs_m2", p[5]);
      Setd("PN_CircularOrbit_GWs_chi1z", p[6]);
      Setd("PN_CircularOrbit_GWs_chi2z", p[7]);
      npar = 7;
      break;
    case 9:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      Setd("PN_CircularOrbit_GWs_chi1x", p[4]);
      Setd("PN_CircularOrbit_GWs_chi1y", p[5]);
      Setd("PN_CircularOrbit_GWs_chi1z", p[6]);
      Setd("PN_CircularOrbit_GWs_chi2x", p[7]);
      Setd("PN_CircularOrbit_GWs_chi2y", p[8]);
      Setd("PN_CircularOrbit_GWs_chi2z", p[9]);
      npar = 9;
      break;
    case 109:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);

      ax = Getd("PN_CircularOrbit_GWs_chi1x");
      ay = Getd("PN_CircularOrbit_GWs_chi1y");
      az = Getd("PN_CircularOrbit_GWs_chi1z");
      a = sqrt(ax*ax + ay*ay + az*az);
      /* p[4] = thetax;  p[5] = phix;  */
      nx_ny_nz_from_thetax_phix(&nx, &ny, &nz, p[4], p[5]);
      Setd("PN_CircularOrbit_GWs_chi1x", a*nx);
      Setd("PN_CircularOrbit_GWs_chi1y", a*ny);
      Setd("PN_CircularOrbit_GWs_chi1z", a*nz);

      ax = Getd("PN_CircularOrbit_GWs_chi2x");
      ay = Getd("PN_CircularOrbit_GWs_chi2y");
      az = Getd("PN_CircularOrbit_GWs_chi2z");
      a = sqrt(ax*ax + ay*ay + az*az);
      /* p[6] = thetax;  p[7] = phix;  */
      nx_ny_nz_from_thetax_phix(&nx, &ny, &nz, p[6], p[7]);
      Setd("PN_CircularOrbit_GWs_chi2x", a*nx);
      Setd("PN_CircularOrbit_GWs_chi2y", a*ny);
      Setd("PN_CircularOrbit_GWs_chi2z", a*nz);

      /* p[8] = thetax;  p[9] = phix;  */
      nx_ny_nz_from_thetax_phix(&nx, &ny, &nz, p[8], p[9]);
      Setd("PN_CircularOrbit_GWs_Lnx", nx);
      Setd("PN_CircularOrbit_GWs_Lny", ny);
      Setd("PN_CircularOrbit_GWs_Lnz", nz);
      npar = 9;
      break;
    case 11:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      Setd("PN_CircularOrbit_GWs_chi1x", p[4]);
      Setd("PN_CircularOrbit_GWs_chi1y", p[5]);
      Setd("PN_CircularOrbit_GWs_chi1z", p[6]);
      Setd("PN_CircularOrbit_GWs_chi2x", p[7]);
      Setd("PN_CircularOrbit_GWs_chi2y", p[8]);
      Setd("PN_CircularOrbit_GWs_chi2z", p[9]);
      /* p[10] = thetax; 
         p[11] = phix;  */
      nx_ny_nz_from_thetax_phix(&nx, &ny, &nz, p[10], p[11]);
      Setd("PN_CircularOrbit_GWs_Lnx", nx);
      Setd("PN_CircularOrbit_GWs_Lny", ny);
      Setd("PN_CircularOrbit_GWs_Lnz", nz);
      npar = 11;
      break;
    case 13:
      Setd("PN_CircularOrbit_GWs_Phi", p[1]);
      Setd("PN_CircularOrbit_GWs_omega", p[2]);
      Setd("PN_CircularOrbit_GWs_D", p[3]);
      Setd("PN_CircularOrbit_GWs_chi1x", p[4]);
      Setd("PN_CircularOrbit_GWs_chi1y", p[5]);
      Setd("PN_CircularOrbit_GWs_chi1z", p[6]);
      Setd("PN_CircularOrbit_GWs_chi2x", p[7]);
      Setd("PN_CircularOrbit_GWs_chi2y", p[8]);
      Setd("PN_CircularOrbit_GWs_chi2z", p[9]);
      /* p[10] = thetax; 
         p[11] = phix;  */
      nx_ny_nz_from_thetax_phix(&nx, &ny, &nz, p[10], p[11]);
      Setd("PN_CircularOrbit_GWs_Lnx", nx);
      Setd("PN_CircularOrbit_GWs_Lny", ny);
      Setd("PN_CircularOrbit_GWs_Lnz", nz);

      Setd("PN_CircularOrbit_GWs_m1", p[12]);
      Setd("PN_CircularOrbit_GWs_m2", p[13]);
      npar = 13;
      break;
    default :
      errorexit("p_format__func_to_minimize_for_numrec has undefined value");
  }
  /* store npar in p[0] */
  p[0] = npar;

  /* print pars */
  for(n=1; n<=npar; n++)  printf(" %.4g", p[n]);

  /* set args for PNPsi4_NRPsi4_totaldiff */
  t1 = Getd("PN_CircularOrbit_GWs_t1");
  t2 = Getd("PN_CircularOrbit_GWs_t2");
  dt = Getd("PN_CircularOrbit_GWs_dt");

  tdiff = PNPsi4_NRPsi4_totaldiff(box__func_to_minimize_for_numrec,
                                  ReNRPsi4mode__func_to_minimize_for_numrec,
                                  ImNRPsi4mode__func_to_minimize_for_numrec,
                                  t1,t2,dt);
  printf(" => %.3g\n", tdiff);
  return tdiff;
}


/* minimize PN-NR differnce */
int minimize_PN_NR_diff(tGrid *grid)
{
  tBox *box = grid->box[0];
  int n1=box->n1;
  int ndata;
  double t1, t2, dt;
  int lmax, s;
  double **ReNRPsi4mode, **ImNRPsi4mode; /* arrays for NR Psi4 modes */
  int Re_sYlmind = Ind("PN_CircularOrbit_GWs_Re_sYlm");
  int Im_sYlmind = Ind("PN_CircularOrbit_GWs_Im_sYlm");
  double tdiff;
  int npar, nparmax = 20;
  double *p;   /* array for pars over which we minimize */
  double **xi; /* matrix with initial directions */
  int iter, i,j;
  double fret;
  double mintol = Getd("PN_CircularOrbit_GWs_mintol");
  int p_format = Getd("PN_CircularOrbit_GWs_min_p_format");

  /* allocate pars and matrix for powell */
  p  = dvector(0,nparmax);  /* Note: p[0] is not used by numrec. I use it to store number of pars */
  xi = dmatrix(1,nparmax, 1,nparmax);

  /* Info */
  printf("minimize_PN-NR_diff:\n");
  s    = Geti("PN_CircularOrbit_GWs_spinweight"); /* spin weight we use */
  lmax = Geti("PN_CircularOrbit_GWs_lmax");
  t1 = Getd("PN_CircularOrbit_GWs_t1");
  t2 = Getd("PN_CircularOrbit_GWs_t2");
  dt = Getd("PN_CircularOrbit_GWs_dt");
  ndata = (t2-t1)/dt + 1;
  printf("lmax=%d  s=%d  t1=%g  t2=%g  dt=%g\n", lmax, s, t1,t2,dt);

  /* allocate arrays for numerical Psi4 modes */
  ReNRPsi4mode = dmatrix(0,n1, 0,ndata);
  ImNRPsi4mode = dmatrix(0,n1, 0,ndata);
  
  /* read numerical Psi4 modes */
  read_sYlmModes_of_NR_Psi4(Gets("PN_CircularOrbit_GWs_NRPsi4file_prefix"),
                                 ReNRPsi4mode, ImNRPsi4mode, lmax,s, t1,t2,dt);
  /* precompute the sYlm */
  ModeComputer_set_sYlm_inbox(box, Re_sYlmind, Im_sYlmind, s, lmax);

  /* initial diff */
  tdiff = PNPsi4_NRPsi4_totaldiff(box, ReNRPsi4mode, ImNRPsi4mode, t1,t2,dt);
  printf("Difference before minimization: tdiff=%g\n", tdiff);

  /* call minimizer */
  /* set initial par guesses for minimizer */
  /* note p[0] is not set here */
  set_p_array_from_sgridpars(p, p_format);

  /* set initial directions for minimizer */
  for(i=1; i<=nparmax; i++)
  for(j=1; j<=nparmax; j++)
  {
    if(i==j) xi[i][j] = 1e-4;
    else     xi[i][j] = 0.0;
  }

  /* global vars for func_to_minimize_for_numrec */
  box__func_to_minimize_for_numrec = box;
  ReNRPsi4mode__func_to_minimize_for_numrec = ReNRPsi4mode;
  ImNRPsi4mode__func_to_minimize_for_numrec = ImNRPsi4mode;


  /* try with 1 par first */
  p_format__func_to_minimize_for_numrec = 1;
  set_p_array_from_sgridpars(p, p_format__func_to_minimize_for_numrec);
 
  /* set npar from p_format__func_to_minimize_for_numrec */
  func_to_minimize_for_numrec(p); /* get p[0] */
  npar = p[0];

  /* call minimizer */  
  powell(p, xi, npar, mintol, &iter, &fret, func_to_minimize_for_numrec);

  tdiff = PNPsi4_NRPsi4_totaldiff(box, ReNRPsi4mode, ImNRPsi4mode, t1,t2,dt);
  printf("Difference after first minimization: tdiff=%g\n", tdiff);

  if(p_format>=3)
  {
    /* try with 3 pars */
    p_format__func_to_minimize_for_numrec = 3;
    set_p_array_from_sgridpars(p, p_format__func_to_minimize_for_numrec);
   
    /* set npar from p_format__func_to_minimize_for_numrec */
    func_to_minimize_for_numrec(p); /* get p[0] */
    npar = p[0];

    /* call minimizer */  
    powell(p, xi, npar, mintol, &iter, &fret, func_to_minimize_for_numrec);

    tdiff = PNPsi4_NRPsi4_totaldiff(box, ReNRPsi4mode, ImNRPsi4mode, t1,t2,dt);
    printf("Difference after second minimization: tdiff=%g\n", tdiff);
  }

  if(p_format>3)
  {
    /* try with all pars */
    p_format__func_to_minimize_for_numrec = p_format;
    set_p_array_from_sgridpars(p, p_format__func_to_minimize_for_numrec);
   
    /* set npar from p_format__func_to_minimize_for_numrec */
    func_to_minimize_for_numrec(p); /* get p[0] */
    npar = p[0];
//Yo(0);
//func_to_minimize_for_numrec(p);
//exit(88);
    /* call minimizer */  
    powell(p, xi, npar, mintol, &iter, &fret, func_to_minimize_for_numrec);
  }

  /* print info at end */
  for(i=1; i<=npar; i++)  printf(" %.4g", p[i]);
  printf(":  iter=%d  fret=%g\n", iter, fret);

  /* final diff */
  tdiff = PNPsi4_NRPsi4_totaldiff(box, ReNRPsi4mode, ImNRPsi4mode, t1,t2,dt);
  printf("Difference after minimization: tdiff=%g\n", tdiff);

  /* make sure sgrid pars are equal to p array */
  func_to_minimize_for_numrec(p);

  /* tell min pars */
  printf("Parameters that minimize the difference are:\n");
  for(i=1; i<=npar; i++)  printf("p[%d] = %.20e\n", i, p[i]);

  printf("PN_CircularOrbit_GWs_Phi = %s\n", Gets("PN_CircularOrbit_GWs_Phi"));
  printf("PN_CircularOrbit_GWs_omega = %s\n", Gets("PN_CircularOrbit_GWs_omega"));
  printf("PN_CircularOrbit_GWs_D = %s\n", Gets("PN_CircularOrbit_GWs_D"));

  printf("PN_CircularOrbit_GWs_chi1x = %s\n", Gets("PN_CircularOrbit_GWs_chi1x"));
  printf("PN_CircularOrbit_GWs_chi1y = %s\n", Gets("PN_CircularOrbit_GWs_chi1y"));
  printf("PN_CircularOrbit_GWs_chi1z = %s\n", Gets("PN_CircularOrbit_GWs_chi1z"));
  printf("PN_CircularOrbit_GWs_chi2x = %s\n", Gets("PN_CircularOrbit_GWs_chi2x"));
  printf("PN_CircularOrbit_GWs_chi2y = %s\n", Gets("PN_CircularOrbit_GWs_chi2y"));
  printf("PN_CircularOrbit_GWs_chi2z = %s\n", Gets("PN_CircularOrbit_GWs_chi2z"));

  printf("PN_CircularOrbit_GWs_Lnx = %s\n", Gets("PN_CircularOrbit_GWs_Lnx"));
  printf("PN_CircularOrbit_GWs_Lny = %s\n", Gets("PN_CircularOrbit_GWs_Lny"));
  printf("PN_CircularOrbit_GWs_Lnz = %s\n", Gets("PN_CircularOrbit_GWs_Lnz"));

  printf("PN_CircularOrbit_GWs_m1 = %s\n", Gets("PN_CircularOrbit_GWs_m1"));
  printf("PN_CircularOrbit_GWs_m2 = %s\n", Gets("PN_CircularOrbit_GWs_m2"));

  printf("# time grid used:\n");
  printf("PN_CircularOrbit_GWs_t1 = %s\n", Gets("PN_CircularOrbit_GWs_t1"));
  printf("PN_CircularOrbit_GWs_t2 = %s\n", Gets("PN_CircularOrbit_GWs_t2"));
  printf("PN_CircularOrbit_GWs_dt = %s\n", Gets("PN_CircularOrbit_GWs_dt"));

  /* free arrays with numerical Psi4 modes */
  free_dmatrix(ReNRPsi4mode, 0,n1, 0,ndata+1);
  free_dmatrix(ImNRPsi4mode, 0,n1, 0,ndata+1);

  /* free pars and matrix for powell */
  free_dvector(p, 0,nparmax);
  free_dmatrix(xi, 1,nparmax, 1,nparmax);

  return 0;
}
