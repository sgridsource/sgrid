/* minimize_PN-NR_diff.c */
/* Wolfgang Tichy 5/2010 */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"

#define STRLEN 4000



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
  char str2[STRLEN];
  int ndata = (t2-t1)/dt + 1;
  double t, time, Re, Im;

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
      /* read from file */
      while(fscanline(in, str)!=EOF)
      {
        if(str[0] == '#' || str[0] == '\n') continue;
        sscanf(str, "%s", str2);  time=atof(str2);
        sscanf(str, "%s", str2);  Re=atof(str2);
        sscanf(str, "%s", str2);  Im=atof(str2);
        for(n=0; n<ndata; n++)
        {
          t = n*dt;
          if(dequal( (time-t)/t, 0.0 ))
          {
            ReNRPsi4mode[i][n] = Re;
            ImNRPsi4mode[i][n] = Im;
          }
        }
      } /* end while */
      fclose(in);
    }
    else /* if file does not exist */
    {
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
    /* now set NRPsi4 on all radial indecies i=imin...imax */
    for(i=imin; i<=imax; i++)
    {
      ijk=Index(i,j,k);
      Re_NRPsi4p[ijk] = R;
      Re_NRPsi4p[ijk] = I;
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

  printf("PNPsi4_NRPsi4_totaldiff: Computing PN-NR diff\n");

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
//    compute_sYlmModes_of_PN_H(box, Re_Psi4ind, Im_Psi4ind,
//                              Re_sYlmind, Im_sYlmind, lmax,
//                              Re_Psi4mind, Im_Psi4mind, +1); 
//    /* output Psi4 modes */
//    output_sYlmModes_of_PN_H(Gets("PN_CircularOrbit_GWs_Psi4file_prefix"), 
//                             time, Re_Psi4modep, Im_Psi4modep,
//                             lmax,s, +1);

//    /* get modes of H = h+ - i hx  <-- sign of Im H is neg */
//    compute_sYlmModes_of_PN_H(box, hpind, hxind, Re_sYlmind, Im_sYlmind, lmax,
//                              Re_Hmind, Im_Hmind, -1); 
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

  /* ??? pick PN pars */

  /* ??? test out diff */
tdiff = PNPsi4_NRPsi4_totaldiff(box, ReNRPsi4mode, ImNRPsi4mode, t1,t2,dt);
printf("tdiff=%g\n", tdiff);
                                                             
  /* free arrays with numerical Psi4 modes */
  free_dmatrix(ReNRPsi4mode, 0,n1, 0,ndata+1);
  free_dmatrix(ImNRPsi4mode, 0,n1, 0,ndata+1);

  return 0;
}
