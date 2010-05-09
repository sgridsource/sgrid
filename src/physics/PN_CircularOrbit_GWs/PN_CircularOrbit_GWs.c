/* PN_CircularOrbit_GWs.c */
/* Wolfgang Tichy 4/2010 */

#include "sgrid.h"
#include "PN_CircularOrbit_GWs.h"

#define STRLEN 4000

/*
We conmpute h+, hx, psi4 at the points of SphericalDF in sgrid. Then  I can use
spec_2dIntegral with SphericalDF to compute any surface integrals I like.
E.g. for each mode I intend to compute, I could add one radial point.
E.g. make a var where at
i=0     we put rpsi*rY_0^{0}  + ipsi*iY_0^{0}
i=1     we put rpsi*rY_1^{-1} + ipsi*iY_1^{-1}
i=2     we put rpsi*rY_1^{0}  + ipsi*iY_1^{0}
i=3     we put rpsi*rY_1^{1}  + ipsi*iY_1^{1}
i=4     we put rpsi*rY_2^{-2} + ipsi*iY_2^{-2}
i=5     we put rpsi*rY_2^{-1} + ipsi*iY_2^{-1}
i=6     we put rpsi*rY_2^{0}  + ipsi*iY_2^{0}
i=7     we put rpsi*rY_2^{1}  + ipsi*iY_2^{1}
Then I can use spec_2dIntegral over this var to compute the real part
of all modes.
*/                         


/* setup initial boxsizes */
int PN_CircularOrbit_GWs_set_boxsize(tGrid *grid)
{
  int n1, n2, n3;
  int lmax, Lmax;

  printf("PN_CircularOrbit_GWs_set_boxsize: setting box size and coordinates used ...\n");

  /* get lmax from par */
  lmax = Geti("PN_CircularOrbit_GWs_lmax");
  Lmax = Geti("PN_CircularOrbit_GWs_sphere_Lmax");
  if(Lmax<lmax) printf("WARNING: you probably want Lmax>lmax !\n");
  n3 = 2*Lmax + 2; /* really need 2*Lmax + 1, but lets make it an even number */
  n2 = 2*n3;

  /* total number of modes */
  n1 = (lmax + 1)*(lmax + 1);

  /* set some box pars */
  Seti("box0_n1", n1);
  Seti("box0_n2", n2);
  Seti("box0_n3", n3);

  Sets("box0_Coordinates", "SphericalDF");

  Sets("box0_basis1", "Fourier");
  Seti("box0_min1", 1);
  Seti("box0_max1", n1+1);

  Sets("box0_basis2", "Fourier");
  Sets("box0_min2", "0");
  Sets("box0_max2", "2*pi");

  Sets("box0_basis3", "Fourier");
  Sets("box0_min3", "0");
  Sets("box0_max3", "2*pi");

  return 0;
}


/* startup routine */
int PN_CircularOrbit_GWs_startup(tGrid *grid)
{
  enablevar(grid, Ind("PN_CircularOrbit_GWs_hplus")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_hcross")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_Hmode")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_Hmode")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_sYlm")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_sYlm")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_Psi4")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_Psi4")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_Psi4mode")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_Psi4mode")); 
  if(Getv("PN_CircularOrbit_GWs_match_NR", "yes"))
  {
    enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_NRPsi4"));
    enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_NRPsi4"));
    enablevar(grid, Ind("PN_CircularOrbit_GWs_PN_NR_diff"));
  }

  return 0;
}


/* compute Modes */
int PN_CircularOrbit_GWs(tGrid *grid)
{
  FILE *out_orb;
  char orbit[STRLEN];
  FILE *out;
  char outname[STRLEN];
  int lmax, l,m,s, i;
  tBox *box = grid->box[0];
  int n1=box->n1;
  double time;
  int hpind = Ind("PN_CircularOrbit_GWs_hplus");
  int hxind = Ind("PN_CircularOrbit_GWs_hcross");
  int Re_Hmind = Ind("PN_CircularOrbit_GWs_Re_Hmode");
  int Im_Hmind = Ind("PN_CircularOrbit_GWs_Im_Hmode");
  int Re_sYlmind = Ind("PN_CircularOrbit_GWs_Re_sYlm");
  int Im_sYlmind = Ind("PN_CircularOrbit_GWs_Im_sYlm");
  int Re_Psi4ind = Ind("PN_CircularOrbit_GWs_Re_Psi4");
  int Im_Psi4ind = Ind("PN_CircularOrbit_GWs_Im_Psi4");
  int Re_Psi4mind = Ind("PN_CircularOrbit_GWs_Re_Psi4mode");
  int Im_Psi4mind = Ind("PN_CircularOrbit_GWs_Im_Psi4mode");
  double *Re_Hmodep = box->v[Re_Hmind];
  double *Im_Hmodep = box->v[Im_Hmind];
  double *Re_Psi4modep = box->v[Re_Psi4mind];
  double *Im_Psi4modep = box->v[Im_Psi4mind];
  double yvec[12];  /* state vector used in odeint */
  double m1, m2, D; /* masses and distance of observer of h+,hx */
  double chi1x, chi1y, chi1z;  /* x,y,z comp of chi1 = S1/m1^2 */
  double chi2x, chi2y, chi2z;  /* x,y,z comp of chi2 = S2/m2^2 */
  double t1, t2, dt;   /* initial, final time, time step */
  double ti, tf;       /* initial and final time for integrator */
  int ImHmodesign;

  /* flip sign of Im H mode in case of plus_cross format (as for NINJA) */
  if(Getv("PN_CircularOrbit_GWs_HmodeOutputFormat", "plus_cross"))
    ImHmodesign=-1;
  else
    ImHmodesign=+1;

  printf("PN_CircularOrbit_GWs: Computing h+,hx\n");
  s    = Geti("PN_CircularOrbit_GWs_spinweight"); /* spin weight we use */
  lmax = Geti("PN_CircularOrbit_GWs_lmax");
  printf("lmax=%d  s=%d\n", lmax, s);

  /* precompute the sYlm */
  ModeComputer_set_sYlm_inbox(box, Re_sYlmind, Im_sYlmind, s, lmax);

  /* open files */
  sprintf(orbit, "%s/%s", Gets("outdir"), Gets("PN_CircularOrbit_GWs_orbitfile"));
  out_orb = fopen(orbit, "w");
  if(!out_orb) errorexits("failed opening %s", orbit);

  /* write header in output files */
  for(l=2; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* open output file */
    sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
    Gets("PN_CircularOrbit_GWs_hfile_prefix"), l,m,s);
    printf("l=%d m=%+d mode:  %s\n", l,m, outname);
    out = fopen(outname, "w");
    if(!out) errorexits("failed opening %s", outname);

    /* write header */
    if(ImHmodesign==-1)
      fprintf(out, "%s\n", "# time                    hplus              hcross");
    else
      fprintf(out, "%s\n", "# time                    Re_Hmode           Im_Hmode");

    fclose(out); 
  }

  /* initial values */
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
  t1 = Getd("PN_CircularOrbit_GWs_t1");
  t2 = Getd("PN_CircularOrbit_GWs_t2");
  dt = Getd("PN_CircularOrbit_GWs_dt");
  if(t1>t2) dt=-dt;

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
/*
  for(i=0; i<=11; i++)
    printf("yvec[%d] = %g\n", i, yvec[i]);
*/

  /* print header in orbit file */
  fprintf(out_orb, "# time  separation  omega  S1x  S1y  S1z  S2x  S2y  S2z"
                   "  Lnx  Lny  Lnz  Phi\n");

  /* compute h+, hx at different times */
  printf("computing h+, hx at different times:\n");
  for(time=t1; (dt>0)&&(time<=t2) || (dt<0)&&(time>=t2); time+=dt)
  {
    printf("time = %.15g\n", time);

    /* see if we want Psi4 as well */
    if(Getv("PN_CircularOrbit_GWs_computePsi4", "yes"))
    {
      /* compute Psi4 and hplus and hcross */
      compute_FDpsi4_and_hplus_hcross_on_sphere(box, Re_Psi4ind, Im_Psi4ind,
                                                hpind, hxind, yvec, D,m1,m2,
                                                0, dt*0.001, 0,n1-1, 1);
      /* get modes of Psi4 */
      compute_sYlmModes_of_H(box, Re_Psi4ind, Im_Psi4ind,
                             Re_sYlmind, Im_sYlmind, lmax,
                             Re_Psi4mind, Im_Psi4mind, +1); 
      /* output Psi4 modes */
      output_sYlmModes_of_PN_H(Gets("PN_CircularOrbit_GWs_Psi4file_prefix"), 
                               time, Re_Psi4modep, Im_Psi4modep,
                               lmax,s, +1);
    }
    else
    {
      /* compute hplus and hcross */
      compute_hplus_hcross_on_sphere(box, hpind, hxind,
                                     yvec, D,m1,m2, 0,n1-1, 1);
    }

    /* get modes of H = h+ - i hx  <-- sign of Im H is neg */
    compute_sYlmModes_of_H(box, hpind, hxind, Re_sYlmind, Im_sYlmind, lmax,
                           Re_Hmind, Im_Hmind, -1); 
    /* output H modes */
    output_sYlmModes_of_PN_H(Gets("PN_CircularOrbit_GWs_hfile_prefix"),
                             time, Re_Hmodep, Im_Hmodep,
                             lmax,s, ImHmodesign);

    /* output orbit file */
    fprintf(out_orb, "%-.16e", time);
    for(i=0; i<=11; i++)
      fprintf(out_orb, "  %+.16e", yvec[i]);
    fprintf(out_orb, "\n");

    /* advance orbit to time+dt */
    ti=time;
    tf=time+dt;
    xodeint(m1, m2, ti, tf, yvec);
  }
  /* set time to last time so that any other output get correct time label */
  grid->time=time-dt;
  fclose(out_orb);
  return 0;
}


/* compute hplus and hcross on sphere, 
   and put them on radial points imin<=i<=imax */ 
void compute_hplus_hcross_on_sphere(tBox *box, int hpind, int hxind,
                                    double *yvec, 
                                    double D, double m1, double m2,
                                    int imin, int imax,
                                    int set_doublecovered_points)
{
  int i,j,k, ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *hpp = box->v[hpind];
  double *hxp = box->v[hxind];
  /* double *Xp = box->v[Ind("X")]; */
  double *Yp = box->v[Ind("Y")];
  double *Zp = box->v[Ind("Z")];

  /* compute hplus and hcross */
  for(k=0; k<n3; k++)  
  for(j=0; j<n2/2; j++)
  {
    double hplus, hcross;
    double theta,phi;

    ijk=Index(0,j,k);

    /* theta, phi */
    theta = Yp[ijk] + PI/((1+n2%2)*n2);
    phi   = Zp[ijk];

    /* compute hplus and hcross, and put them on radial points imin<=i<=imax */
    compute_hcross_hplus(yvec, &hcross, &hplus, D, theta, phi, m1, m2);
    if(imax>=n1) errorexit("compute_hplus_hcross_on_sphere: imax>=n1");
    for(i=imin; i<=imax; i++)
    {
      hpp[ijk+i] = hplus;
      hxp[ijk+i] = hcross;
    }
  }
  /* set double covered points for h+, hx */
  if(set_doublecovered_points)
  {
    copy_to_doubleCoveredPoints_SphericalDF(box, hpind);
    copy_to_doubleCoveredPoints_SphericalDF(box, hxind);
  }
}

/* use 2nd order fin diff in time to find psi4 from H = h+ - i hx
   and h+, hx */
void compute_FDpsi4_and_hplus_hcross_on_sphere(tBox *box,
                                               int Rpsi4ind, int Ipsi4ind, 
                                               int hpind, int hxind, double *yin, 
                                               double D, double m1, double m2,
                                               double t, double dt,
                                               int imin, int imax,
                                               int set_doublecovered_points)
{
  double yvec[12];
  int i,j,k, ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *hpp = box->v[hpind];
  double *hxp = box->v[hxind];
  double *Rpsi4p = box->v[Rpsi4ind];
  double *Ipsi4p = box->v[Ipsi4ind];

  /* copy yin, (we never change yin) */
  for(i=0; i<12; i++) yvec[i]=yin[i]; 

  /* use 2nd order fin diff in time */
  /* put H(t) = h+ - i hx at i=1 */ 
  compute_hplus_hcross_on_sphere(box, hpind, hxind, yvec, D,m1,m2, 1,1, 0);
  
  /* set yvec at -dt */
  xodeint(m1, m2, t, t-dt, yvec); 
  /* put H(t-dt) = h+ - i hx at i=0 */ 
  compute_hplus_hcross_on_sphere(box, hpind, hxind, yvec, D,m1,m2, 0,0, 0);

  /* set yvec at +dt */
  for(i=0; i<12; i++) yvec[i]=yin[i]; 
  xodeint(m1, m2, t, t+dt, yvec);
  /* put H(t+dt) = h+ - i hx at i=2 */ 
  compute_hplus_hcross_on_sphere(box, hpind, hxind, yvec, D,m1,m2, 2,2, 0);
    
  /* use 2nd order fin diff in time to find psi4
     and set psi4 and h+, hx for all imin<=i<=imax from i=1 */
  for(k=0; k<n3; k++)  
  for(j=0; j<n2/2; j++)
  {
    double hplus, hcross,  hplus_m, hcross_m,  hplus_p, hcross_p; 
    double Repsi4, Impsi4;

    ijk=Index(0,j,k);
    hplus_m  = hpp[ijk+0];
    hcross_m = hxp[ijk+0];
    hplus    = hpp[ijk+1];
    hcross   = hxp[ijk+1];
    hplus_p  = hpp[ijk+2];
    hcross_p = hxp[ijk+2];
    Repsi4 = (hplus_p + hplus_m - 2.0*hplus)/(dt*dt);
    Impsi4 =-(hcross_p + hcross_m - 2.0*hcross)/(dt*dt);

    if(imax>=n1) errorexit("compute_hplus_hcross_on_sphere: imax>=n1");
    for(i=imin; i<=imax; i++)
    {
      hpp[ijk+i] = hplus;
      hxp[ijk+i] = hcross;
      Rpsi4p[ijk+i] = Repsi4;
      Ipsi4p[ijk+i] = Impsi4;
    }
  }

  /* set double covered points for psi4 and h+, hx */
  if(set_doublecovered_points)
  {
    copy_to_doubleCoveredPoints_SphericalDF(box, Rpsi4ind);
    copy_to_doubleCoveredPoints_SphericalDF(box, Ipsi4ind);
    copy_to_doubleCoveredPoints_SphericalDF(box, hpind);
    copy_to_doubleCoveredPoints_SphericalDF(box, hxind);
  }
}

/* compute the modes of any complex function H = ReH + i ImH * ImHsign */
void compute_sYlmModes_of_H(tBox *box, int ReHind, int ImHind,
                            int Re_sYlmind, int Im_sYlmind, int lmax,
                            int Re_Hmind, int Im_Hmind, 
                            int ImHsign)
{
  int l,m, i,j,k, ijk;
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *ReHp = box->v[ReHind];
  double *ImHp = box->v[ImHind];
  double *Re_Hmodep = box->v[Re_Hmind];
  double *Im_Hmodep = box->v[Im_Hmind];
  double *Re_sYlmp = box->v[Re_sYlmind];
  double *Im_sYlmp = box->v[Im_sYlmind];
  double *Xp = box->v[Ind("X")];
  /* double *Yp = box->v[Ind("Y")];
     double *Zp = box->v[Ind("Z")]; */

  /* set integrands */
  i=0;
  for(l=0; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    for(k=0; k<n3; k++)
    for(j=0; j<n2/2; j++)
    {
      double R,I, RsYlm,IsYlm, r; /* ,theta,phi */
      ijk=Index(i,j,k);

      /* get r, theta, phi */
      r     = Xp[ijk];
      /* theta = Yp[ijk] + PI/((1+n2%2)*n2); */
      /* phi   = Zp[ijk]; */

      /* get spin-weighted spherical harmonic sYlm */
      RsYlm = Re_sYlmp[ijk];
      IsYlm = Im_sYlmp[ijk];

      /* get Re and Im part of H */
      R=+ReHp[ijk];
      I=+ImHp[ijk]*ImHsign;

      /* compute modes of H */
      Re_Hmodep[ijk] = RsYlm * R + IsYlm * I;
      Im_Hmodep[ijk] = RsYlm * I - IsYlm * R;

      /* devide by r^2 since spec_sphericalDF2dIntegral multiplies by r^2 */
      Re_Hmodep[ijk] = Re_Hmodep[ijk]/(r*r);
      Im_Hmodep[ijk] = Im_Hmodep[ijk]/(r*r);
    }
    i++;
  }
  
  /* set double covered points of modes */
  copy_to_doubleCoveredPoints_SphericalDF(box, Re_Hmind);
  copy_to_doubleCoveredPoints_SphericalDF(box, Im_Hmind);

  /* integrate over spheres */
  spec_sphericalDF2dIntegral(box, Re_Hmodep, Re_Hmodep);
  spec_sphericalDF2dIntegral(box, Im_Hmodep, Im_Hmodep);
}

/* output modes ReHmode, ImHmode * ImHmodeSign */
void output_sYlmModes_of_PN_H(char *prefix, double time,
                              double *Re_Hmodep, double *Im_Hmodep,
                              int lmax, int s, int ImHmodeSign)
{
  int i, l,m;
  FILE *out;
  char outname[STRLEN];

  /* output */
  i=4;  /* (1+1)^2 */
  for(l=2; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* open output file */
    sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
            prefix, l,m,s);
    out = fopen(outname, "a");
    if(!out) errorexits("failed opening %s", outname);

    /* write time Re and Im part of mode */      
    //printf("%.13g  %.13g  %.13g\n", time1, Re_Hmodep[i], Im_Hmodep[i]);
    fprintf(out, "%-.16e  %+.10e  %+.10e\n",
            time, Re_Hmodep[i], ImHmodeSign*Im_Hmodep[i]);

    fclose(out); 
    i++;
  }
}
