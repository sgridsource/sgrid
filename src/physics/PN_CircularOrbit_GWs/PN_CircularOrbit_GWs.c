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
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_var")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_var")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_mode")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_mode")); 
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Re_sYlm"));
  enablevar(grid, Ind("PN_CircularOrbit_GWs_Im_sYlm"));

  return 0;
}


/* compute Modes */
int PN_CircularOrbit_GWs(tGrid *grid)
{
//  FILE *in1;
//  char Re_file[STRLEN];
  FILE *out;
  char outname[STRLEN];
  int lmax, l,m,s, i,j,k, ijk;
  tBox *box = grid->box[0];
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double time;
  int hpind = Ind("PN_CircularOrbit_GWs_hplus");
  int hxind = Ind("PN_CircularOrbit_GWs_hcross");
  int Re_Hmind = Ind("PN_CircularOrbit_GWs_Re_Hmode");
  int Im_Hmind = Ind("PN_CircularOrbit_GWs_Im_Hmode");
  int Re_sYlmind = Ind("PN_CircularOrbit_GWs_Re_sYlm");
  int Im_sYlmind = Ind("PN_CircularOrbit_GWs_Im_sYlm");
  double *hpp = box->v[hpind];
  double *hxp = box->v[hxind];
  double *Re_Hmodep = box->v[Re_Hmind];
  double *Im_Hmodep = box->v[Im_Hmind];
  double *Re_sYlmp = box->v[Re_sYlmind];
  double *Im_sYlmp = box->v[Im_sYlmind];
  double *Xp = box->v[Ind("X")];
  double *Yp = box->v[Ind("Y")];
  double *Zp = box->v[Ind("Z")];
  double yvec[12];  /* state vector used in odeint */
  double m1, m2, D; /* masses and distance of observer of h+,hx */
  double hi1, hi2;  /* is this chi1/2 ??? Is it used at all by anything??? */
  double ti, tf, dt;   /* initial, final time, time step */
  double t1, t2;       /* initial and final time for integrator */

  printf("PN_CircularOrbit_GWs: Computing h+,hx\n");
  s    = Geti("PN_CircularOrbit_GWs_spinweight"); /* spin weight we use */
  lmax = Geti("PN_CircularOrbit_GWs_lmax");
  printf("lmax=%d  s=%d\n", lmax, s);

  /* precompute the sYlm */
  ModeComputer_set_sYlm_inbox(box, Re_sYlmind, Im_sYlmind, s, lmax);

  /* open files */
  //sprintf(Re_file, "%s", Gets("PN_CircularOrbit_GWs_Re_sphere_data"));
  //in1 = fopen(Re_file, "rb");
  //if(!in1) errorexits("failed opening %s", Re_file);

  /* write header in output files */
  for(l=0; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* open output file */
    sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
    Gets("PN_CircularOrbit_GWs_outfile_prefix"), l,m,s);
    printf("l=%d m=%+d mode:  %s\n", l,m, outname);
    out = fopen(outname, "w");
    if(!out) errorexits("failed opening %s", outname);

    /* write header */
    fprintf(out, "%s\n", "# time             Re_Hmode           Im_Hmode");

    fclose(out); 
  }

  /* initial values */
  m1 = 0.5;
  m2 = 0.5;
  hi1 = 1.0;
  hi2 = 1.0;
  ti = 0.0;
  tf = 10.0;
  D = 1.0;
  yvec[0] = 0.0;  // you don't need to initialize this field it will be used as the storage of the separation
  yvec[1] = 0.01; // initial orbital frequency
  yvec[2] = 0.1;  // S1x not Scap1! S1 = hi1*m1^2*Scap1
  yvec[3] = 0.2;  // S1y
  yvec[4] = 0.3;  // S1z
  yvec[5] =-0.1;  // S2x
  yvec[6] = 0.2;  // S2y
  yvec[7] = 0.3;  // S2z 
  yvec[8] = 0.0;  // Lnx
  yvec[9] = 0.0;  // Lny
  yvec[10]= 1.0;  // Lnz
  yvec[11]= 0.0;  // Phi orbital phase         

  /* compute h+, hx at different times */
  printf("computing h+, hx at different times:\n");
  for(time=ti; time<=tf; time+=dt)
  {
    printf("time = %g\n", time);

    /* compute orbit at time */
    t1=time-dt;
    t2=time;
    xodeint(m1, m2, hi1, hi2, t1, t2, yvec);

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

      /* compute hplus and hcross, and put them all radial points */
      /* FIXME: maybe store only at i=0 */
      compute_hcross_hplus(yvec, &hcross, &hplus, D, theta, phi, m1, m2);
      for(i=0; i<n1; i++)
      {
        hpp[ijk+i] = hplus;
        hxp[ijk+i] = hcross;
      }
    }
    /* set double covered points for h+, hx */
    copy_to_doubleCoveredPoints_SphericalDF(box, hpind);
    copy_to_doubleCoveredPoints_SphericalDF(box, hpind);

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

        /* get Re and Im part of H = h+ - i hx */
        R=+hpp[ijk];
        I=-hxp[ijk];
        
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
            
    /* output */
    i=0;
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      /* open output file */
      sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
      Gets("PN_CircularOrbit_GWs_outfile_prefix"), l,m,s);
      out = fopen(outname, "a");
      if(!out) errorexits("failed opening %s", outname);

      /* write time Re and Im part of mode */      
      //printf("%.13g  %.13g  %.13g\n", time1, Re_Hmodep[i], Im_Hmodep[i]);
      fprintf(out, "%-15.8g  %16.10e  %16.10e\n", time, Re_Hmodep[i], Im_Hmodep[i]);

      fclose(out); 
      i++;
    }
  }

//  fclose(in1);
  return 0;
}
