/* ModeComputer.c */
/* Wolfgang Tichy 2/2010 */

#include "sgrid.h"
#include "ModeComputer.h"

#define STRLEN 4000

/*
If we have output at the points of SphericalDF in sgrid I can use
spec_2dIntegral with SphericalDF to compute any surface integrals I like.
E.g. for each mode I intend to compute, I could add one radial point.
E.g. make a var where at
i=0     we put rpsi
i=1     we put rpsi*rY_1^{-1} + ipsi*iY_1^{-1}
i=2     we put rpsi*rY_1^{0}  + ipsi*iY_1^{0}
i=3     we put rpsi*rY_1^{1}  + ipsi*iY_1^{1}
i=4     we put rpsi*rY_2^{-2} + ipsi*iY_2^{-2}
i=5     we put rpsi*rY_2^{-1} + ipsi*iY_2^{-1}
i=6     we put rpsi*rY_2^{0}  + ipsi*iY_2^{0}
i=6     we put rpsi*rY_2^{1}  + ipsi*iY_2^{1}
This way I can use spec_2dIntegral over this var to compute all modes.
*/                         


/* setup initial boxsizes */
int ModeComputer_set_boxsize(tGrid *grid)
{
  FILE *in1;
  char Re_file[STRLEN];
  char gridtype[STRLEN];
  int n1, n2, n3, ntheta, nphi;
  int lmax;

  printf("ModeComputer_set_boxsize: setting box size and coordinates used ...\n");
    
  /* open files */
  sprintf(Re_file, "%s", Gets("ModeComputer_Re_sphere_data"));
  in1 = fopen(Re_file, "rb");
  if(!in1) errorexits("failed opening %s", Re_file);

  /* get grid type */
  ModeComputer_read_spheregrid(in1, &n2, &n3, gridtype, &ntheta, &nphi);

  /* compute lmax from number of Fourier modes in theta dir */
  /* lmax = (n2-1)/2; */
  /* get lmax from par */
  lmax = Geti("ModeComputer_lmax");
  
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

  fclose(in1);
  return 0;
}


/* startup routine */
int ModeComputer_startup(tGrid *grid)
{
  enablevar(grid, Ind("ModeComputer_Re_var")); 
  enablevar(grid, Ind("ModeComputer_Im_var")); 
  enablevar(grid, Ind("ModeComputer_Re_mode")); 
  enablevar(grid, Ind("ModeComputer_Im_mode")); 

  return 0;
}


/* read n2, n3 and gridtype from file pointer fp */
int ModeComputer_read_spheregrid(FILE *fp, int *n2, int *n3, char *gridtype,
                                 int *ntheta, int *nphi)
{
  char str[STRLEN];

  /* check for label SphericalDF: */
  if(fgotonext(fp, "SphericalDF:")==EOF)
    errorexit("ModeComputer_read_spheregrid: could not find SphericalDF:");

  /* get n2, n3 and gridtype */
  fgetparameter(fp, "n2", str);  *n2=atoi(str);
  fgetparameter(fp, "n3", str);  *n3=atoi(str);
  fgetparameter(fp, "gridtype", gridtype);
  fgetparameter(fp, "ntheta", str);  *ntheta=atoi(str);
  fgetparameter(fp, "nphi", str);    *nphi=atoi(str);
  
  return 0;
}


/* read data for time t from file pointer fp into var with index vind */
int ModeComputer_read_spheredata(tGrid *grid, FILE *fp, char *gridtype,
                                 int ntheta, int nphi,
                                 double *t, int vind)
{
  tBox *box = grid->box[0];
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double *var = box->v[vind];
  int i,j,k;
  char str[STRLEN];
  double val;
  int bitantsym, rotantsym, octantsym;

  /* check if there is one more iteration */
  if(fgotonext(fp, "iteration")==EOF) return EOF;

  /* get time t */
  fgetparameter(fp, "time", str);  *t=atof(str);
  
  /* read data into var array */
  for(k=0; k<nphi; k++)
  for(j=0; j<ntheta; j++)
  {
    fscanf(fp, "%s", str);
    val = atof(str);
    for(i=0; i<n1; i++) var[Index(i,j,k)]=val; 
  }

  /* deal with symmetries in bam */
  if(strstr(VarName(vind), "_Re_")!=NULL)
  {
    bitantsym = Geti("ModeComputer_Re_bitantsym");
    rotantsym = Geti("ModeComputer_Re_rotantsym");
    octantsym = Geti("ModeComputer_Re_octantsym");
  }
  else if(strstr(VarName(vind), "_Im_")!=NULL)
  {
    bitantsym = Geti("ModeComputer_Im_bitantsym");
    rotantsym = Geti("ModeComputer_Im_rotantsym");
    octantsym = Geti("ModeComputer_Im_octantsym");
  }
  if(strcmp(gridtype,"bitant")==0 ||
     strcmp(gridtype,"quadrant")==0 ||
     strcmp(gridtype,"octant")==0)
  { /* copy into region with theta>PI/2 i.e. j>=ntheta */
    for(k=0; k<nphi; k++)
    for(j=ntheta; j<n2/2; j++)
    for(i=0; i<n1; i++)
      var[Index(i,j,k)] = bitantsym * var[Index(i,n2/2-j-1,k)];
  }
  if(strcmp(gridtype,"octant")==0)
  { /* copy into region with PI/2<=phi<PI i.e. nphi<=k<n3/2 */
    for(k=nphi; k<n3/2; k++)
    for(j=0; j<n2/2; j++)
    for(i=0; i<n1; i++)
      var[Index(i,j,k)] = octantsym * var[Index(i,j,n3/2-k)];
  }
  if(strcmp(gridtype,"rotant")==0 ||
     strcmp(gridtype,"quadrant")==0 ||
     strcmp(gridtype,"octant")==0)
  { /* copy into region with phi>=PI i.e. k>=n3/2 */
    for(k=nphi; k<n3; k++)
    for(j=0; j<n2/2; j++)
    for(i=0; i<n1; i++)
      var[Index(i,j,k)] = rotantsym * var[Index(i,j,k-n3/2)];
  }

  /* copy var into double covered regions */
  copy_to_doubleCoveredPoints_SphericalDF(box, vind);
  
  return 0;
}


/* compute Modes */
int ModeComputer(tGrid *grid)
{
  FILE *in1;
  FILE *in2;
  FILE *out;
  char Re_file[STRLEN];
  char Im_file[STRLEN];
  char outname[STRLEN];
  char gridtype[STRLEN];
  int ntheta, nphi;
  int lmax, l,m,s, i,j,k, ijk;
  tBox *box = grid->box[0];
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  double time1, time2;
  int Re_vind = Ind("ModeComputer_Re_var");
  int Im_vind = Ind("ModeComputer_Im_var");
  int Re_mind = Ind("ModeComputer_Re_mode");
  int Im_mind = Ind("ModeComputer_Im_mode");
  double *Re_varp = box->v[Re_vind];
  double *Im_varp = box->v[Im_vind];
  double *Re_modep = box->v[Re_mind];
  double *Im_modep = box->v[Im_mind];
  double *Xp = box->v[Ind("X")];
  double *Yp = box->v[Ind("Y")];
  double *Zp = box->v[Ind("Z")];
  
  printf("ModeComputer: Computing Modes\n");
  s    = Geti("ModeComputer_spinweight"); /* spin weight we use */
  lmax = Geti("ModeComputer_lmax");

  /* open files */
  sprintf(Re_file, "%s", Gets("ModeComputer_Re_sphere_data"));
  sprintf(Im_file, "%s", Gets("ModeComputer_Im_sphere_data"));
  in1 = fopen(Re_file, "rb");
  if(!in1) errorexits("failed opening %s", Re_file);
  in2 = fopen(Im_file, "rb");
  if(!in2) errorexits("failed opening %s", Im_file);

  /* write header in output files */
  for(l=0; l<=lmax; l++)
  for(m=-l; m<=l; m++)
  {
    /* open output file */
    sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
    Gets("ModeComputer_outfile_prefix"), l,m,s);
    printf("making mode file %s\n", outname);
    out = fopen(outname, "w");
    if(!out) errorexits("failed opening %s", outname);

    /* write header */
    fprintf(out, "%s\n", "# time             Re_mode            Im_mode");

    fclose(out); 
  }

  /* get grid type */
  ModeComputer_read_spheregrid(in1, &j, &k, gridtype, &ntheta, &nphi);

  /* step from time label to time label and read data */
  while(ModeComputer_read_spheredata(grid, in1, gridtype, ntheta, nphi, 
                                     &time1, Re_vind)!=EOF)
  {
    ModeComputer_read_spheredata(grid, in2, gridtype, ntheta, nphi, 
                                 &time2, Im_vind);
    if(time1!=time2) errorexit("times in files for Re and Im part disagree");
    printf("time = %g\n", time1);

    i=0;
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      for(k=0; k<n3; k++)
      for(j=0; j<n2/2; j++)
      {
        double R,I, r,theta,phi, costheta, RsYlm,IsYlm;
        ijk=Index(i,j,k);

        /* get r, theta, phi */
        r     = Xp[ijk];
        theta = Yp[ijk] + PI/((1+n2%2)*n2);
        phi   = Zp[ijk];
        costheta = cos(theta);

        /* get spin-weighted spherical harmonic sYlm */
        RsYlm = Re_sYlm(l,m,s, costheta, phi);
        IsYlm = Im_sYlm(l,m,s, costheta, phi);
        
        /* get Re and Im part of data */
        R=Re_varp[ijk];
        I=Im_varp[ijk];
        
        /* There is a choice of sign here: define the inner product by
           (f,g) = int f^* g
           and define
           psi_mode = (Y, psi)  */
        Re_modep[ijk] = RsYlm * R + IsYlm * I;
        Im_modep[ijk] = RsYlm * I - IsYlm * R;

        /* devide by r^2 since spec_sphericalDF2dIntegral multiplies by r^2 */
        Re_modep[ijk] = Re_modep[ijk]/(r*r);
        Im_modep[ijk] = Im_modep[ijk]/(r*r);
      }
      i++;
    }
    
    /* set double covered points of modes */
    copy_to_doubleCoveredPoints_SphericalDF(box, Re_mind);
    copy_to_doubleCoveredPoints_SphericalDF(box, Im_mind);

    /* integrate over spheres */
    spec_sphericalDF2dIntegral(box, Re_modep, Re_modep);
    spec_sphericalDF2dIntegral(box, Im_modep, Im_modep);
            
    /* output */
    i=0;
    for(l=0; l<=lmax; l++)
    for(m=-l; m<=l; m++)
    {
      /* open output file */
      sprintf(outname, "%s/%sl%dm%d_s%d.t", Gets("outdir"),
      Gets("ModeComputer_outfile_prefix"), l,m,s);
      out = fopen(outname, "a");
      if(!out) errorexits("failed opening %s", outname);

      /* write time Re and Im part of mode */      
      //printf("%.13g  %.13g  %.13g\n", time1, Re_modep[i], Im_modep[i]);
      fprintf(out, "%-15.8g  %16.10e  %16.10e\n", time1, Re_modep[i], Im_modep[i]);

      fclose(out); 
      i++;
    }
  }

  fclose(in1);
  fclose(in2);
  return 0;
}
