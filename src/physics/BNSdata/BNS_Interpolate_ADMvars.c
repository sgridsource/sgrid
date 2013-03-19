/* BNS_Interpolate_ADMvars.c */
/* Wolfgang Tichy 10/2009 */
/* interpolate onto points which we read from a file */

#include "sgrid.h"
#include "BNSdata.h"



/* position filepointer below the label str */
int position_fileptr_below_label(FILE *in, char *str)
{
  char line[1000];
  
  while(fgets(line, 999, in)!=NULL)
  {
    if(strstr(line, str)!=NULL) return 1; //break;
  }
  return EOF;
}

/* read xyz of next point */
int read_next_xyz_from_pointsfile(FILE *in, double *x, double *y, double *z)
{
  double xyz[3];
  int n;

  n=fread(xyz, sizeof(double), 3, in);
  if(n!=3) n=EOF;
  *x=xyz[0];
  *y=xyz[1];
  *z=xyz[2];
  return n;
}


/* interpolate ADM initial data onto points listed in pointsfile */
/* Note: here we check extra pars that are not introduced in sgrid_BNSdata.c :
   BNSdata_Interpolate_pointsfile = some_file_name
   BNSdata_Interpolate_verbose = no # yes
   BNSdata_Interpolate_make_finer_grid2_forXYZguess = yes # no
   BNSdata_Interpolate_max_xyz_diff = 1e-8
   The last 3 are just for debugging. Their values are set e.g. by bam. */
int BNS_Interpolate_ADMvars(tGrid *grid)
{
  int pr=GetvLax("BNSdata_Interpolate_verbose", "yes");
  FILE *in, *out;
  char *pointsfile;
  char *outfile;
  tGrid *grid2=NULL;
  double x,y,z, val;
  double gxx=-1e300;
  double X,Y,Z;
  int ind,j,b;
  tVarList *vlu;
  tVarList *vlc;
  int corot1 = Getv("BNSdata_rotationstate1","corotation");
  int corot2 = Getv("BNSdata_rotationstate2","corotation");

  if(GetsLax("BNSdata_Interpolate_pointsfile")==0) return 0;
  prdivider(0);
  printf("BNS_Interpolate_ADMvars:\n");

  /* allocate varlists */
  vlu = vlalloc(grid);
  vlc = vlalloc(grid);

  /* add all vars to vlu */
  vlpush(vlu, Ind("alpha"));
  vlpush(vlu, Ind("BNSdata_Bx"));
  vlpush(vlu, Ind("gxx"));
  vlpush(vlu, Ind("Kxx"));
  vlpush(vlu, Ind("BNSdata_q"));
  vlpush(vlu, Ind("BNSdata_VRx"));
  /* Those we do not need:
    vlpush(vlu, Ind("BNSdata_Sigmax"));
    vlpush(vlu, Ind("BNSdata_wBx"));
    vlpush(vlu, Ind("rho"));
    vlpush(vlu, Ind("jx"));
    vlpush(vlu, Ind("Sxx"));
  */

  /* now duplicate vlu to get vlc */  
  vlc = AddDuplicateEnable(vlu, "_c");

  /* write coeffs of vlu in all boxes into vlc */
  forallboxes(grid, b)
    spec_Coeffs_varlist(grid->box[b], vlu, vlc);

  /* pick grid2 used to get guess for X,Y,Z */
  if(GetvLax("BNSdata_Interpolate_make_finer_grid2_forXYZguess", "no"))
    grid2 = grid;
  else
  {
    /* make a finer grid2 so that nearest_b_XYZ_of_xyz_inboxlist used
       with grid2 finds points that are closer to the correct point. */
    printf(" making finer grid2 to get a good guess for X,Y,Z ...\n");
    fflush(stdout);
    /* use 40 points in A,B and leave point number in other directions. */
    grid2 = make_grid_with_sigma_pm(grid, 40, 
                                    grid->box[1]->n3, grid->box[5]->n1);
  }

  /* filenames */
  pointsfile = Gets("BNSdata_Interpolate_pointsfile");
  outfile    = Gets("BNSdata_Interpolate_output");
  printf(" BNSdata_Interpolate_pointsfile = %s\n", pointsfile);
  printf(" BNSdata_Interpolate_output = %s\n", outfile);
  fflush(stdout);

  /* open both files */
  in = fopen(pointsfile, "rb");
  if(!in) errorexits("failed opening %s", pointsfile);
  out = fopen(outfile, "wb");
  if(!out) errorexits("failed opening %s", outfile);


  /* write header info */  
  fprintf(out, "%s", "#");
  for(j=0; j<vlu->n; j++)
    fprintf(out, " %s", VarName(vlu->index[j]));
  fprintf(out, "%s", "\n");
  fprintf(out, "%s\n", "$BEGIN_data:");
  
  j=position_fileptr_below_label(in, "$BEGIN_data:");
  if(j==EOF) errorexits("could not find $BEGIN_data: in %s", pointsfile);
  while(read_next_xyz_from_pointsfile(in, &x,&y,&z)!=EOF)
  {
    int blist[6];

    if(pr) printf("(x,y,z)=(%g,%g,%g)\n", x,y,z);

    /* initial guess for X,Y,Z, b: */
    if(x>=0.0) { blist[0]=1;  blist[1]=0;   blist[2]=5; }
    else       { blist[0]=2;  blist[1]=3;   blist[2]=4; }
    nearest_b_XYZ_of_xyz_inboxlist(grid2, blist,3, &b, &ind, &X,&Y,&Z, x,y,z);
    if(b<4)
    {
      if(dequal(Y, 0.0)) Y=0.01;
      if(dequal(Y, 1.0)) Y=1.0-0.01;
      Z = Arg_plus(y,z);
    }
    else
    {
      X=x;
      Y=y;
      Z=z;
    }
    if(pr) printf("guess:  b=%d (X,Y,Z)=(%g,%g,%g)  nearest ind=%d\n", b, X,Y,Z, ind);

    /* get X,Y,Z, b of x,y,z */
    b=BNSgrid_Get_BoxAndCoords_of_xyz(grid, &X,&Y,&Z,b, x,y,z);
    if(pr) printf("actual: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z);
    if(b<0)
    {
      printf("point: (x,y,z)=(%g,%g,%g)\n", x,y,z);
      printf("error: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z); 
      errorexit("could not find point");
    }
    if( GetsLax("BNSdata_Interpolate_max_xyz_diff")!=0 )
    {
      double tol= Getd("BNSdata_Interpolate_max_xyz_diff");
      tBox *box = grid->box[b];
      double xs,ys,zs, diff;

      if(box->x_of_X[1]==NULL) { xs=X;  ys=Y;  zs=Z; }
      else
      {
        xs = box->x_of_X[1]((void *) box, -1, X,Y,Z);
        ys = box->x_of_X[2]((void *) box, -1, X,Y,Z);
        zs = box->x_of_X[3]((void *) box, -1, X,Y,Z);
      }
      diff = sqrt( (x-xs)*(x-xs) + (y-ys)*(y-ys) + (z-zs)*(z-zs) );
      if(diff>tol)
      {
        printf("error: b=%d {X,Y,Z}={%.15g,%.15g,%.15g}\n", b, X,Y,Z);
        printf("=> {xs,ys,zs}={x(X,Y,Z), y(X,Y,Z), z(X,Y,Z)} with:\n");
        printf(" {xs,ys,zs}={%.15g,%.15g,%.15g}\n", xs,ys,zs);
        printf("BUT {x,y,z}={%.15g,%.15g,%.15g}\n", x,y,z);
        printf(" diff=%.15g\n", diff);
        fflush(stdout);
        errorexit("x(X,Y,Z), y(X,Y,Z), z(X,Y,Z) are wrong!");
      }
    }

    /* interpolate vlu (using coeffs in vlc) to X,Y,Z in box b */
    for(j=0; j<vlc->n; j++)
    {
      tBox *box = grid->box[b];
      double *c = box->v[vlc->index[j]];

      /* HACK: don't interpolate for some vars that are zero
         or already known */
      if( strcmp(VarName(vlu->index[j]),"gxx")==0 ) 
        gxx = val = spec_interpolate(box, c, X,Y,Z);
      else if( strcmp(VarName(vlu->index[j]),"gyy")==0 ||
               strcmp(VarName(vlu->index[j]),"gzz")==0   )
        val = gxx;
      else if( strcmp(VarName(vlu->index[j]),"gxy")==0 ||
               strcmp(VarName(vlu->index[j]),"gxz")==0 ||
               strcmp(VarName(vlu->index[j]),"gyz")==0   )
        val=0.0;
      else if( strcmp(VarName(vlu->index[j]),"BNSdata_q")==0 )
      {
        if(b==1 || b==2) val=0.0;
        else             val=spec_interpolate(box, c, X,Y,Z);
      }
      else if( strcmp(VarName(vlu->index[j]),"BNSdata_VRx")==0 ||
               strcmp(VarName(vlu->index[j]),"BNSdata_VRy")==0 ||
               strcmp(VarName(vlu->index[j]),"BNSdata_VRz")==0   )
      {
        if(b==1 || b==2) val=0.0;
        else if((b==0 || b==5) && corot1 || (b==3 || b==4) && corot2) val=0.0;
        else val=spec_interpolate(box, c, X,Y,Z);
      }
      else val = spec_interpolate(box, c, X,Y,Z);
      /* if we always interpolate we need:
      val = spec_interpolate(box, c, X,Y,Z); */
      if(!finite(val))
      {
        printf("point:  (x,y,z)=(%g,%g,%g)\n", x,y,z);
        printf("NAN at: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z); 
        errorexit("spec_interpolate returned NAN, probably (X,Y,Z) was bad!");
      }
      if(pr) printf("%s=%g\n", VarName(vlu->index[j]), val);
      fwrite(&val, sizeof(double), 1, out);
    }
  }
  
  /* remove grid2, if it does not point to grid itself */
  if(grid2 != grid) free_grid(grid2);

  /* close files */
  fclose(out);
  fclose(in);
  
  /* free var lists */
  vlfree(vlu);
  VLDisableFree(vlc);

  printf("BNS_Interpolate_ADMvars: finished interpolations "
         "and wrote all output.\n");
  prdivider(0);
  fflush(stdout);
     
  return 0;
}
