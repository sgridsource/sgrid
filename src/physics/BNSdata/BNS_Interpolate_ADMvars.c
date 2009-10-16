/* BNS_Interpolate_ADMvars.c */
/* Wolfgang Tichy 10/2009 */
/* interpolate onto points which we read from a file */

#include "sgrid.h"
#include "BNSdata.h"



/* position filepointer after the string str */
int position_fileptr_after_str(FILE *in, char *str)
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
int BNS_Interpolate_ADMvars(tGrid *grid)
{
  FILE *in, *out;
  char *pointsfile;
  char *outfile;
  double x,y,z, val;
  double X,Y,Z;
  int ind,j,b;
  tVarList *vlu;
  tVarList *vlc;

  if(GetsLax("BNSdata_Interpolate_pointsfile")==0) return 0;
  printf("BNS_Interpolate_ADMvars:\n");

  /* allocate varlists */
  vlu = vlalloc(grid);
  vlc = vlalloc(grid);

  /* add all vars to vlu */
  vlpush(vlu, Ind("gxx"));
  vlpush(vlu, Ind("alpha"));
  vlpush(vlu, Ind("betax"));
  vlpush(vlu, Ind("Kxx"));
  vlpush(vlu, Ind("rho"));
  vlpush(vlu, Ind("jx"));
  vlpush(vlu, Ind("Sxx"));

  /* now duplicate vlu to get vlc */  
  vlc = AddDuplicateEnable(vlu, "_c");

  /* write coeffs of vlu in all boxes into vlc */
  forallboxes(grid, b)
    spec_Coeffs_varlist(grid->box[b], vlu, vlc);

  /* filenames */
  pointsfile = Gets("BNSdata_Interpolate_pointsfile");
  outfile    = Gets("BNSdata_Interpolate_output");
  printf("BNSdata_Interpolate_pointsfile=%s\n", pointsfile);
  printf("BNSdata_Interpolate_output=%s\n", outfile);

  /* open both files */
  in = fopen(pointsfile, "rb");
  if(!in) errorexits("failed opening %s", pointsfile);
  out = fopen(outfile, "a");
  if(!out) errorexits("failed opening %s", outfile);

  /* write header info */  
  fprintf(out, "%s", "#");
  for(j=0; j<vlu->n; j++)
    fprintf(out, " %s", VarName(vlu->index[j]));
  fprintf(out, "%s", "\n");
  fprintf(out, "%s\n", "$BEGIN_data:");
  
  position_fileptr_after_str(in, "$BEGIN_data:\n");
  while(read_next_xyz_from_pointsfile(in, &x,&y,&z)!=EOF)
  {
    /* initial guess for X,Y,Z, b: */
    nearest_b_XYZ_of_xyz(grid, &b, &ind, &X,&Y,&Z, x,y,z);

    /* get X,Y,Z, b of x,y,z */
    b=BNSgrid_Get_BoxAndCoords_of_xyz(grid, &X,&Y,&Z,b, x,y,z);

    /* interpolate vlu (using coeffs in vlc) to X,Y,Z in box b */
    for(j=0; j<vlc->n; j++)
    {
      tBox *box = grid->box[b];
      double *c = box->v[vlc->index[j]];
      val = spec_interpolate(box, c, X,Y,Z);
      fwrite(&val, sizeof(double), 1, out);
    }
  }

  /* close files */
  fclose(out);
  fclose(in);
  
  /* free var lists */
  vlfree(vlu);
  VLDisableFree(vlc);
      
  return 0;
}
