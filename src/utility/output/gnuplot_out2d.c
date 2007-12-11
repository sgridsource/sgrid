/* gnuplot_out2d.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "output.h"



void write_plane(tBox *box, FILE *fp, int normal, int plane, int iv)
{
  int ix;
  int iy;
  int iz;
  double *p1;
  double *p2;
  double *pv = box->v[iv];
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;
  int index;
  char Xvarname[1000];
  char Yvarname[1000];
  char Zvarname[1000];

  if(pv==NULL) return;

  fprintf(fp, "# \"time = %.16g\"", box->grid->time);
  
  /* get pointers to X, Y, Z */
  snprintf(Xvarname, 999, "outputReplaceXby_box%d", box->b);
  snprintf(Yvarname, 999, "outputReplaceYby_box%d", box->b);
  snprintf(Zvarname, 999, "outputReplaceZby_box%d", box->b);
  ix = Ind(Gets(Xvarname));
  iy = Ind(Gets(Yvarname));
  iz = Ind(Gets(Zvarname));
  
  imin = jmin = kmin=0;
  imax = n1 - 1;
  jmax = n2 - 1;
  kmax = n3 - 1;
  
  if(normal==1)
  {
    p1 = box->v[iy];
    p2 = box->v[iz];
    imin = imax = plane;
    fprintf(fp, ", i=%d, X=%.16g\n", plane, box->v[Ind("X")][Index(plane,0,0)]);
  }
  else if(normal==2)
  {
    p1 = box->v[ix];
    p2 = box->v[iz];
    jmin = jmax = plane;
    fprintf(fp, ", j=%d, Y=%.16g\n", plane, box->v[Ind("Y")][Index(0,plane,0)]);
  }
  else
  {
    p1 = box->v[ix];
    p2 = box->v[iy];
    kmin = kmax = plane;
    fprintf(fp, ", k=%d, Z=%.16g\n", plane, box->v[Ind("Z")][Index(0,0,plane)]);
  }

  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
  {
    for(j=jmin; j<=jmax; j++)
    {
      for(i=imin; i<=imax; i++)
      {
        index=Index(i,j,k);
        fprintf(fp, "%.16g %.16g %.16g\n", p1[index], p2[index], pv[index]);
      }
      if(normal==3) fprintf(fp, "\n");
    }
    if(normal==2 || normal==1) fprintf(fp, "\n");
  }
  fprintf(fp,"\n");
}


/* 2d gnuplot */
void gnuplot_out2d_boxvar(tBox *box, char *name)
{
  FILE *fXY, *fXZ, *fYZ;
  char XYfilename[1000];
  char XZfilename[1000];
  char YZfilename[1000];
  char str[1000];

  if(box->v[Ind(name)]==NULL) return;

  snprintf(XYfilename, 999, "%s/%s.XY%d", Gets("outdir"), name, box->b);
  snprintf(XZfilename, 999, "%s/%s.XZ%d", Gets("outdir"), name, box->b);
  snprintf(YZfilename, 999, "%s/%s.YZ%d", Gets("outdir"), name, box->b);

  /* open file */
  fXY = fopen(XYfilename, "a");
  if (!fXY) errorexits("failed opening %s", XYfilename);
  fXZ = fopen(XZfilename, "a");
  if (!fXZ) errorexits("failed opening %s", XZfilename);
  fYZ = fopen(YZfilename, "a");
  if (!fYZ) errorexits("failed opening %s", YZfilename);

  /* XY-plane:  Z = Z0 */
  snprintf(str, 999, "outputZ0_box%d", box->b);
  write_plane(box, fXY, 3 , find_ind_closest_to_Z0(box,Getd(str)), Ind(name));

  /* XZ-plane:  Y = Y0 */
  snprintf(str, 999, "outputY0_box%d", box->b);
  write_plane(box, fXZ, 2 , find_ind_closest_to_Y0(box,Getd(str)), Ind(name));

  /* YZ-plane:  X = X0 */
  snprintf(str, 999, "outputX0_box%d", box->b);
  write_plane(box, fYZ, 1 , find_ind_closest_to_X0(box,Getd(str)), Ind(name));
    
  /* close files */
  fclose(fXY);
  fclose(fXZ);
  fclose(fYZ);
}
