/* gnuplot_out2d.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"


void write_plane(tBox *box, FILE *fp, int normal, int plane, int iv)
{
  int ix = Ind("X");
  int iy = Ind("Y");
  int iz = Ind("Z");
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
  
  imin = jmin = kmin=0;
  imax = n1 - 1;
  jmax = n2 - 1;
  kmax = n3 - 1;
  
  if(normal==1)
  {
    p1 = box->v[iy];
    p2 = box->v[iz];
    imin = imax = plane;
  }
  else if(normal==2)
  {
    p1 = box->v[ix];
    p2 = box->v[iz];
    jmin = jmax = plane;
  }
  else
  {
    p1 = box->v[ix];
    p2 = box->v[iy];
    kmin = kmax = plane;
  }

  fprintf(fp, "# time = %.16g\"\n", box->grid->time);
              
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
}


/* 2d gnuplot */
void gnuplot_out2d_boxvar(tBox *box, char *name)
{
  FILE *fx, *fy, *fz;
  char xfilename[1000];
  char yfilename[1000];
  char zfilename[1000];
  
  snprintf(xfilename, 999, "%s/%s.XY%d", Gets("outdir"), name, box->b);
  snprintf(yfilename, 999, "%s/%s.XZ%d", Gets("outdir"), name, box->b);
  snprintf(zfilename, 999, "%s/%s.YZ%d", Gets("outdir"), name, box->b);

  /* open file */
  fx = fopen(xfilename, "a");
  if (!fx) errorexits("failed opening %s", xfilename);
  fy = fopen(yfilename, "a");
  if (!fy) errorexits("failed opening %s", yfilename);
  fz = fopen(zfilename, "a");
  if (!fz) errorexits("failed opening %s", zfilename);

  /* xy-plane:  z = 0 */
  write_plane(box, fx, 3 , box->n3/2, Ind(name));
  write_plane(box, fy, 2 , box->n2/2, Ind(name));
  write_plane(box, fz, 1 , box->n1/2, Ind(name));

  /* close file */
  fprintf(fx,"\n"); fclose(fx);
  fprintf(fy,"\n"); fclose(fy);
  fprintf(fz,"\n"); fclose(fz);
}
