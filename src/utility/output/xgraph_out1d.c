/* xgraph_out1d.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"


void write_line(tBox *box, FILE *fp, int line, int iv)
{
  int ix = Ind("X");
  int iy = Ind("Y");
  int iz = Ind("Z");
  double *px;
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
  
  if(line==1)
  {
    px = box->v[ix];
    jmin = jmax = n2/2;
    kmin = kmax = n3/2;
  }
  else if(line==2)
  {
    px = box->v[iy];
    imin = imax = n1/2;
    kmin = kmax = n3/2;
  }
  else
  {
    px = box->v[iz];
    jmin = jmax = n2/2;
    imin = imax = n1/2;
  }

  if(pv==NULL) return;

  fprintf(fp, "\"time = %.16g\"\n", box->grid->time);
              
  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
    for(j=jmin; j<=jmax; j++)
      for(i=imin; i<=imax; i++)
      {
        index=Index(i,j,k);
        //printf("i=%d j=%d k=%d index=%d\n", i,j,k, index);
        fprintf(fp, "%.16g %.16g\n", px[index], pv[index]);
      }
}


/* 1d xgraph plots */
void xgraph_out1_boxvar(tBox *box, char *name)
{
  FILE *fx, *fy, *fz;
  char xfilename[1000];
  char yfilename[1000];
  char zfilename[1000];
  
  snprintf(xfilename, 999, "%s/%s.X%d", Gets("outdir"), name, box->b);
  snprintf(yfilename, 999, "%s/%s.Y%d", Gets("outdir"), name, box->b);
  snprintf(zfilename, 999, "%s/%s.Z%d", Gets("outdir"), name, box->b);

  /* open file */
  fx = fopen(xfilename, "a");
  if (!fx) errorexits("failed opening %s", xfilename);
  fy = fopen(yfilename, "a");
  if (!fy) errorexits("failed opening %s", yfilename);
  fz = fopen(zfilename, "a");
  if (!fz) errorexits("failed opening %s", zfilename);

  /* lines */
  write_line(box, fx, 1 , Ind(name));  /* x-line */
  write_line(box, fy, 2 , Ind(name));  /* y-line */
  write_line(box, fz, 3 , Ind(name));  /* z-line */

  /* close file */
  fprintf(fx,"\n"); fclose(fx);
  fprintf(fy,"\n"); fclose(fy);
  fprintf(fz,"\n"); fclose(fz);
}

