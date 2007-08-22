/* xgraph_out1d.c */
/* Wolfgang Tichy, April 2005 */

#include "sgrid.h"
#include "output.h"


void write_line(tBox *box, FILE *fp, int line, int iv)
{
  int ix;
  int iy;
  int iz;
  double *px;
  double *pv = box->v[iv];
  int n1=box->n1;
  int n2=box->n2;
  int n3=box->n3;
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;
  int index;
  char str[1000];
  char Xvarname[1000];
  char Yvarname[1000];
  char Zvarname[1000];

  if(pv==NULL) return;

  fprintf(fp, "\"time = %.16g\" ", box->grid->time);

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
  
  if(line==1)
  {
    px = box->v[ix];

    snprintf(str, 999, "outputY0_box%d", box->b);
    jmin = jmax = find_ind_closest_to_Y0(box,Getd(str));

    snprintf(str, 999, "outputZ0_box%d", box->b);
    kmin = kmax = find_ind_closest_to_Z0(box,Getd(str));

    fprintf(fp, ", j=%d, k=%d, Y=%g, Z=%g\n", jmin, kmin, 
            box->v[Ind("Y")][Index(0,jmin,kmin)],
            box->v[Ind("Z")][Index(0,jmin,kmin)]);
  }
  else if(line==2)
  {
    px = box->v[iy];

    snprintf(str, 999, "outputX0_box%d", box->b);
    imin = imax = find_ind_closest_to_X0(box,Getd(str));

    snprintf(str, 999, "outputZ0_box%d", box->b);
    kmin = kmax = find_ind_closest_to_Z0(box,Getd(str));

    fprintf(fp, ", i=%d, k=%d, X=%g, Z=%g\n", imin, kmin, 
            box->v[Ind("X")][Index(imin,0,kmin)],
            box->v[Ind("Z")][Index(imin,0,kmin)]);
  }
  else
  {
    px = box->v[iz];

    snprintf(str, 999, "outputX0_box%d", box->b);
    imin = imax = find_ind_closest_to_X0(box,Getd(str));

    snprintf(str, 999, "outputY0_box%d", box->b);
    jmin = jmax = find_ind_closest_to_Y0(box,Getd(str));

    fprintf(fp, ", i=%d, j=%d, X=%g, Y=%g\n", imin, jmin, 
            box->v[Ind("X")][Index(imin,jmin,0)],
            box->v[Ind("Y")][Index(imin,jmin,0)]);
  }

              
  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
    for(j=jmin; j<=jmax; j++)
      for(i=imin; i<=imax; i++)
      {
        index=Index(i,j,k);
        //printf("i=%d j=%d k=%d index=%d\n", i,j,k, index);
        fprintf(fp, "%.16g %.16g\n", px[index], pv[index]);
      }
  fprintf(fp,"\n");
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

  /* close files */
  fclose(fx);
  fclose(fy);
  fclose(fz);
}

