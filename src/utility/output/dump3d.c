/* dump_out3d.c */
/* Wolfgang Tichy, Sept 2008 */

#include "sgrid.h"
#include "output.h"



void dump3d_write_box(tBox *box, FILE *fp, int iv)
{
  double *pv = box->v[iv];
  int ijk;

  if(pv==NULL) return;

  fprintf(fp, "# \"time = %.16g\" , box%d: n1=%d n2=%d n3=%d , "
          "N=%d numbers , type=float\n",
          box->grid->time, box->b,box->n1,box->n2,box->n3, box->nnodes);
  
  forallpoints(box,ijk)
  {
    float fl;
    fl=pv[ijk];
    fwrite(&fl, sizeof(fl), 1, fp);
  }
  fprintf(fp,"\n\n");
}


/* 3d dump of box */
void dump3d_boxvar(tBox *box, char *name)
{
  FILE *fp;
  char filename[1000];

  if(box->v[Ind(name)]==NULL) return;

  snprintf(filename, 999, "%s/%s.XYZ%d", Gets("outdir"), name, box->b);

  /* open file */
  fp = fopen(filename, "a");
  if (!fp) errorexits("failed opening %s", filename);

  dump3d_write_box(box, fp, Ind(name));
    
  /* close file */
  fclose(fp);
}
