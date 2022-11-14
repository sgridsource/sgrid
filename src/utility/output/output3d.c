/* output3d.c */
/* Wolfgang Tichy 1/2013 */

#include "sgrid.h"
#include "output.h"



/* 3d output of one variable output */
void write3d_boxvar(tBox *box, char *name)
{
  tGrid *grid = box->grid; 
  int text       = Getv("3dformat", "text");
  int binary     = Getv("3dformat", "binary");
  int vtk        = Getv("3dformat", "vtk");
  int fakepoints = Getv("3dformat", "fakepoints");
  int arrange_as_1d = Getv("3dformat", "arrange_as_1d");
  int addpoints  = Getv("3dformat", "addpoints");
  int dump       = Getv("3dformat", "dump");
  int flt        = Getv("3dformat", "float");
  int dbl        = Getv("3dformat", "double");
  FILE *fp;
  int nseries;
  double *pX = box->v[Ind("X")];
  double *pY = box->v[Ind("Y")];
  double *pZ = box->v[Ind("Z")];
  double *px = box->v[Ind("x")];
  double *py = box->v[Ind("y")];
  double *pz = box->v[Ind("z")];
  double *pV = box->v[Ind(name)];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;

  /* return if var has no memory */
  if(box->v[Ind(name)]==NULL) return;  

  /* parameter defaults */
  if (!flt && !dbl) flt = 1;
  if (!text && !binary) binary = 1;

  /* a number that counts, the output */
  nseries = timeforoutput_di_dt(grid, Geti("3doutiter"), Getd("3douttime"));

  /* simply dump output, in sgrid's own (non-standard) format */
  if(dump) dump3d_boxvar(box, name);

  /* VTK */
  /* output one file per time step in separate subdirectories */
  if(vtk)
  {
    if(fakepoints) /* put data on a fake grid with uniform grid spacings dX,dY,dZ */
    {
      double X0 = box->bbox[0];
      double Y0 = box->bbox[2];
      double Z0 = box->bbox[4];
      double X1 = box->bbox[1];
      double Y1 = box->bbox[3];
      double Z1 = box->bbox[5];
      double dX = fabs(X1-X0)/(n1);
      double dY = fabs(Y1-Y0)/(n2);
      double dZ = fabs(Z1-Z0)/(n3);

      if(arrange_as_1d) /* pretend that all points are along X-dir */
      {
        n1 = n1*n2*n3;
        n2 = n3 = 1;
        X0 = Y0 = Z0 = 0.0;
        dX = dY = dZ = 1.0;
      }

      /* open file (returns non-null file pointer) */
      fp = fopen_vtk(name, "XYZ", box->b, nseries-1);

      /* write header */
      fprintf(fp, "# vtk DataFile Version 2.0\n");
      fprintf(fp, "variable %s, box %d, time %.16g, "
              "Note: uniform grid spacings below are fake\n", 
              name, box->b, grid->time);
      fprintf(fp, binary ? "BINARY" : "ASCII\n");
      fprintf(fp, "\n");
      fprintf(fp, "DATASET STRUCTURED_POINTS\n");
      fprintf(fp, "DIMENSIONS %d %d %d\n", n1, n2, n3);
      fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", X0, Y0, Z0);
      fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", dX, dY, dZ);
      fprintf(fp, "\n");
      fprintf(fp, "POINT_DATA %d\n", n1*n2*n3);
      fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
      fprintf(fp, "LOOKUP_TABLE default\n");
      /* write data,
         has to be in the file right after the \n of the last header line */
      write_raw_vtk_data(fp, pV, n1*n2*n3,1,0, dbl, flt, text, binary);

      fclose(fp);
    }
    else
    {
      if(addpoints) /* add grid in x,y,z coords. */
      {
        /* open file (returns non-null file pointer) */
        fp = fopen_vtk(name, "xyz", box->b, nseries-1);
        /* here we use the lower case Cartesian x,y,z */

        fprintf(fp, "# vtk DataFile Version 2.0\n");
        fprintf(fp, "variable %s, box %d, time %.16g\n",
                name, box->b, grid->time);
        fprintf(fp, binary ? "BINARY" : "ASCII\n");
        fprintf(fp, "\n");
        fprintf(fp, "DATASET STRUCTURED_GRID\n");
        fprintf(fp, "DIMENSIONS %d %d %d\n", n1, n2, n3);
        fprintf(fp, "POINTS %d %s\n", n1*n2*n3, dbl ? "double" : "float");
        write_raw_vtk_points(fp, px, py, pz, n1*n2*n3,1,0, dbl, flt, text, binary);

        fprintf(fp, "\n\n");
        fprintf(fp, "POINT_DATA %d\n", n1*n2*n3);
        fprintf(fp, "SCALARS %s %s\n", name, dbl ? "double" : "float");
        fprintf(fp, "LOOKUP_TABLE default\n");
        /* write data,
           has to be in the file right after the \n of the last header line */
        write_raw_vtk_data(fp, pV, n1*n2*n3,1,0, dbl, flt, text, binary);

        fclose(fp);
      }
      else /* use RECTILINEAR_GRID */
      {
        /* open file (returns non-null file pointer) */
        fp = fopen_vtk(name, "XYZ", box->b, nseries-1);

        /* write header */
        fprintf(fp, "# vtk DataFile Version 2.0\n");
        fprintf(fp, "variable %s, box %d, time %.16g\n", 
                name, box->b, grid->time);
        fprintf(fp, binary ? "BINARY" : "ASCII\n");
        fprintf(fp, "\n");
        fprintf(fp, "DATASET RECTILINEAR_GRID\n");
        fprintf(fp, "DIMENSIONS %d %d %d\n", n1, n2, n3);
        fprintf(fp, "X_COORDINATES %d %s\n", n1, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pX, n1,1,0, 1, 0, text, binary);
        fprintf(fp, "\n\n");
        fprintf(fp, "Y_COORDINATES %d %s\n", n2, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pY, n2,n1,0, 1, 0, text, binary);
        fprintf(fp, "\n\n");
        fprintf(fp, "Z_COORDINATES %d %s\n", n3, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pZ, n3,n1*n2,0, 1, 0, text, binary);
        fprintf(fp, "\n\n");
        fprintf(fp, "POINT_DATA %d\n", n1*n2*n3);
        fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
        fprintf(fp, "LOOKUP_TABLE default\n");
        /* write data,
           has to be in the file right after the \n of the last header line */
        write_raw_vtk_data(fp, pV, n1*n2*n3,1,0, dbl, flt, text, binary);

        fclose(fp);
      }
    }
  } /* end if(vtk) */
}
