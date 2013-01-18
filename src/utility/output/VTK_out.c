/* VTK_out.c */
/* Wolfgang Tichy 1/2013 */

#include "sgrid.h"
#include "output.h"



/* About binary formats:

   The VTK lagacy data format used here is by default big endian binary data.
   This used to be the default on "workstations", while x86 "PCs"
   typically use little endian.

   There exists /usr/include/endian.h, but how standard is this?
   also see /usr/include/byteswap.h, which does not do 8 bytes for non-gcc?
   see include/bits/byteswap.h for gcc/x86 optimized code
*/
#include <endian.h>
#if  __BYTE_ORDER == __LITTLE_ENDIAN
#define BYTE_ORDER_LITTLE 1
#else
#define BYTE_ORDER_LITTLE 0
#endif


/* in 3d and 2d it usually comes down to this: write numbers 
   efficiency shouldn't matter much since usually the disk is still slower
   defaults are based on legacy:
   - default is double if float is not specified
   - text/binary has no default (error exit)
*/
void write_raw_vtk_data(FILE *fp, double *buffer, int n, int stride, int offset,
		        int dbl, int flt, int text, int binary)
{
  int swap = BYTE_ORDER_LITTLE;  // if little endian, then swap since
                                 // the vtk default is big endian
  float xfloat;
  double xdouble;
  int i;

  if(!text && !binary)
    errorexit("write_raw_vtk_data: pick text or binary format");
  if(text && binary)
    errorexit("write_raw_vtk_data: pick either text or binary format");
  if(dbl && flt)
    errorexit("write_raw_vtk_data: pick either double or float format");

  if(text)
  {
    if(dbl) 
      for (i = 0; i < n; i++)
	fprintf(fp, "%.16g\n", buffer[stride*i+offset]);
    else
      for (i = 0; i < n; i++)
	fprintf(fp, "%.7g\n", (float) buffer[stride*i+offset]);
  }

  if(binary)
  {
    if(sizeof(char) != 1)
      errorexit("write_raw_vtk_data: size of char is not 1");

    if(!swap)
    {
      if(dbl)
      {
	double xdouble;
	for(i = 0; i < n; i++)
	{
	  xdouble = buffer[stride*i+offset];
	  fwrite(&xdouble, sizeof(double), 1, fp);
	}
      }
      if(flt)
      {
	float xfloat;
	for (i = 0; i < n; i++)
	{
	  xfloat = buffer[stride*i+offset];
	  fwrite(&xfloat, sizeof(float), 1, fp);
	}
      }
    }

    if(swap)
    {
      if(dbl && sizeof(double) == 8)
      {
	double xdouble;
	char c[8], *x;
	for (i = 0; i < n; i++)
	{
	  xdouble = buffer[stride*i+offset];
	  x = (char *) &xdouble;
	  c[0] = x[7];
	  c[1] = x[6];
	  c[2] = x[5];
	  c[3] = x[4];
	  c[4] = x[3];
	  c[5] = x[2];
	  c[6] = x[1];
	  c[7] = x[0];
	  fwrite(c, sizeof(char), 8, fp);
	}
      }
      else if (flt && sizeof(float) == 4)
      {
	float xfloat;
	char c[8], *x;
	for (i = 0; i < n; i++)
	{
	  xfloat = buffer[stride*i+offset];
	  x = (char *) &xfloat;
	  c[0] = x[3];
	  c[1] = x[2];
	  c[2] = x[1];
	  c[3] = x[0];
	  fwrite(c, sizeof(char), 4, fp);
	}
      }
      else
	errorexit("write_raw_vtk_data: size of float/double is not 4/8");
    }
  }
}

/* open file for vtk writing */
FILE *fopen_vtk(char *varname, char *suffix, int b, int n)
{
  char filename[1000];
  char *outdir = Gets("outdir");
  FILE *fp;

  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_vtk", outdir, varname, suffix, b);
  fp = fopen(filename, "r");
  if(!fp)
    mkdir(filename, 0777);  
  else
    fclose(fp);
  
  /* open file */
  snprintf(filename, 1000, "%s/%s.%s%d_vtk/%s.%s%d_%04d.vtk", 
	   outdir, varname, suffix, b,
	   varname, suffix, b, n);
  fp = fopen(filename, "wb");
  if(!fp) 
    errorexits("failed opening %s", filename);
  
  /* return non-null file pointer */
  return fp;
}
