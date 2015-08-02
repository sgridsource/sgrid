/* endianIO.c */
/* Wolfgang Tichy, August 2015 */

#include "sgrid.h"

/* About binary formats:
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


/* use fwrite to write an array of doubles to a file in little endian format */
size_t fwrite_double_little(const double *buf, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */ 

  /* just use fwrite if we are little endian */
  if(little)
    return fwrite(buf, sizeof(double), nmemb, fp);

  /* sanity check */
  if(sizeof(char) != 1)
    errorexit("fwrite_double_little: size of char is not 1");

  /* if we get here we assume big endian */
  if(sizeof(double) == 8)
  {
    double xdouble;
    char c[8], *x;
    size_t i, count=0;

    for(i=0; i<nmemb; i++)
    {
      xdouble = buf[i];
      x = (char *) &xdouble;
      c[0] = x[7];
      c[1] = x[6];
      c[2] = x[5];
      c[3] = x[4];
      c[4] = x[3];
      c[5] = x[2];
      c[6] = x[1];
      c[7] = x[0];
      count = count + fwrite(c, sizeof(char), 8, fp);
    }
    return count/8;
  }
  errorexit("fwrite_double_little: size of double is not 8");
  return -1; /* hopefully we never get here */
}


/* use fread to read an array of doubles from a file in little endian format */
size_t fread_double_little(double *buf, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */ 

  /* just use fread if we are little endian */
  if(little)
    return fread(buf, sizeof(double), nmemb, fp);

  /* sanity check */
  if(sizeof(char) != 1)
    errorexit("fread_double_little: size of char is not 1");

  /* if we get here we assume big endian */
  if(sizeof(double) == 8)
  {
    double xdouble;
    char c[8], *x;
    size_t i, count=0;

    for(i=0; i<nmemb; i++)
    {
      count = count + fread(c, sizeof(char), 8, fp);
      x = (char *) &xdouble;
      x[0] = c[7];
      x[1] = c[6];
      x[2] = c[5];
      x[3] = c[4];
      x[4] = c[3];
      x[5] = c[2];
      x[6] = c[1];
      x[7] = c[0];
      buf[i] = xdouble;
    }
    return count/8;
  }
  errorexit("fread_double_little: size of double is not 8");
  return -1; /* hopefully we never get here */
}
