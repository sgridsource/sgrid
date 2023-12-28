/* endianIO.c */
/* Wolfgang Tichy, Feb. 2019 */

/* FIXME: According to Rob Pike the code below is stupid.
   From https://commandcenter.blogspot.com/2012/04/byte-order-fallacy.html :
   Let's say your data stream has a little-endian-encoded 32-bit integer.
   Here's how to extract it (assuming unsigned bytes):
   i = (data[0]<<0) | (data[1]<<8) | (data[2]<<16) | (data[3]<<24);
   If it's big-endian, here's how to extract it:
   i = (data[3]<<0) | (data[2]<<8) | (data[1]<<16) | (data[0]<<24);
   On the other hand what I have is general and works! */


#include "sgrid.h"


/* We can check the byte order with this function.
   Without inline and already with -O1 this will be turned into:
        mov    $0x1,%eax
        retq
   But it will also be inlined, and thus will cause no overhead at all! */
static inline unsigned int SGRID_byte_order_is_little(void)
{
  unsigned int ui = 1;
  unsigned char *s = (unsigned char *) &ui;
  return s[0]; /* returns first byte in ui, which is 1 for little-endian */
}


/* About binary formats:
   There is <endian.h> or "/usr/include/endian.h" on Linux. This is not
   standard yet. On BSD this can be in <sys/endian.h>, or on crazy systems
   like macos in <machine/endian.h>.
   On Linux see also "/usr/include/byteswap.h", and
   include/bits/byteswap.h for gcc/x86 optimized code */

/* Directives to check check if BYTE_ORDER_LITTLE was set in MyConfig.
   If not, use info from <endian.h> to set it: */
/*
#ifndef BYTE_ORDER_LITTLE
#include <endian.h>
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define BYTE_ORDER_LITTLE 1
#else
#define BYTE_ORDER_LITTLE 0
#endif
#endif
*/

/* instead of the "defines" above for BYTE_ORDER_LITTLE we simply use the
   inline-function SGRID_byte_order_is_little to just test endianness */
#define BYTE_ORDER_LITTLE SGRID_byte_order_is_little()
/* this should cause no overhead! */




/* use fwrite to write an array, but swap byte order */
size_t SGRID_fwrite_swapbytes(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  const char *buf = ptr;
  size_t i, count;
  char c[size];
  unsigned int b;

  for(count=0, i=0; i<nmemb; i++)
  {
    for(b=0; b<size; b++) c[b] = buf[i*size + size-1-b];
    count = count + fwrite(c, sizeof(char), size, fp);
  }
  return count/size;
}

/* use fread to write an array, but swap byte order */
size_t SGRID_fread_swapbytes(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  char *buf = ptr;
  size_t i, count;
  char c[size];
  unsigned int b;

  for(count=0, i=0; i<nmemb; i++)
  {
    count += fread(c, sizeof(char), size, fp);
    for(b=0; b<size; b++) buf[i*size + b] = c[size-1-b];
  }
  return count/size;
}

/* use fwrite to write to a file in little endian format */
size_t SGRID_fwrite_little(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fwrite if we are little endian */
  if(little)
    return fwrite(ptr, size, nmemb, fp);
  else
    return SGRID_fwrite_swapbytes(ptr, size, nmemb, fp);
}


/* use fread to read from a file in little endian format */
size_t SGRID_fread_little(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are little endian */
  if(little)
    return fread(ptr, size, nmemb, fp);
  else
    return SGRID_fread_swapbytes(ptr, size, nmemb, fp);
}

/* use fwrite to write to a file in big endian format */
size_t SGRID_fwrite_big(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fwrite if we are big endian */
  if(!little)
    return fwrite(ptr, size, nmemb, fp);
  else
    return SGRID_fwrite_swapbytes(ptr, size, nmemb, fp);
}


/* use fread to read from a file in big endian format */
size_t SGRID_fread_big(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int little = BYTE_ORDER_LITTLE; /* endianess */

  /* just use fread if we are big endian */
  if(!little)
    return fread(ptr, size, nmemb, fp);
  else
    return SGRID_fread_swapbytes(ptr, size, nmemb, fp);
}


/* return value of BYTE_ORDER_LITTLE */
int SGRID_return_BYTE_ORDER_LITTLE(void)
{
  return BYTE_ORDER_LITTLE; /* endianess */
}

/* print what byte order we have */
int SGRID_print_endian_info(tGrid *grid)
{
  printf("BYTE_ORDER_LITTLE = %d\n", SGRID_return_BYTE_ORDER_LITTLE());
  return 0;
}


/******************************************************/
/* special funcs just for doubles, used in checkpoint */
/******************************************************/

/* write an array of doubles to a file in little endian format */
size_t fwrite_double_little(const double *buf, size_t nmemb, FILE *fp)
{
  /* sanity checks, that we probably do not need */
  if(sizeof(char) != 1)
    errorexit("fwrite_double_little: size of char is not 1");
  if(sizeof(double) != 8)
    errorexit("fwrite_double_little: size of double is not 8");

  return SGRID_fwrite_little(buf, sizeof(double), nmemb, fp);

}

/* read an array of doubles from a file in little endian format */
size_t fread_double_little(double *buf, size_t nmemb, FILE *fp)
{
  /* sanity checks, that we probably do not need */
  if(sizeof(char) != 1)
    errorexit("fread_double_little: size of char is not 1");
  if(sizeof(double) != 8)
    errorexit("fread_double_little: size of double is not 8");

  return SGRID_fread_little(buf, sizeof(double), nmemb, fp);
}
