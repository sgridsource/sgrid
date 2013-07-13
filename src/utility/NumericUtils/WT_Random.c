/* WT_Random.c contains functions that return random numbers */
/* (c) Wolfgang Tichy 13.7.2013 */


#include <stdlib.h>

/* On Linux both rand and random return a random integer
   between 0 and RAND_MAX, and both use the same random number generator */
/* BUT: rand conforms to C89, C99.  random only to 4.3BSD, POSIX.1-2001. */
/* YET on some non-Linux systems lower-order bits of rand are much 
   less random  than  the  higher-order bits. random does not suffer 
   from this problem. */
/* for now I use rand */
#define RAND_FUNC rand

/* return random double in [0,1] */
double RND(void)
{
  double r = RAND_FUNC();
  return r/RAND_MAX;
}
