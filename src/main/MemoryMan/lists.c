/* lists.c */
/* Wolfgang Tichy, May 2017 */

#include <stdlib.h>
#include <stdio.h>

/* compile list functions with lists of type int
   e.g. intList *alloc_intList(void);  */
#define TYP int
#include "list_templates.c"
#undef TYP
/* to use them in a file we need to add
#define TYP double
#include "list_templates.h"
#undef TYP
to this file. */


/* if we need the same lists but with double entries, do this:
#define TYP double
#include "list_templates.c"
#undef TYP
*/
